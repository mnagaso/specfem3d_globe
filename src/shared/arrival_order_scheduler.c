/*=====================================================================
 *
 *  Barrier-free arrival-order I/O reservation scheduler (single-node).
 *
 *=====================================================================*/

#define _GNU_SOURCE
#include "arrival_order_scheduler.h"

#include <errno.h>
#include <fcntl.h>
#include <pthread.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

typedef struct arrival_order_state {
  unsigned int magic;
  int version;
  pthread_mutex_t mutex;
  int initialized;

  int current_checkpoint_id;
  int num_groups;
  int next_sequence;
  int seen[AOS_MAX_GROUPS];
  long long last_reserved_start_ns;
  int last_reserved_valid;

  int owner_pid;
  char run_key[AOS_MAX_RUN_KEY];
} arrival_order_state_t;

static arrival_order_state_t *g_state = NULL;
static int g_shm_fd = -1;
static size_t g_shm_size = 0;
static char g_shm_path[512];
static char g_last_error[512];
static int g_num_groups = 0;
static char g_run_key[AOS_MAX_RUN_KEY];

static void aos_set_error(const char *msg) {
  snprintf(g_last_error, sizeof(g_last_error), "%s", msg ? msg : "unknown");
}

const char *aos_last_error(void) {
  return g_last_error;
}

long long aos_now_ns(void) {
  struct timespec ts;
  if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0) {
    return -1;
  }
  return (long long)ts.tv_sec * 1000000000LL + (long long)ts.tv_nsec;
}

void aos_compute_reservation(
    long long arrival_ns,
    long long last_reserved_start_ns,
    int is_first,
    double spacing_sec,
    long long *scheduled_start_ns,
    double *assigned_delay_sec) {

  long long spacing_ns;
  long long scheduled;

  if (spacing_sec < 0.0) {
    spacing_sec = 0.0;
  }
  spacing_ns = (long long)(spacing_sec * 1.0e9 + 0.5);

  if (is_first) {
    scheduled = arrival_ns;
  } else {
    long long earliest = last_reserved_start_ns + spacing_ns;
    scheduled = arrival_ns > earliest ? arrival_ns : earliest;
  }

  if (scheduled_start_ns) {
    *scheduled_start_ns = scheduled;
  }
  if (assigned_delay_sec) {
    long long delay_ns = scheduled - arrival_ns;
    if (delay_ns < 0) {
      delay_ns = 0;
    }
    *assigned_delay_sec = (double)delay_ns * 1.0e-9;
  }
}

static void aos_build_path(char *out, size_t out_len,
                           const char *state_dir, const char *run_key) {
  const char *dir = (state_dir && state_dir[0]) ? state_dir : "/dev/shm";
  snprintf(out, out_len, "%s/specfem_aos_%s", dir, run_key);
}

static int aos_owner_alive(int owner_pid) {
  if (owner_pid <= 1) {
    return 0;
  }
  if (kill(owner_pid, 0) == 0) {
    return 1;
  }
  /* EPERM means the process exists but we cannot signal it. */
  return (errno == EPERM) ? 1 : 0;
}

static int aos_state_ready(const arrival_order_state_t *st,
                           const char *run_key,
                           int num_groups) {
  if (!st) {
    return 0;
  }
  if (!st->initialized || st->magic != AOS_MAGIC || st->version != AOS_VERSION) {
    return 0;
  }
  if (st->num_groups != num_groups) {
    return 0;
  }
  if (strncmp(st->run_key, run_key, AOS_MAX_RUN_KEY) != 0) {
    return 0;
  }
  if (!aos_owner_alive(st->owner_pid)) {
    return 0;
  }
  return 1;
}

static void aos_unmap_local(void) {
  if (g_state != NULL) {
    munmap(g_state, g_shm_size);
    g_state = NULL;
  }
  if (g_shm_fd >= 0) {
    close(g_shm_fd);
    g_shm_fd = -1;
  }
}

static int aos_init_mutex(pthread_mutex_t *mutex) {
  pthread_mutexattr_t attr;
  int rc;

  rc = pthread_mutexattr_init(&attr);
  if (rc != 0) {
    return rc;
  }
  rc = pthread_mutexattr_setpshared(&attr, PTHREAD_PROCESS_SHARED);
  if (rc != 0) {
    pthread_mutexattr_destroy(&attr);
    return rc;
  }
  rc = pthread_mutexattr_setrobust(&attr, PTHREAD_MUTEX_ROBUST);
  if (rc != 0) {
    pthread_mutexattr_destroy(&attr);
    return rc;
  }
  rc = pthread_mutex_init(mutex, &attr);
  pthread_mutexattr_destroy(&attr);
  return rc;
}

static int aos_lock(void) {
  int rc = pthread_mutex_lock(&g_state->mutex);
  if (rc == EOWNERDEAD) {
    /* Previous owner died while holding the lock. State may be corrupt;
       mark as needing reset and recover the mutex. */
    rc = pthread_mutex_consistent(&g_state->mutex);
    if (rc != 0) {
      aos_set_error("pthread_mutex_consistent failed after EOWNERDEAD");
      return AOS_ERR_MUTEX;
    }
    g_state->initialized = 0;
    g_state->current_checkpoint_id = -1;
    g_state->next_sequence = 0;
    g_state->last_reserved_valid = 0;
    g_state->last_reserved_start_ns = 0;
    memset(g_state->seen, 0, sizeof(g_state->seen));
    return AOS_OK;
  }
  if (rc != 0) {
    snprintf(g_last_error, sizeof(g_last_error),
             "pthread_mutex_lock failed: %d", rc);
    return AOS_ERR_MUTEX;
  }
  return AOS_OK;
}

static int aos_unlock(void) {
  int rc = pthread_mutex_unlock(&g_state->mutex);
  if (rc != 0) {
    snprintf(g_last_error, sizeof(g_last_error),
             "pthread_mutex_unlock failed: %d", rc);
    return AOS_ERR_MUTEX;
  }
  return AOS_OK;
}

static void aos_reset_checkpoint(arrival_order_state_t *st, int checkpoint_id) {
  int i;
  st->current_checkpoint_id = checkpoint_id;
  st->next_sequence = 0;
  st->last_reserved_valid = 0;
  st->last_reserved_start_ns = 0;
  for (i = 0; i < AOS_MAX_GROUPS; ++i) {
    st->seen[i] = 0;
  }
}

int aos_init(
    const char *run_key,
    const char *state_dir,
    int num_groups,
    int world_rank,
    int create_new) {

  int fd;
  arrival_order_state_t *mapped;

  (void)world_rank;

  if (!run_key || !run_key[0] || num_groups <= 0 || num_groups > AOS_MAX_GROUPS) {
    aos_set_error("aos_init: invalid arguments");
    return AOS_ERR_ARGS;
  }

  if (g_state != NULL) {
    aos_set_error("aos_init: already initialized");
    return AOS_ERR_ARGS;
  }

  snprintf(g_run_key, sizeof(g_run_key), "%s", run_key);
  g_num_groups = num_groups;
  aos_build_path(g_shm_path, sizeof(g_shm_path), state_dir, run_key);
  g_shm_size = sizeof(arrival_order_state_t);

  if (create_new) {
    /* Replace any leftover object from a previous crashed job. */
    unlink(g_shm_path);
    fd = open(g_shm_path, O_RDWR | O_CREAT | O_EXCL, 0600);
    if (fd < 0) {
      snprintf(g_last_error, sizeof(g_last_error),
               "aos_init: create shm failed: %s", strerror(errno));
      return AOS_ERR_SHM;
    }

    if (ftruncate(fd, (off_t)g_shm_size) != 0) {
      snprintf(g_last_error, sizeof(g_last_error),
               "aos_init: ftruncate failed: %s", strerror(errno));
      close(fd);
      unlink(g_shm_path);
      return AOS_ERR_SHM;
    }

    mapped = (arrival_order_state_t *)mmap(NULL, g_shm_size,
                                           PROT_READ | PROT_WRITE,
                                           MAP_SHARED, fd, 0);
    if (mapped == MAP_FAILED) {
      snprintf(g_last_error, sizeof(g_last_error),
               "aos_init: mmap failed: %s", strerror(errno));
      close(fd);
      unlink(g_shm_path);
      return AOS_ERR_SHM;
    }

    g_shm_fd = fd;
    g_state = mapped;
    memset(g_state, 0, sizeof(*g_state));
    g_state->magic = AOS_MAGIC;
    g_state->version = AOS_VERSION;
    g_state->num_groups = num_groups;
    g_state->owner_pid = (int)getpid();
    snprintf(g_state->run_key, sizeof(g_state->run_key), "%s", run_key);
    g_state->current_checkpoint_id = -1;
    if (aos_init_mutex(&g_state->mutex) != 0) {
      aos_set_error("aos_init: mutex init failed");
      aos_unmap_local();
      unlink(g_shm_path);
      return AOS_ERR_MUTEX;
    }
    /* Publish last so attachers never observe a half-initialized object. */
    __sync_synchronize();
    g_state->initialized = 1;
    aos_set_error("");
    return AOS_OK;
  }

  /*
   * Non-creator attach protocol:
   * A leftover shm from a crashed job keeps initialized=1. If we open that
   * inode once and keep it, we race with the creator's unlink+recreate and
   * can reserve against stale seen[] bits (AOS_ERR_DOUBLE_RESERVE).
   * Always re-open until the mapped object is ready AND its owner is alive.
   */
  {
    int spins;
    int saw_bad = 0;
    for (spins = 0; spins < 5000; ++spins) {
      fd = open(g_shm_path, O_RDWR, 0600);
      if (fd < 0) {
        usleep(1000);
        continue;
      }
      mapped = (arrival_order_state_t *)mmap(NULL, g_shm_size,
                                             PROT_READ | PROT_WRITE,
                                             MAP_SHARED, fd, 0);
      if (mapped == MAP_FAILED) {
        close(fd);
        usleep(1000);
        continue;
      }
      g_shm_fd = fd;
      g_state = mapped;

      if (aos_state_ready(g_state, run_key, num_groups)) {
        aos_set_error("");
        return AOS_OK;
      }

      saw_bad = 1;
      aos_unmap_local();
      usleep(1000);
    }

    if (saw_bad) {
      aos_set_error("aos_init: stale or corrupt shared state detected");
      return AOS_ERR_CORRUPT;
    }
    snprintf(g_last_error, sizeof(g_last_error),
             "aos_init: timed out waiting for shm %s", g_shm_path);
    return AOS_ERR_SHM;
  }
}

int aos_reserve(
    int checkpoint_id,
    int group_id,
    double spacing_sec,
    long long *arrival_ns,
    long long *scheduled_start_ns,
    int *arrival_sequence,
    double *assigned_delay_sec) {

  long long now;
  long long scheduled = 0;
  double delay = 0.0;
  int seq = 0;
  int is_first;
  int rc;

  if (g_state == NULL) {
    aos_set_error("aos_reserve: not initialized");
    return AOS_ERR_NOT_INIT;
  }
  if (group_id < 0 || group_id >= g_num_groups) {
    aos_set_error("aos_reserve: invalid group_id");
    return AOS_ERR_ARGS;
  }

  now = aos_now_ns();
  if (now < 0) {
    aos_set_error("aos_reserve: clock_gettime failed");
    return AOS_ERR_ARGS;
  }

  rc = aos_lock();
  if (rc != AOS_OK) {
    return rc;
  }

  if (g_state->magic != AOS_MAGIC || g_state->version != AOS_VERSION) {
    aos_unlock();
    aos_set_error("aos_reserve: corrupt shared state");
    return AOS_ERR_CORRUPT;
  }

  if (!g_state->initialized ||
      g_state->current_checkpoint_id != checkpoint_id) {
    aos_reset_checkpoint(g_state, checkpoint_id);
    g_state->initialized = 1;
  }

  if (g_state->seen[group_id]) {
    aos_unlock();
    snprintf(g_last_error, sizeof(g_last_error),
             "aos_reserve: double reserve checkpoint=%d group=%d",
             checkpoint_id, group_id);
    return AOS_ERR_DOUBLE_RESERVE;
  }

  is_first = !g_state->last_reserved_valid;
  aos_compute_reservation(
      now,
      g_state->last_reserved_start_ns,
      is_first,
      spacing_sec,
      &scheduled,
      &delay);

  seq = g_state->next_sequence;
  g_state->next_sequence += 1;
  g_state->seen[group_id] = 1;
  g_state->last_reserved_start_ns = scheduled;
  g_state->last_reserved_valid = 1;

  rc = aos_unlock();
  if (rc != AOS_OK) {
    return rc;
  }

  if (arrival_ns) {
    *arrival_ns = now;
  }
  if (scheduled_start_ns) {
    *scheduled_start_ns = scheduled;
  }
  if (arrival_sequence) {
    *arrival_sequence = seq;
  }
  if (assigned_delay_sec) {
    *assigned_delay_sec = delay;
  }

  aos_set_error("");
  return AOS_OK;
}

int aos_sleep_until(long long scheduled_start_ns) {
  struct timespec ts;
  int rc;

  if (scheduled_start_ns < 0) {
    aos_set_error("aos_sleep_until: invalid time");
    return AOS_ERR_ARGS;
  }

  ts.tv_sec = (time_t)(scheduled_start_ns / 1000000000LL);
  ts.tv_nsec = (long)(scheduled_start_ns % 1000000000LL);

  do {
    rc = clock_nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME, &ts, NULL);
  } while (rc == EINTR);

  if (rc != 0) {
    snprintf(g_last_error, sizeof(g_last_error),
             "aos_sleep_until: clock_nanosleep failed: %d", rc);
    return AOS_ERR_SLEEP;
  }
  return AOS_OK;
}

int aos_finalize(const char *run_key, const char *state_dir, int world_rank) {
  char path[512];

  aos_unmap_local();

  /* Only world rank 0 unlinks the shared object. */
  if (world_rank == 0) {
    if (run_key && run_key[0]) {
      aos_build_path(path, sizeof(path), state_dir, run_key);
    } else {
      snprintf(path, sizeof(path), "%s", g_shm_path);
    }
    if (path[0]) {
      unlink(path);
    }
  }

  g_shm_path[0] = '\0';
  g_num_groups = 0;
  g_run_key[0] = '\0';
  aos_set_error("");
  return AOS_OK;
}

int aos_get_hostname(char *buf, int buflen) {
  if (!buf || buflen <= 1) {
    return -1;
  }
  if (gethostname(buf, (size_t)buflen) != 0) {
    return -1;
  }
  buf[buflen - 1] = '\0';
  return 0;
}
