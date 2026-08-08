/*=====================================================================
 * Unit tests for arrival-order reservation algorithm and shm scheduler.
 *
 * Build:
 *   gcc -O2 -Wall -Wextra -o test_arrival_order_scheduler \
 *     arrival_order_scheduler.c test_arrival_order_scheduler.c -lpthread
 *=====================================================================*/

#include "arrival_order_scheduler.h"

#include <fcntl.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/wait.h>
#include <unistd.h>

static int g_failures = 0;

static void expect_near(const char *name, double got, double want, double tol) {
  if (fabs(got - want) > tol) {
    fprintf(stderr, "FAIL %s: got=%.9f want=%.9f\n", name, got, want);
    g_failures++;
  } else {
    printf("PASS %s\n", name);
  }
}

static void expect_ll(const char *name, long long got, long long want) {
  if (got != want) {
    fprintf(stderr, "FAIL %s: got=%lld want=%lld\n", name, got, want);
    g_failures++;
  } else {
    printf("PASS %s\n", name);
  }
}

static void expect_int(const char *name, int got, int want) {
  if (got != want) {
    fprintf(stderr, "FAIL %s: got=%d want=%d\n", name, got, want);
    g_failures++;
  } else {
    printf("PASS %s\n", name);
  }
}

static void expect_ok(const char *name, int rc) {
  if (rc != AOS_OK) {
    fprintf(stderr, "FAIL %s: rc=%d err=%s\n", name, rc, aos_last_error());
    g_failures++;
  } else {
    printf("PASS %s\n", name);
  }
}

static void test_close_arrivals(void) {
  long long arrivals[] = {0, 100000000LL, 200000000LL, 300000000LL};
  long long last = 0;
  long long scheduled;
  double delay;
  int i;
  double expect_start[] = {0.0, 0.5, 1.0, 1.5};
  double expect_delay[] = {0.0, 0.4, 0.8, 1.2};
  char name[64];

  for (i = 0; i < 4; ++i) {
    aos_compute_reservation(arrivals[i], last, i == 0, 0.5, &scheduled, &delay);
    snprintf(name, sizeof(name), "t1_start%d", i);
    expect_near(name, (double)scheduled * 1e-9, expect_start[i], 1e-9);
    snprintf(name, sizeof(name), "t1_delay%d", i);
    expect_near(name, delay, expect_delay[i], 1e-9);
    last = scheduled;
  }
}

static void test_already_spaced(void) {
  long long arrivals[] = {0, 100000000LL, 2000000000LL, 3000000000LL};
  long long last = 0;
  long long scheduled;
  double delay;
  int i;
  double expect_start[] = {0.0, 0.5, 2.0, 3.0};
  double expect_delay[] = {0.0, 0.4, 0.0, 0.0};
  char name[64];

  for (i = 0; i < 4; ++i) {
    aos_compute_reservation(arrivals[i], last, i == 0, 0.5, &scheduled, &delay);
    snprintf(name, sizeof(name), "t2_start%d", i);
    expect_near(name, (double)scheduled * 1e-9, expect_start[i], 1e-9);
    snprintf(name, sizeof(name), "t2_delay%d", i);
    expect_near(name, delay, expect_delay[i], 1e-9);
    last = scheduled;
  }
}

static void test_reverse_order(void) {
  struct {
    int gid;
    long long arrival;
    double expect_start;
    double expect_delay;
  } cases[] = {
      {1, 0LL, 0.0, 0.0},
      {0, 800000000LL, 0.8, 0.0},
      {3, 900000000LL, 1.3, 0.4},
      {2, 2000000000LL, 2.0, 0.0},
  };
  long long last = 0;
  long long scheduled;
  double delay;
  int i;
  char name[64];

  for (i = 0; i < 4; ++i) {
    aos_compute_reservation(cases[i].arrival, last, i == 0, 0.5,
                            &scheduled, &delay);
    snprintf(name, sizeof(name), "t3_g%d_start", cases[i].gid);
    expect_near(name, (double)scheduled * 1e-9, cases[i].expect_start, 1e-9);
    snprintf(name, sizeof(name), "t3_g%d_delay", cases[i].gid);
    expect_near(name, delay, cases[i].expect_delay, 1e-9);
    last = scheduled;
  }
}

static void test_n_zero(void) {
  long long arrivals[] = {0, 100000000LL, 250000000LL, 400000000LL};
  long long last = 0;
  long long scheduled;
  double delay;
  int i;
  char name[64];

  for (i = 0; i < 4; ++i) {
    aos_compute_reservation(arrivals[i], last, i == 0, 0.0, &scheduled, &delay);
    snprintf(name, sizeof(name), "t4_start%d", i);
    expect_ll(name, scheduled, arrivals[i]);
    snprintf(name, sizeof(name), "t4_delay%d", i);
    expect_near(name, delay, 0.0, 1e-12);
    last = scheduled;
  }
}

static void test_shm_rollover_and_double(void) {
  const char *key = "unittest_aos_ck";
  const char *dir = "/dev/shm";
  char path[256];
  int rc;
  long long a, s;
  int seq;
  double d;

  snprintf(path, sizeof(path), "%s/specfem_aos_%s", dir, key);
  unlink(path);

  rc = aos_init(key, dir, 4, 0, 1);
  expect_ok("t5_init", rc);

  rc = aos_reserve(0, 0, 0.5, &a, &s, &seq, &d);
  expect_ok("t5_res_ckpt0_g0", rc);
  expect_int("t5_seq0", seq, 0);

  rc = aos_reserve(0, 1, 0.5, &a, &s, &seq, &d);
  expect_ok("t5_res_ckpt0_g1", rc);
  expect_int("t5_seq1", seq, 1);

  rc = aos_reserve(1, 2, 0.5, &a, &s, &seq, &d);
  expect_ok("t5_res_ckpt1_g2", rc);
  expect_int("t5_seq_reset", seq, 0);

  rc = aos_reserve(1, 2, 0.5, &a, &s, &seq, &d);
  expect_int("t6_double_reserve", rc, AOS_ERR_DOUBLE_RESERVE);

  aos_finalize(key, dir, 0);
}

static void test_stale_shm(void) {
  const char *key = "unittest_aos_stale";
  const char *dir = "/dev/shm";
  char path[256];
  int fd;
  int rc;
  char junk[64];

  snprintf(path, sizeof(path), "%s/specfem_aos_%s", dir, key);
  unlink(path);
  fd = open(path, O_RDWR | O_CREAT | O_TRUNC, 0600);
  if (fd < 0) {
    fprintf(stderr, "FAIL t7_open_stale\n");
    g_failures++;
    return;
  }
  memset(junk, 0xAB, sizeof(junk));
  if (write(fd, junk, sizeof(junk)) != (ssize_t)sizeof(junk)) {
    fprintf(stderr, "FAIL t7_write_stale\n");
    g_failures++;
    close(fd);
    return;
  }
  close(fd);

  rc = aos_init(key, dir, 4, 0, 0);
  expect_int("t7_detect_stale", rc, AOS_ERR_CORRUPT);

  rc = aos_init(key, dir, 4, 0, 1);
  expect_ok("t7_reinit", rc);
  aos_finalize(key, dir, 0);
}

typedef struct {
  unsigned int magic;
  int version;
  pthread_mutex_t mutex;
  int initialized;
  int current_checkpoint_id;
  int num_groups;
  int next_sequence;
  int seen[64];
  long long last_reserved_start_ns;
  int last_reserved_valid;
  int owner_pid;
  char run_key[128];
} aos_state_layout_t;

static void lock_and_die(const char *path) {
  int fd = open(path, O_RDWR);
  aos_state_layout_t *st;
  if (fd < 0) {
    _exit(2);
  }
  st = mmap(NULL, sizeof(*st), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
  if (st == MAP_FAILED) {
    _exit(3);
  }
  pthread_mutex_lock(&st->mutex);
  _exit(0);
}

static void test_owner_death(void) {
  const char *key = "unittest_aos_owner";
  const char *dir = "/dev/shm";
  char path[256];
  pid_t pid;
  int status;
  int rc;
  long long a, s;
  int seq;
  double d;

  snprintf(path, sizeof(path), "%s/specfem_aos_%s", dir, key);
  unlink(path);

  rc = aos_init(key, dir, 4, 0, 1);
  expect_ok("t8_init", rc);

  /* Unmap locally but keep shm file for child. */
  aos_finalize(key, dir, /*world_rank=*/1);

  rc = aos_init(key, dir, 4, 0, 0);
  expect_ok("t8_reattach", rc);

  pid = fork();
  if (pid == 0) {
    lock_and_die(path);
  }
  waitpid(pid, &status, 0);

  rc = aos_reserve(0, 0, 0.1, &a, &s, &seq, &d);
  expect_ok("t8_recover_after_owner_death", rc);

  aos_finalize(key, dir, 0);
}

static void test_dead_owner_stale_then_replace(void) {
  const char *key = "unittest_aos_dead_replace";
  const char *dir = "/dev/shm";
  char path[256];
  int rc;
  long long a, s;
  int seq;
  double d;
  pid_t pid;
  int status;

  snprintf(path, sizeof(path), "%s/specfem_aos_%s", dir, key);
  unlink(path);

  rc = aos_init(key, dir, 4, 0, 1);
  expect_ok("t9_init", rc);
  rc = aos_reserve(1, 2, 0.5, &a, &s, &seq, &d);
  expect_ok("t9_reserve_group2", rc);
  /* Keep leftover shm with seen[2]=1 but drop local mapping without unlink. */
  aos_finalize(key, dir, /*world_rank=*/1);

  /* Poison owner_pid so attachers treat the leftover as dead/stale. */
  {
    int fd = open(path, O_RDWR);
    aos_state_layout_t *st;
    if (fd < 0) {
      fprintf(stderr, "FAIL t9_open_poison\n");
      g_failures++;
      return;
    }
    st = mmap(NULL, sizeof(*st), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (st == MAP_FAILED) {
      fprintf(stderr, "FAIL t9_mmap_poison\n");
      g_failures++;
      close(fd);
      return;
    }
    st->owner_pid = 1; /* init is never a valid SPECFEM owner */
    munmap(st, sizeof(*st));
    close(fd);
  }

  pid = fork();
  if (pid == 0) {
    /* Attacher starts first and must not stick to the poisoned inode. */
    usleep(50000);
    rc = aos_init(key, dir, 4, 1, 0);
    if (rc != AOS_OK) {
      _exit(10);
    }
    rc = aos_reserve(1, 2, 0.5, &a, &s, &seq, &d);
    aos_finalize(key, dir, 1);
    _exit(rc == AOS_OK ? 0 : 11);
  }

  usleep(100000);
  rc = aos_init(key, dir, 4, 0, 1);
  expect_ok("t9_creator_replace", rc);
  waitpid(pid, &status, 0);
  if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
    fprintf(stderr, "FAIL t9_attacher_child status=%d\n", status);
    g_failures++;
  } else {
    printf("PASS t9_attacher_avoids_stale\n");
  }

  /* Fresh object: another group can still reserve. */
  rc = aos_reserve(1, 0, 0.5, &a, &s, &seq, &d);
  expect_ok("t9_creator_reserve_group0", rc);
  aos_finalize(key, dir, 0);
}

int main(void) {
  printf("=== arrival-order scheduler unit tests ===\n");
  test_close_arrivals();
  test_already_spaced();
  test_reverse_order();
  test_n_zero();
  test_shm_rollover_and_double();
  test_stale_shm();
  test_owner_death();
  test_dead_owner_stale_then_replace();

  if (g_failures == 0) {
    printf("ALL TESTS PASSED\n");
    return 0;
  }
  printf("%d FAILURES\n", g_failures);
  return 1;
}
