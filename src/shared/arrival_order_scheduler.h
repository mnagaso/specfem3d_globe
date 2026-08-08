/*=====================================================================
 *
 *  Barrier-free arrival-order I/O reservation scheduler (single-node).
 *  Shared state lives in /dev/shm via POSIX shared memory.
 *
 *=====================================================================*/

#ifndef ARRIVAL_ORDER_SCHEDULER_H
#define ARRIVAL_ORDER_SCHEDULER_H

#ifdef __cplusplus
extern "C" {
#endif

#define AOS_OK                 0
#define AOS_ERR_ARGS          -1
#define AOS_ERR_SHM           -2
#define AOS_ERR_MUTEX         -3
#define AOS_ERR_MULTIHOST     -4
#define AOS_ERR_DOUBLE_RESERVE -5
#define AOS_ERR_CORRUPT       -6
#define AOS_ERR_NOT_INIT      -7
#define AOS_ERR_SLEEP         -8
#define AOS_ERR_FALLBACK      -9

#define AOS_MAX_GROUPS        64
#define AOS_MAX_RUN_KEY       128
#define AOS_MAX_HOSTNAME      128
#define AOS_MAGIC             0xA0510DE5u
#define AOS_VERSION           1

/* Pure algorithm helper (also used by unit tests). */
void aos_compute_reservation(
    long long arrival_ns,
    long long last_reserved_start_ns,
    int is_first,
    double spacing_sec,
    long long *scheduled_start_ns,
    double *assigned_delay_sec);

int aos_init(
    const char *run_key,
    const char *state_dir,
    int num_groups,
    int world_rank,
    int create_new);

int aos_reserve(
    int checkpoint_id,
    int group_id,
    double spacing_sec,
    long long *arrival_ns,
    long long *scheduled_start_ns,
    int *arrival_sequence,
    double *assigned_delay_sec);

int aos_sleep_until(long long scheduled_start_ns);

long long aos_now_ns(void);

int aos_finalize(const char *run_key, const char *state_dir, int world_rank);

/* Expose last error message (static buffer). */
const char *aos_last_error(void);

/* Hostname helper for single-node checks. */
int aos_get_hostname(char *buf, int buflen);

#ifdef __cplusplus
}
#endif

#endif /* ARRIVAL_ORDER_SCHEDULER_H */
