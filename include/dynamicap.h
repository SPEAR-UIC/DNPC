#ifdef __cplusplus
extern "C" {
#endif

/* include guard */
#ifndef DNPC_H
#define DNPC_H

#include "PoLiMEr.h"
#include "papi.h"


/* platform-specific constants */
#ifdef BEBOP
#define TDP         120
#define MIN_CAP     50
#define MAX_CAP     120
#define MAX_FREQ    1500.0
#define MAX_IPC     1.5
#define RAPL_BUFFER 10.0
#define MAX_PAPI_EVENTS 5
#endif

#ifdef THETA
#define TDP           215.0
#define MIN_CAP       75.0
#define MAX_CAP       215.0
#define MAX_FREQ      1400.0
#define MAX_IPC       1.5
#define RAPL_BUFFER   10.0
#define MAX_PAPI_EVENTS 5
#endif

/* struct for signal handler data */
struct data_point {
  int valid; // useful for not printing potentially bad measurements
  double timestamp;
  double tss; // time since start
  double pkg_energy;
  double pkg_power;
  double dram_energy;
  double dram_power;
  double poli_frequency;
  double papi_frequency;
  double power_cap;
  double est_normal_perf;    
  double est_current_perf;   
  double est_relative_perf;  
  long long int *papi_values;
};

/* interface */
int dnpc_init(void);
int dnpc_finalize(void);

#ifdef AIX
void dnpc_init_fort(void);
void dnpc_finalize_fort(void);
#else
void dnpc_init_fort_(void);
void dnpc_finalize_fort_(void);
#endif

/* algorithm function signatures */
double algorithm_static_cap(void);
double algorithm_frequency_ipc(double, double, double, double);
double algorithm_freq_ipc_perf(double, double, double, double, double);
double algorithm_freq_uops_perf(double, double, double, double, double, double);
double algorithm_uops_freq_perf(double, double, double, double, double, double);
double algorithm_perf_ipc(double, double, double, double, double, double); // alg 4
double algorithm_perf_uops(double, double, double, double, double, double); // alg 5
double algorithm_new_ipc(double, double, double, double, double, double); // alg 6
double algorithm_new_uops(double, double, double, double, double, double); // alg 7

#endif
#ifdef __cplusplus
}
#endif
