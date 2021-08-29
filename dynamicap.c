#include <time.h>
#include <sys/time.h>
#ifndef _NOMPI
#include <mpi.h>
#endif
#ifndef _NOOMP
#include <omp.h>
#endif
#include "PoLiMEr.h"
#include "papi.h"
#include "dnpc.h"

/* helper functions */
static void custom_handler(int);
static void start_timer(void);
static void stop_timer(void);
static int initialize_papi(void);
static int initialize_polimer(void);
static int finalize_papi(void);
static int finalize_polimer(void);
static int write_output_files(void);

/* globals */
// environment variables
static char *EXP;
static char *CONF;
static char *TRIAL;
static char *PAPI_EVENTS;
static int ALGORITHM;
static int DEBUG_FLAG;
static int MAX_VALUES;
static int SAMPLE_SIZE;
static int MY_POLL_INTERVAL;     // given in ms, so multiply by 1000 to get us
static double IPC_STEP;          // how much to increase/decrease boost sensitivity
static double FREQ_STEP;         // how much to increase/decrease degradation sensitivity
static double UOPS_STEP;         // how much of a change in uops_ratio is needed
static double POWER_STEP;        // how much to increase/decrease power cap change
static double PERF_DEG_THRESHOLD;// threshold of performance degradation allowed
// global timer variables (needed by start and stop)
static struct sigaction *act;
static struct itimerval *timer;
static double timer_start_time;
// global papi variables (needed by init and finalize)
//static char *event_names[MAX_PAPI_EVENTS];
static char **event_names;
static int num_events;
static int eventset;
// handler variables (need to save state between handler calls)
static int *poll_thread;
static struct energy_reading *current_energy;
static struct data_point *dcap_data;
static double *uops_ratio_samples;
static double *power_samples;
static double *freq_samples;
static double *ipc_samples;
static double prev_uops_ratio_avg; // assume 1:1 normal instruction
static double prev_power_avg;      // assume max power
static double prev_power_cap;      // assume max power
static double prev_freq_avg;       // assume max frequency
static double prev_papi_ins;        
static double prev_papi_cyc;
static double prev_ipc_avg;        // assume no ipc
static int measurements;
static int timer_index;
static int steps_from_last_valid;
// algorithm variables (may need to save state between calls)
static double power_change;      // how much to increase/decrease power cap (start conservative)
static double ipc_give_or_take; 
static double uops_give_or_take; 
static double freq_give_or_take; // make comparisons with a give or take percentage


/* custom handler */
void
custom_handler(int sig)
{
  if (DEBUG_FLAG) {
    fprintf(stderr, "Custom handler called\n");
  }
  double timestamp;
  double tss; 
  double pkg_energy;
  double pkg_power;
  double dram_energy;
  double dram_power;
  double poli_frequency;
  double papi_frequency;
  double power_cap;
  double ipc = 0;
  double papi_cyc = 0;
  double papi_ins = 0; 
  double papi_uops = 0; 
  double uops_ratio = 0;
  double pkg_energy_diff, dram_energy_diff, time_diff;
  double est_normal_perf, est_current_perf, est_relative_perf;
  double papi_ins_diff, papi_cyc_diff, papi_uops_diff; 
  int last_valid, valid = 1; 
  int i, rv; 
  
  /* only the chosen poll thread can use handler */
  #ifndef _NOOMP
	if (*poll_thread != omp_get_thread_num()) {
		return;
	}
  #endif

  // measure and compute power data (maybe not aysnc-safe)

  if ((rv = poli_get_current_energy(current_energy)) != 0) {
    fprintf(stderr, "ERROR: POLI couldn't read current power\n");
    valid = 0;
  } 
  if ((rv = poli_get_current_frequency(&poli_frequency)) != 0) {
    fprintf(stderr, "ERROR: POLI couldn't read current frequency \n");
    valid = 0;
  }
  if ((rv = poli_get_power_cap(&power_cap)) != 0) {
    fprintf(stderr, "ERROR: POLI couldn't get current power cap\n");
    valid = 0;
  }

  // fill in PAPI values for this data point
  long long int *temp_values = calloc(num_events, sizeof(long long int));
  //if ((rv = PAPI_read(eventset, dcap_data[timer_index].papi_values)) != PAPI_OK) {
  if ((rv = PAPI_read(eventset, temp_values)) != PAPI_OK) {
    fprintf(stderr, "PAPI_ERROR: %d = %s \n", rv, PAPI_strerror(rv));
    valid = 0;
  }
  
  // assign values to be used in further calculations (async-safe)
  #ifndef _NOMPI
  timestamp   = MPI_Wtime();
  #else
  struct timeval tv;
  gettimeofday(&tv, NULL);
  timestamp = ((double) tv.tv_sec) + ((double) tv.tv_usec / 1000000);
  #endif
  tss         = timestamp - timer_start_time;
  pkg_energy  = current_energy->rapl_energy.package;
  dram_energy = current_energy->rapl_energy.dram;
  papi_cyc    = temp_values[0]; // assume first event given is cycles
  papi_ins    = temp_values[1]; // assume second event given is instructions
  papi_uops   = temp_values[3]; // assume second event given is instructions
  // calculate last valid index for next poll's calculations
  //last_valid = timer_index - steps_from_last_valid;
  last_valid = timer_index - 1; // no longer updates timer_index, so led to segfaults

  if (timer_index > 0) {
    time_diff = tss - dcap_data[last_valid].tss;
    /* error checking */
	if (DEBUG_FLAG) {
	  fprintf(stderr, "about to check papi poll \n");
	}
    for (i = 0; i < num_events; i++) {
      if (temp_values[i] <= 0) {
        valid = 0;
        if (DEBUG_FLAG) {
          fprintf(stderr, "ERROR: bad papi poll \n");
        }
      }
    }
      if (DEBUG_FLAG) {
        fprintf(stderr, "about to check if energy is negative \n");
      }
    if (pkg_energy < 0 || dram_energy < 0) {
      valid = 0;
      if (DEBUG_FLAG) {
        fprintf(stderr, "ERROR: negative energy\n");
      }
    }
	  if (DEBUG_FLAG) {
		fprintf(stderr, "about to check if time difference is too small \n");
	  }
    if (time_diff < (((double) MY_POLL_INTERVAL) / 2.0) * 0.000001) {
      valid = 0;
      if (DEBUG_FLAG) {
        fprintf(stderr, "ERROR: time_diff too small\n");
      }
    }
    if (DEBUG_FLAG) {
	#ifndef _NOOMP
      fprintf(stderr, "Poll %d, Thread %d, Time since start: %f, Valid?: %d, steps_from_last_valid: %d\n", 
		  timer_index, omp_get_thread_num(), tss, valid, steps_from_last_valid);
	#else
      fprintf(stderr, "Poll %d, Time since start: %f, Valid?: %d\n", timer_index, tss, valid);
	#endif
    }
    if (valid) { 
      if (DEBUG_FLAG) {
        fprintf(stderr, "about to copy papi values over\n");
      }
      steps_from_last_valid = 1;
	  // copy valid temporary papi values into permanent record for this poll
      memcpy(dcap_data[timer_index].papi_values, temp_values, (sizeof(long long int) * num_events));
      free(temp_values);
      if (DEBUG_FLAG) {
        fprintf(stderr, "copied papi values over\n");
      }
    } else {
      steps_from_last_valid++;
      free(temp_values);
      // store dummy results (to not corrupt memory) except for valid field
      // everything is already 0 from calloc (and 0 in the valid field means invalid)
      // by returning, we leave the data at that index with the zero'd out values
      return;
    }
    if (steps_from_last_valid > SAMPLE_SIZE) {
      // TODO: might mess up average calculations
    }
  } else {
    // on first poll, just gotta hope the data is valid :)
    steps_from_last_valid = 1;
    memcpy(dcap_data[timer_index].papi_values, temp_values, (sizeof(long long int) * num_events));
    free(temp_values);
  }

  // need to check if first poll, since first data point has nothing to compare to
  if (timer_index > 0) {
    pkg_energy_diff  = pkg_energy  - dcap_data[last_valid].pkg_energy;
    dram_energy_diff = dram_energy - dcap_data[last_valid].dram_energy;
    pkg_power        = (pkg_energy_diff  / time_diff);
    dram_power       = (dram_energy_diff / time_diff);
    papi_cyc_diff    = papi_cyc    - dcap_data[last_valid].papi_values[0];
    papi_ins_diff    = papi_ins    - dcap_data[last_valid].papi_values[1];
    papi_uops_diff   = papi_uops   - dcap_data[last_valid].papi_values[3];
    papi_frequency   = (double) ((papi_cyc_diff / time_diff) / 1000000);
    ipc              = (double) (papi_ins_diff  / papi_cyc_diff);
    uops_ratio       = (double) (papi_uops_diff / papi_ins_diff);
    // if frequencies are the same, then perf. deg.
    // otherwise, papi_freq showing something from no cap profile
    if (poli_frequency - papi_frequency < (0.05 * MAX_FREQ)) { 
      est_normal_perf    = dcap_data[last_valid].est_normal_perf + MAX_FREQ;
    } else { 
      est_normal_perf    = dcap_data[last_valid].est_normal_perf + papi_frequency;
    }
    est_current_perf     = dcap_data[last_valid].est_current_perf + papi_frequency; 
  } else {
    uops_ratio     = (double) (papi_uops / papi_ins);
    ipc            = (double) (papi_ins / papi_cyc);
    papi_frequency = (double) ((papi_cyc / tss) / 1000000);
    pkg_power      = (double) (pkg_energy / tss); // just an assumption
    dram_power     = (double) (dram_energy / tss); 
    est_normal_perf    = papi_frequency;
    est_current_perf   = papi_frequency; 
  }
  est_relative_perf = 1 + ((est_normal_perf - est_current_perf) / est_normal_perf);


  // assign more values (async-safe)
  dcap_data[timer_index].valid       = valid;
  dcap_data[timer_index].timestamp   = timestamp;
  dcap_data[timer_index].tss         = tss;
  dcap_data[timer_index].pkg_energy  = pkg_energy;
  dcap_data[timer_index].pkg_power   = pkg_power;
  dcap_data[timer_index].dram_energy = dram_energy;
  dcap_data[timer_index].dram_power  = dram_power;
  dcap_data[timer_index].power_cap   = power_cap;
  dcap_data[timer_index].poli_frequency = poli_frequency;
  dcap_data[timer_index].papi_frequency = papi_frequency;

  // assign performance values (async-safe)
  dcap_data[timer_index].est_normal_perf   = est_normal_perf; 
  dcap_data[timer_index].est_current_perf  = est_current_perf; 
  dcap_data[timer_index].est_relative_perf = est_relative_perf;

  if (DEBUG_FLAG) {
    fprintf(stderr, "assigned values to data struct \n");
  }

  /* Store basic metric calculations for each algorithm*/
  uops_ratio_samples[measurements] = uops_ratio;
  power_samples[measurements] = dcap_data[timer_index].pkg_power;
  freq_samples[measurements]  = poli_frequency; // use poli_frequency for state machine
  ipc_samples[measurements]   = ipc;
  measurements += 1;

  if (DEBUG_FLAG) {
	fprintf(stderr, "assigned measurement values into sample \n");
  }

  // once we reach x samples, call algorithm for cap
  if (measurements == SAMPLE_SIZE) {
    if (DEBUG_FLAG) {
      fprintf(stderr, "collected %d measurements, calling algorithm\n", measurements);
    }

    double uops_ratio_avg;
    double power_avg;
    double freq_avg;
    double ipc_avg;   

    // calculate averages for different metrics
    uops_ratio_avg = 0;
    power_avg = 0;
    freq_avg  = 0;
    ipc_avg   = 0;
    for (i = 0; i < SAMPLE_SIZE; i++) {
      uops_ratio_avg = uops_ratio_avg + uops_ratio_samples[i];
      power_avg = power_avg + power_samples[i];
      freq_avg  = freq_avg  + freq_samples[i];
      ipc_avg   = ipc_avg   + ipc_samples[i];
    }
    uops_ratio_avg = uops_ratio_avg / SAMPLE_SIZE;
    power_avg = power_avg / SAMPLE_SIZE;
    freq_avg  = freq_avg  / SAMPLE_SIZE;
    ipc_avg   = ipc_avg   / SAMPLE_SIZE;

    // call algorithm for power cap
    // only after first period to get rid of startup errors (initial power assumptions)
    if (timer_index > SAMPLE_SIZE) { 
      switch(ALGORITHM) {
        case 100: 
          algorithm_static_cap();
          break;
        case 0:
          algorithm_frequency_ipc(timestamp, power_avg, freq_avg, ipc_avg);
          break;
        case 1:
          algorithm_freq_ipc_perf(timestamp, power_avg, freq_avg, ipc_avg, dcap_data[timer_index].est_relative_perf);
          break;
        case 2:
          algorithm_freq_uops_perf(timestamp, power_avg, freq_avg, ipc_avg, uops_ratio_avg, dcap_data[timer_index].est_relative_perf);
          break;
        case 3:
          algorithm_uops_freq_perf(timestamp, power_avg, freq_avg, ipc_avg, uops_ratio_avg, dcap_data[timer_index].est_relative_perf);
          break;
        case 4:
          algorithm_perf_ipc(timestamp, power_avg, freq_avg, ipc_avg, uops_ratio_avg, dcap_data[timer_index].est_relative_perf);
          break;
        case 5:
          algorithm_perf_uops(timestamp, power_avg, freq_avg, ipc_avg, uops_ratio_avg, dcap_data[timer_index].est_relative_perf);
          break;
        case 6:
          algorithm_new_ipc(timestamp, power_avg, freq_avg, ipc_avg, uops_ratio_avg, dcap_data[timer_index].est_relative_perf);
          break;
        case 7:
          algorithm_new_uops(timestamp, power_avg, freq_avg, ipc_avg, uops_ratio_avg, dcap_data[timer_index].est_relative_perf);
          break;
        default:
          if (DEBUG_FLAG) {
            fprintf(stderr, "No power capping algorithm selected\n");
          }
      }
    }
    // reset measurement counter
    measurements = 0;
  } else {
	  if (DEBUG_FLAG) {
		fprintf(stderr, "measurement %d taken\n", measurements);
	  }
  }
  // update index for next poll
  timer_index += 1;

}

/* timer functions */
void
start_timer()
{
  if (DEBUG_FLAG) {
    fprintf(stderr, "start_timer called\n");
  }

  timer_index = 0;
#ifndef _NOMPI
  timer_start_time = MPI_Wtime();
#else
  struct timeval tv;
  gettimeofday(&tv, NULL);
  timer_start_time = ((double) tv.tv_sec) + ((double) tv.tv_usec / 1000000);
#endif

  //act    = (struct sigaction *) malloc(sizeof(struct sigaction));
  //timer  = (struct itimerval *) malloc(sizeof(struct itimerval));
  act    = malloc(sizeof(struct sigaction));
  timer  = malloc(sizeof(struct itimerval));

  act->sa_handler = &custom_handler;
  sigaction(SIGALRM, act, NULL);
  
  // using microseconds instead of seconds to get <1 second polling
  timer->it_value.tv_sec     = 0;
  timer->it_value.tv_usec    = MY_POLL_INTERVAL;
  timer->it_interval.tv_sec  = 0;
  timer->it_interval.tv_usec = MY_POLL_INTERVAL;

  setitimer(ITIMER_REAL, timer, NULL);
  
  if (DEBUG_FLAG) {
    fprintf(stderr, "start_timer finished \n");
  }

  return;
}

void
stop_timer()
{
  if (DEBUG_FLAG) {
    fprintf(stderr, "stop_timer called \n");
  }
  timer->it_value.tv_sec  = 0;
  timer->it_value.tv_usec = 0;
  setitimer(ITIMER_REAL, timer, NULL);

  free(act);
  free(timer);

  if (DEBUG_FLAG) {
    fprintf(stderr, "stop_timer finished \n");
  }

  return;
}

int
initialize_polimer()
{
  if (poli_init() != 0) {
    fprintf(stderr, "PoliMEr_ERROR: could not initialize\n");
    return 1;
  }
  return 0;
}

int
initialize_papi()
{
  char *temp;
  int eventcode; 
  int rv, i;

  // start PAPI
  num_events = 0;
  eventset  = PAPI_NULL;
  event_names = calloc(MAX_PAPI_EVENTS, sizeof(char*));
  for (i = 0; i < MAX_PAPI_EVENTS; i++) {
    event_names[i] = calloc(50, sizeof(char));
  }

  if ((rv = PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT) {
    fprintf(stderr,"PAPI_ERROR: couldnt initialize\n");
    fprintf(stderr,"PAPI_ERROR: %d = %s \n", rv, PAPI_strerror(rv));
    return 1;
  }
  if ((rv = PAPI_create_eventset(&eventset)) != PAPI_OK) {
    fprintf(stderr,"PAPI_ERROR: couldnt create eventset\n");
    fprintf(stderr,"PAPI_ERROR: %d = %s \n", rv, PAPI_strerror(rv));
    return 1;
  }

  // check and parse events from environment variable
  if (PAPI_EVENTS == NULL) {
    printf("PAPI_ERROR: could not get any papi events from environment variable\n");
    return 1;
  }

  // parse events and store in allocated string array
  temp = strtok(PAPI_EVENTS, ",");
  do {
    strcpy(event_names[num_events], temp);
    num_events++;
  } while ((temp = strtok(NULL, ",")) != NULL);
  

  // now we know how many events passed in, malloc value array for each data point
  if (DEBUG_FLAG) {
    fprintf(stderr, "about to malloc papi value array in each data point in data array\n");
  }
  fprintf(stderr,"PAPI: %d events found in environment variable \n", num_events);
  for (i = 0; i < MAX_VALUES; i++) {
    //dcap_data[i].papi_values = (long long int *) calloc(num_events, sizeof(long long int));
    dcap_data[i].papi_values = calloc(num_events, sizeof(long long int));
  }

  // add events to event set
  for (i = 0; i < num_events; i++) {
    if ((rv = PAPI_event_name_to_code(event_names[i], &eventcode)) != PAPI_OK) {
      fprintf(stderr,"PAPI_ERROR: couldn't find eventcode for event %s \n", event_names[i]);
      fprintf(stderr,"PAPI_ERROR: %d = %s \n", rv, PAPI_strerror(rv));
    } else if((rv = PAPI_add_event(eventset, eventcode)) != PAPI_OK) {
      fprintf(stderr,"PAPI_ERROR: couldn't add eventcode %d for event %s \n", eventcode, event_names[i]);
      fprintf(stderr,"PAPI_ERROR: %d = %s \n", rv, PAPI_strerror(rv));
    } else {
      fprintf(stderr,"PAPI: added event %d = %s successfully to event set \n", eventcode, event_names[i]);
    }
  }

  // start counting events here
  if (DEBUG_FLAG) {
    fprintf(stderr, "about to start counting from PAPI event set\n");
  }
  if ((rv = PAPI_start(eventset)) != PAPI_OK) {
    fprintf(stderr,"PAPI_ERROR: couldn't start reading from eventset\n");
    fprintf(stderr,"PAPI_ERROR: %d = %s \n", rv, PAPI_strerror(rv));
    return 1;
  }

  return 0;
}

int 
dnpc_init() 
{
  /* initialize globals from environemt */
  EXP                = getenv("EXPERIMENT");
  CONF               = getenv("CONFIGURATION");
  TRIAL              = getenv("TRIAL");
  IPC_STEP           = (atof(getenv("IPC_GIVE_OR_TAKE"))  / 100.0) * MAX_IPC;  
  FREQ_STEP          = (atof(getenv("FREQ_GIVE_OR_TAKE")) / 100.0) * MAX_FREQ; 
  UOPS_STEP          = (atof(getenv("UOPS_GIVE_OR_TAKE")) / 100.0);                    
  ALGORITHM          = atoi(getenv("ALGORITHM")); 
  POWER_STEP         = (atof(getenv("POWER_STEP")) / 100.0)        * TDP;      
  MAX_VALUES         = atoi(getenv("MAX_VALUES"));
  DEBUG_FLAG         = atoi(getenv("DEBUG_FLAG"));
  SAMPLE_SIZE        = atoi(getenv("SAMPLE_SIZE"));
  PAPI_EVENTS        = getenv("PAPI_EVENTS");
  MY_POLL_INTERVAL   = 1000 * atoi(getenv("MY_POLL_INTERVAL")); 
  PERF_DEG_THRESHOLD = (atof(getenv("PERF_DEG_THRESHOLD")) / 100.0);   
  // algorithm variables
  power_change       = RAPL_BUFFER;      
  ipc_give_or_take   = IPC_STEP;  
  freq_give_or_take  = FREQ_STEP; 
  uops_give_or_take  = UOPS_STEP; 
  // handler variables
  prev_uops_ratio_avg = 1;   // assume 1:1 normal instruction
  prev_power_avg = TDP;      // assume max power
  prev_power_cap = TDP;      // assume max power
  prev_freq_avg  = MAX_FREQ; // assume max frequency
  prev_papi_ins  = 0;        
  prev_papi_cyc  = 0;
  prev_ipc_avg   = 0;        // assume no ipc
  measurements   = 0;
  timer_index    = 0;
  steps_from_last_valid = 1;


  /* allocate data structure memory */
  if (DEBUG_FLAG) {
    fprintf(stderr, "about to malloc energy, power, and data array\n");
  }
  // note: need to cast mallocs in C++ but not in C
  dcap_data          = calloc(MAX_VALUES,  sizeof(struct data_point));
  ipc_samples        = calloc(SAMPLE_SIZE, sizeof(double));
  freq_samples       = calloc(SAMPLE_SIZE, sizeof(double));
  power_samples      = calloc(SAMPLE_SIZE, sizeof(double));
  current_energy     = calloc(1,           sizeof(struct energy_reading));
  uops_ratio_samples = calloc(SAMPLE_SIZE, sizeof(double));

  /* assign only one thread to use handler */
  #ifndef _NOOMP
	poll_thread = calloc(1, sizeof(int));
	*poll_thread = omp_get_thread_num();
  #endif

  if (DEBUG_FLAG) {
    fprintf(stderr, "malloc'd energy, power, and data array successfully\n");
  }

  /* start PAPI */
  if (DEBUG_FLAG) {
    fprintf(stderr, "about to initialize PAPI\n");
  }
  initialize_papi();

  /* start PoLiMEr */
  if (DEBUG_FLAG) {
    fprintf(stderr, "about to initialize PoLiMEr\n");
  }
  initialize_polimer();

  /* start periodic data collection */
  start_timer();

  return 0;
}

int
finalize_polimer()
{
  poli_finalize();
  return 0;
}

int
finalize_papi()
{
  int rv;
  long long int *temp_values;

  temp_values = calloc(num_events, sizeof(long long int));

  if ((rv = PAPI_stop(eventset, temp_values)) != PAPI_OK) {
    fprintf(stderr, "PAPI_ERROR: PAPI couldn't stop counting the event set and return the counted values\n");
    fprintf(stderr, "PAPI_ERROR: %d = %s \n", rv, PAPI_strerror(rv));
    return 1;
  }

  // clean up
  free(temp_values);
  return 0;
}

int
write_output_files()
{
  char *poli_filename;
  char *papi_filename;
  char temp[50];
  FILE *poli_fp;
  FILE *papi_fp;
  int i, j;

  // create output files
  poli_filename = calloc(50, sizeof(char));
  papi_filename = calloc(50, sizeof(char));
  sprintf(temp, "data/exp%s_conf%s_trial%s_POLI.txt", EXP, CONF, TRIAL);
  strcpy(poli_filename, temp);
  sprintf(temp, "data/exp%s_conf%s_trial%s_PAPI.txt", EXP, CONF, TRIAL);
  strcpy(papi_filename, temp);
  if (DEBUG_FLAG) {
    fprintf(stderr, "about to write data struct to files\n");
    fprintf(stderr, "printing to %s\n", poli_filename); 
    fprintf(stderr, "printing to %s\n", papi_filename); 
  }
  
  // open PAPI data file and POLI data file
  poli_fp = fopen(poli_filename, "w");
  if (poli_fp == NULL) {
    fprintf(stderr, "PAPI_ERROR: couldn't open poli counter output file for writing \n");
  }
  papi_fp = fopen(papi_filename, "w");
  if (papi_fp == NULL) {
    fprintf(stderr, "PAPI_ERROR: couldn't open papi counter output file for writing \n");
  }

  // write column header to PAPI data file (dynamic event names and number of events)
  if (DEBUG_FLAG) {
    fprintf(stderr, "about to write column headers to data files\n");
  }
  fprintf(papi_fp, "%s, %s,", "Valid", "Time Since Start (s)");
  for (i = 0; i < num_events; i++) {
    fprintf(papi_fp, "%s", event_names[i]);
    if (i < (num_events - 1)) {
      fprintf(papi_fp, ", ");
    } else {
      fprintf(papi_fp, "\n");
    }
  }
  // write column header to POLI data file
  fprintf(poli_fp, "%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n", 
    "Valid",
    "Time Since Start (s)", 
    "RAPL PKG Energy (J)", 
    "RAPL PKG Power (W)", 
    "RAPL DRAM Energy (J)", 
    "RAPL DRAM Power (W)", 
    "POLI CPU Frequency (MHz)",
    "PAPI CPU Frequency (MHz)",
    "PKG Power Cap (W)",
    "Normal Performance",
    "Current Performance",
    "Relative Performance"
    );

  if (DEBUG_FLAG) {
    fprintf(stderr, "about to write data to data files\n");
  }

  // write each data values to both files
  for (i = 0; i < timer_index; i++) {
    // write valid datapoint bit for both files
    fprintf(papi_fp, "%d, ", dcap_data[i].valid);
    fprintf(poli_fp, "%d, ", dcap_data[i].valid);
    // write time since start for both files
    fprintf(papi_fp, "%f, ", dcap_data[i].tss);
    fprintf(poli_fp, "%f, ", dcap_data[i].tss);
    // write papi values to papi file
    for (j = 0; j < num_events; j++) {
      fprintf(papi_fp, "%lld", dcap_data[i].papi_values[j]);
      if (j < (num_events - 1)) {
        fprintf(papi_fp, ", ");
      } else {
        fprintf(papi_fp, "\n");
      }
    }
    // write poli values to poli file
    fprintf(poli_fp, "%f, %f, %f, %f, %f, %f, %f, ", 
      dcap_data[i].pkg_energy, 
      dcap_data[i].pkg_power,
      dcap_data[i].dram_energy, 
      dcap_data[i].dram_power,
      dcap_data[i].poli_frequency,
      dcap_data[i].papi_frequency,
      dcap_data[i].power_cap
      );
    // write performance values to same line in poli file
    fprintf(poli_fp, "%f, %f, %f\n", 
      dcap_data[i].est_normal_perf, 
      dcap_data[i].est_current_perf,
      dcap_data[i].est_relative_perf
      );
  }

  // close out files to disk
  fclose(papi_fp);
  fclose(poli_fp);

  // clean up
  //free(papi_filename);
  //free(poli_filename);
  for (i = 0; i < MAX_PAPI_EVENTS; i++) {
    free(event_names[i]);
  }
  free(event_names);
  return 0;
}

int 
dnpc_finalize() 
{
  int i;

  /* stop data collection */
  if (DEBUG_FLAG) {
    fprintf(stderr, "about to stop timer \n");
  }
  stop_timer(); 
  if (DEBUG_FLAG) {
    fprintf(stderr, "about to stop PoLiMEr\n");
  }
  finalize_polimer();
  if (DEBUG_FLAG) {
    fprintf(stderr, "about to stop counting from PAPI event set \n");
  }
  finalize_papi();

  /* write to files */
  if (DEBUG_FLAG) {
    fprintf(stderr, "about to close data files\n");
  }
  write_output_files();

  /* clean up */
  #ifndef _NOOMP
	free(poll_thread);
  #endif
  if (DEBUG_FLAG) {
    fprintf(stderr, "about to free papi values for each data point\n");
  }
  for (i = 0; i < MAX_VALUES; i++) {
    free(dcap_data[i].papi_values);
  }
  if (DEBUG_FLAG) {
    fprintf(stderr, "about to free uops_ratio_samples\n");
  }
  free(uops_ratio_samples);
  if (DEBUG_FLAG) {
    fprintf(stderr, "about to free current_energy\n");
  }
  free(current_energy);
  if (DEBUG_FLAG) {
    fprintf(stderr, "about to free power_samples\n");
  }
  free(power_samples);
  if (DEBUG_FLAG) {
    fprintf(stderr, "about to free freq_samples\n");
  }
  free(freq_samples);
  if (DEBUG_FLAG) {
    fprintf(stderr, "about to free ipc_samples\n");
  }
  free(ipc_samples);
  if (DEBUG_FLAG) {
    fprintf(stderr, "about to free dcap_data\n");
  }
  free(dcap_data);
  return 0;
}


/* Algorithms Definitions */
double 
algorithm_frequency_ipc(double timestamp, double power_avg, double freq_avg, double ipc_avg)
{
  if (DEBUG_FLAG) {
    fprintf(stderr, "called frequency ipc algorithm \n");
  }

  /* variables */
  double delta_ipc_avg  = ipc_avg  - prev_ipc_avg; 
  double power_cap      = prev_power_cap;
  
  /* actual state machine */
  if (freq_avg >= MAX_FREQ - freq_give_or_take) { // State 1: candidate for lower power cap
    power_cap = power_avg; 
  } else {
    if (delta_ipc_avg < 0 - ipc_give_or_take) {   // State 2: phase change from lower power to higher power
      power_cap = TDP;
    } else {                                      // State 3: no phase change; we just capped too low
      power_cap = power_cap + POWER_STEP;
    }
  }

  /* check power cap bounds then set cap */
  if (power_cap < MIN_CAP) {
    power_cap = MIN_CAP;
  } else if (power_cap > MAX_CAP) {
    power_cap = MAX_CAP;
  }
  poli_set_power_cap(power_cap);
  
  /* store these stats as previous */
  prev_freq_avg   = freq_avg;
  prev_ipc_avg    = ipc_avg;
  prev_power_avg  = power_avg;
  prev_power_cap  = power_cap;

  return power_cap;
}

double 
algorithm_freq_ipc_perf(double timestamp, double power_avg, double freq_avg, double ipc_avg, double est_relative_perf)
{
  if (DEBUG_FLAG) {
    fprintf(stderr, "called frequency ipc performance algorithm \n");
  }

  /* statistics */
  double delta_ipc_avg  = ipc_avg  - prev_ipc_avg; 
  double power_cap      = prev_power_cap;
  double perf_power_cap;

  /* online performance */
  if (est_relative_perf > 1 + PERF_DEG_THRESHOLD) { 
    power_change = RAPL_BUFFER + POWER_STEP; // 15 W above measured should not affect measured power
    perf_power_cap = power_avg + power_change;
  } else if (est_relative_perf > 1 + (PERF_DEG_THRESHOLD * 0.50)) { 
    /* be more conservative in capping */
    power_change  += POWER_STEP; 
    perf_power_cap = power_avg + power_change;
    if (power_change > RAPL_BUFFER) { // don't want to increase above rapl power cap buffer 
      power_change = RAPL_BUFFER; 
      perf_power_cap = power_avg + power_change;
    }
  } else { 
    /* be more aggressive in capping */
    power_change  -= POWER_STEP; 
    perf_power_cap = power_avg + power_change;
    if (perf_power_cap <= MIN_CAP) { // don't want to decrease below minimum cap on consecutive min cap settings
      power_change += POWER_STEP; 
      perf_power_cap = power_avg + power_change;
    }
  }
  /* actual state machine */
  if (freq_avg >= MAX_FREQ - freq_give_or_take) { // state 1: candidate for lower power cap
    power_cap = perf_power_cap;
  } else {
    if (delta_ipc_avg < 0 - ipc_give_or_take) {   // state 2: phase change from lower power to higher power
      power_cap = TDP;
    } else {                                      // state 3: no phase change; we just capped too low
      power_cap = power_cap + POWER_STEP;
    }
  }

  /* check power cap bounds then set cap */
  if (power_cap < MIN_CAP) {
    power_cap = MIN_CAP;
  } else if (power_cap > MAX_CAP) {
    power_cap = MAX_CAP;
  }
  poli_set_power_cap(power_cap);
  
  
  /* store current stats as previous */
  prev_freq_avg   = freq_avg;
  prev_ipc_avg    = ipc_avg;
  prev_power_avg  = power_avg;
  prev_power_cap  = power_cap;

  return power_cap;
}

double 
algorithm_freq_uops_perf(double timestamp, double power_avg, double freq_avg, double ipc_avg, double uops_ratio_avg, double est_relative_perf)
{
  if (DEBUG_FLAG) {
    fprintf(stderr, "called frequency uops performance algorithm \n");
  }

  /* statistics */
  double delta_uops_ratio_avg = uops_ratio_avg - prev_uops_ratio_avg; 
  double power_cap      = prev_power_cap;

  /* online performance */
  if (est_relative_perf >= 1 + PERF_DEG_THRESHOLD) { 
    // be very conservative in capping
    power_change   = RAPL_BUFFER + POWER_STEP; 
  } else if (est_relative_perf > 1 + (PERF_DEG_THRESHOLD * 0.50)) { 
    // be more conservative in capping
    power_change  += POWER_STEP; 
    if (power_change > RAPL_BUFFER) { 
      power_change = RAPL_BUFFER; 
    }
  } else { 
    // be more aggressive in capping 
    power_change  -= POWER_STEP; 
  }
  
  /* actual state machine */
  if (freq_avg >= MAX_FREQ - freq_give_or_take) { // state 1: candidate for lower power cap
    power_cap = power_avg + power_change;
  } else {
    if (delta_uops_ratio_avg < 0 - uops_give_or_take) {   // state 2: phase change from lower power to higher power
      power_cap = TDP;
    } else {                                      // state 3: no phase change; we just capped too low
      power_cap = prev_power_cap + POWER_STEP;
    }
  }

  /* check power cap bounds then set cap */
  if (power_cap < MIN_CAP) {
    power_cap = MIN_CAP;
  } else if (power_cap > MAX_CAP) {
    power_cap = MAX_CAP;
  }
  poli_set_power_cap(power_cap);
  
  /* store current stats as previous */
  prev_uops_ratio_avg = uops_ratio_avg;
  prev_freq_avg   = freq_avg;
  prev_ipc_avg    = ipc_avg;
  prev_power_avg  = power_avg;
  prev_power_cap  = power_cap;

  return power_cap;
}



double 
algorithm_uops_freq_perf(double timestamp, double power_avg, double freq_avg, double ipc_avg, double uops_ratio_avg, double est_relative_perf)
{
  if (DEBUG_FLAG) {
    fprintf(stderr, "called uops frequency performance algorithm \n");
  }

  /* statistics */
  double delta_uops_ratio_avg = uops_ratio_avg - prev_uops_ratio_avg; 
  double power_cap      = prev_power_cap;

  /* online performance */
  if (est_relative_perf >= 1 + PERF_DEG_THRESHOLD) { 
    // be very conservative in capping
    power_change   = RAPL_BUFFER + POWER_STEP; 
  } else if (est_relative_perf > 1 + (PERF_DEG_THRESHOLD * 0.50)) { 
    // be more conservative in capping
    power_change  += POWER_STEP; 
    if (power_change > RAPL_BUFFER) { 
      power_change = RAPL_BUFFER; 
    }
  } else { 
    // be more aggressive in capping 
    power_change  -= POWER_STEP; 
  }
  
  /* actual state machine */
  if (delta_uops_ratio_avg < 0 - uops_give_or_take) { // state 2: phase change from lower power to higher power
      power_cap = TDP;
      power_change = POWER_STEP; 
  } else {
    if (freq_avg < MAX_FREQ - freq_give_or_take) {    // state 3: no phase change; we just capped too low
      power_cap = prev_power_cap + POWER_STEP;
    } else {                                          // state 1: candidate for lower power cap 
      power_cap = power_avg + power_change;
    }
  }

  /* check power cap bounds then set cap */
  if (power_cap < MIN_CAP) {
    power_cap = MIN_CAP;
  } else if (power_cap > MAX_CAP) {
    power_cap = MAX_CAP;
  }
  poli_set_power_cap(power_cap);
  
  /* store current stats as previous */
  prev_uops_ratio_avg = uops_ratio_avg;
  prev_freq_avg   = freq_avg;
  prev_ipc_avg    = ipc_avg;
  prev_power_avg  = power_avg;
  prev_power_cap  = power_cap;

  return power_cap;
}

double
algorithm_static_cap() 
{
  double power_cap;

  power_cap = atof(getenv("STATIC_POWER_CAP"));
  if (!power_cap) {
    power_cap = TDP;
  }

  /* check power cap bounds then set cap */
  if (power_cap < MIN_CAP) {
    power_cap = MIN_CAP;
  } else if (power_cap > MAX_CAP) {
    power_cap = MAX_CAP;
  }
  poli_set_power_cap(power_cap);
  
  return power_cap;
}

double 
algorithm_perf_uops(double timestamp, double power_avg, double freq_avg, double ipc_avg, double uops_ratio_avg, double est_relative_perf)
{
  if (DEBUG_FLAG) {
    fprintf(stderr, "called performance uops algorithm \n");
  }

  /* statistics */
  double delta_uops_ratio_avg = uops_ratio_avg - prev_uops_ratio_avg; 
  double power_cap            = prev_power_cap;

  /* online performance */
  if (est_relative_perf >= 1 + PERF_DEG_THRESHOLD) { 
    // be very conservative in capping
    power_change   = RAPL_BUFFER + POWER_STEP; 
  } else if (est_relative_perf > 1 + (PERF_DEG_THRESHOLD * 0.50)) { 
    // be more conservative in capping
    power_change  += POWER_STEP; 
    if (power_change > RAPL_BUFFER) { 
      power_change = RAPL_BUFFER; 
    }
  } else { 
    // be more aggressive in capping 
    power_change  -= POWER_STEP; 
  }
  
  /* actual state machine */
  if (est_relative_perf < 1 + PERF_DEG_THRESHOLD) { // state 1: candidate for lower power cap
    power_cap = power_avg + power_change;
  } else {
    if (delta_uops_ratio_avg < 0 - uops_give_or_take) { // state 2: phase change from lower to higher power
      power_cap = TDP;
    } else { // state 3: no phase change; we just capped too low
      power_cap = prev_power_cap + power_change;
    }
  }

  /* check power cap bounds then set cap */
  if (power_cap < MIN_CAP) {
    power_cap = MIN_CAP;
  } else if (power_cap > MAX_CAP) {
    power_cap = MAX_CAP;
  }
  poli_set_power_cap(power_cap);
  
  /* store current stats as previous */
  prev_uops_ratio_avg = uops_ratio_avg;
  prev_freq_avg   = freq_avg;
  prev_ipc_avg    = ipc_avg;
  prev_power_avg  = power_avg;
  prev_power_cap  = power_cap;

  return power_cap;
}


double 
algorithm_perf_ipc(double timestamp, double power_avg, double freq_avg, double ipc_avg, double uops_ratio_avg, double est_relative_perf)
{
  if (DEBUG_FLAG) {
    fprintf(stderr, "called performance ipc algorithm \n");
  }

  /* statistics */
  double delta_uops_ratio_avg = uops_ratio_avg - prev_uops_ratio_avg; 
  double delta_ipc_avg        = ipc_avg  - prev_ipc_avg; 
  double power_cap            = prev_power_cap;

  /* online performance */
  if (est_relative_perf >= 1 + PERF_DEG_THRESHOLD) { 
    // be very conservative in capping
    power_change   = RAPL_BUFFER + POWER_STEP; 
  } else if (est_relative_perf > 1 + (PERF_DEG_THRESHOLD * 0.50)) { 
    // be more conservative in capping
    power_change  += POWER_STEP; 
    if (power_change > RAPL_BUFFER) { 
      power_change = RAPL_BUFFER; 
    }
  } else { 
    // be more aggressive in capping 
    power_change  -= POWER_STEP; 
  }
  
  /* actual state machine */
  if (est_relative_perf < 1 + PERF_DEG_THRESHOLD) { 
	// state 1: candidate for lower power cap
    power_cap = power_avg + power_change;
  } else {
    if (delta_ipc_avg < 0 - ipc_give_or_take) { 
	  // state 2: phase change from lower to higher power
      power_cap = TDP;
    } else { 
	  // state 3: no phase change; we just capped too low
      power_cap = prev_power_cap + power_change;
    }
  }

  /* check power cap bounds then set cap */
  if (power_cap < MIN_CAP) {
    power_cap = MIN_CAP;
  } else if (power_cap > MAX_CAP) {
    power_cap = MAX_CAP;
  }
  poli_set_power_cap(power_cap);
  
  /* store current stats as previous */
  prev_uops_ratio_avg = uops_ratio_avg;
  prev_freq_avg   = freq_avg;
  prev_ipc_avg    = ipc_avg;
  prev_power_avg  = power_avg;
  prev_power_cap  = power_cap;

  return power_cap;
}

double 
algorithm_new_ipc(double timestamp, double power_avg, double freq_avg, double ipc_avg, double uops_ratio_avg, double est_relative_perf)
{
  if (DEBUG_FLAG) {
    fprintf(stderr, "called new ipc algorithm \n");
  }

  /* statistics */
  double delta_uops_ratio_avg = uops_ratio_avg - prev_uops_ratio_avg; 
  double delta_ipc_avg        = ipc_avg  - prev_ipc_avg; 
  double power_cap            = prev_power_cap;

  if (est_relative_perf >= (1 + PERF_DEG_THRESHOLD)) {
	// state 1: Cap that doesn't go below power
	power_cap = TDP;
  } else {
    if (delta_ipc_avg < 0 - ipc_give_or_take) {
	  // state 2: phase change from lower to higher power
      power_cap = TDP;
    } else { 
	  if (est_relative_perf >= (1 + (PERF_DEG_THRESHOLD * 0.5))) {
		// state 3: moderate degradation; cautious cap 
		power_cap = power_avg + POWER_STEP;
	  } else {
	    // state 4: no degradation; cap lower 
		power_cap = power_avg - POWER_STEP;
	  }
    }
  }

  /* check power cap bounds then set cap */
  if (power_cap < MIN_CAP) {
    power_cap = MIN_CAP;
  } else if (power_cap > MAX_CAP) {
    power_cap = MAX_CAP;
  }
  poli_set_power_cap(power_cap);
  
  /* store current stats as previous */
  prev_uops_ratio_avg = uops_ratio_avg;
  prev_freq_avg   = freq_avg;
  prev_ipc_avg    = ipc_avg;
  prev_power_avg  = power_avg;
  prev_power_cap  = power_cap;

  return power_cap;
}

double 
algorithm_new_uops(double timestamp, double power_avg, double freq_avg, double ipc_avg, double uops_ratio_avg, double est_relative_perf)
{
  if (DEBUG_FLAG) {
    fprintf(stderr, "called new uops algorithm \n");
  }

  /* statistics */
  double delta_uops_ratio_avg = uops_ratio_avg - prev_uops_ratio_avg; 
  double delta_ipc_avg        = ipc_avg  - prev_ipc_avg; 
  double power_cap            = prev_power_cap;

  /* actual state machine */
  if (est_relative_perf >= (1 + PERF_DEG_THRESHOLD)) {
	// state 1: Cap that doesn't go below power
	power_cap = TDP;
  } else {
    if (delta_uops_ratio_avg < 0 - uops_give_or_take) {
	  // state 2: phase change from lower to higher power
      power_cap = TDP;
    } else { 
	  if (est_relative_perf >= (1 + (PERF_DEG_THRESHOLD * 0.5))) {
		// state 3: moderate degradation; cautious cap 
		power_cap = power_avg + POWER_STEP;
	  } else {
	    // state 4: no degradation; cap lower 
		power_cap = power_avg - POWER_STEP;
	  }
    }
  }

  /* check power cap bounds then set cap */
  if (power_cap < MIN_CAP) {
    power_cap = MIN_CAP;
  } else if (power_cap > MAX_CAP) {
    power_cap = MAX_CAP;
  }
  poli_set_power_cap(power_cap);
  
  /* store current stats as previous */
  prev_uops_ratio_avg = uops_ratio_avg;
  prev_freq_avg   = freq_avg;
  prev_ipc_avg    = ipc_avg;
  prev_power_avg  = power_avg;
  prev_power_cap  = power_cap;

  return power_cap;
}


/* Fortran Wrapper Functions */
void
#ifdef AIX
dnpc_init_fort
#else
dnpc_init_fort_
#endif
() {  
	dnpc_init();
}

void
#ifdef AIX
dnpc_finalize_fort
#else
dnpc_finalize_fort_
#endif
() {  
	dnpc_finalize();
}
