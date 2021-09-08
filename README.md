# DNPC: A Dynamic Node-Level Power Capping Library for HPC Applications

## Overview
DNPC is a user-level dynamic power capping library for dynamic power and energy management. 
It works by combining the energy monitoring and hardware counter monitoring of PoLiMEr and PAPI, respectively.
It employs one of the built-in or user-supplied dynamic power capping methods during run-time to set power-caps on the application.
DNPC works with C and C++ applications through simply including the library and calling the appropriate DNPC functions.
It also exposes a set of Fortran wrapper functions that allow DNPC to work with Fortran applications as well.

While DNPC is a user-level library, PoLiMEr internally uses RAPL for power capping, so the user should contact their system administrator for privileges to read/write certain MSRs as outlined in PoLiMEr's documentation.

## Using DNPC
### Dependencies
- PAPI
	- Terpstra, D., Jagode, H., You, H., Dongarra, J. "Collecting Performance Data with PAPI-C", Tools for High Performance Computing 2009, Springer Berlin / Heidelberg, 3rd Parallel Tools Workshop, Dresden, Germany, pp. 157-173, 2010. 
- PoLiMEr (Modified Version Required)
	- Modified Repository: https://github.com/SPEAR-IIT/PoLiMEr-dnpc
	- Original Repository: https://xgitlab.cels.anl.gov/open-source/PoLiMEr
	- Ivana Marincic, Venkatram Vishwanath, and Henry Hoffmann. 2017. PoLiMEr: An Energy Monitoring and Power Limiting Interface for HPC Applications. In Proceedings of the 5th International Workshop on Energy Efficient Supercomputing (E2SC'17). Association for Computing Machinery, New York, NY, USA, Article 7, 1–8. DOI:https://doi.org/10.1145/3149412.3149419

### Adding DNPC To Applications 
A minimal working example of a C/C++ application instrumented with DNPC is shown below.
```C
#include "dnpc.h"

int
main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	dnpc_init(algorithm);

	// ...
	// app code here
	// ...

	dnpc_finalize();
	MPI_Finalize();
	return 0;
}
```
If MPI is not being used, then the MPI statements can be ommitted.

## Configuring DNPC
DNPC uses environment variables to control its behavior.
Before running, it is important to set the value of these variables.
- Controlling Capping Behavior
	- `PERF_DEG_THRESHOLD`
		- Value should be a number within 0 and 100
		- Represents the percentage of degradation the algorithm should allow
	- `ALGORITHM`
		- Value should be a number corresponding to an algorithm
		- Represents the algorithm used to decide capping behavior
		- The default supported algorithms included are the following:
			- Value: `-1`, measures power and performance, but performs no capping 
			- Value: `100`, Sets a single, static power cap throughout execution 
			- Value: `6`, Checks change in IPC, with performance checking
			- Value: `7`, Checks change in Micro-Operations, with performance checking
	- `STATIC_POWER_CAP`
		- Value should be a number within minimum and maximum power cap values
		- When the appropriate algorithm is set, this variable is used to decide the value of the static power cap to be used throughout the entire execution of the application
	- `PAPI_EVENTS`
		- Value should be a comma-separated list of PAPI event names
		- ex: `PAPI_EVENTS="PAPI_TOT_CYC,PAPI_TOT_INS,LLC_MISSES,UOPS_RETIRED"`
			- Note: this configuration of PAPI events is necessary for all the default algorithms to work correctly
		- Represents the list of PAPI counters to measure during execution. Some algorithms may require certain PAPI counters to work correctly.
- Tuning Capping Behavior
	- `FREQ_GIVE_OR_TAKE`
	- `IPC_GIVE_OR_TAKE`
	- `UOPS_GIVE_OR_TAKE`
	- `POWER_STEP`
	- `MAX_VALUES`
	- `SAMPLE_SIZE`
	- `MY_POLL_INTERVAL`
- Bookkeeping
	- `EXPERIMENT`
	- `CONFIGURATION`
	- `TRIAL`
	- `DEBUG_FLAG`

## Linking to DNPC
Include the `dnpc.h` file in your applications.
When compiling the application code, use dynamic linking, as shown below:

```bash
gcc -I<DNPC_PATH>/include -Wl,-rpath=<DNPC_PATH>/lib <APP>.c -o <APP> -L<DNPC_PATH>/lib -ldnpc
```

With `<DNPC_PATH>` set to the location of the DNPC repository on your system.

## Building DNPC
Clone the repository. 
Afterwards, edit the Makefile to include the appropriate include and library paths to PAPI and PoLiMEr on your system. 
Make sure that PoLiMEr was compiled with the `TIMER_OFF` flag in order to not conflict with DNPC's timer.
Once those are complete, run the following commands from the directory's root.
```
mkdir bin; mkdir lib
make
```

## Citing DNPC 
Note: If you use DNPC in your work, please cite the following work:

Sahil Sharma, Zhiling Lan, Xingfu Wu, and Valerie Taylor, “A Dynamic Power Capping Library for HPC Applications”, IEEE Cluster 2021 (research poster)
[Cluster21Poster.pdf](https://github.com/SPEAR-IIT/DNPC/files/7131672/Cluster21Poster.pdf)]
