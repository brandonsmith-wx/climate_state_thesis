#!//bin/bash
# example file for CESM setup - testing 2 day run with 16 processors
# processor/thread configuration
./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val 16
./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val 1
./xmlchange -file env_mach_pes.xml -id ROOTPE_ATM -val 0
./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val 16
./xmlchange -file env_mach_pes.xml -id NTHRDS_LND -val 1
./xmlchange -file env_mach_pes.xml -id ROOTPE_LND -val 0
./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val 16
./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val 1
./xmlchange -file env_mach_pes.xml -id ROOTPE_ICE -val 0
./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val 16
./xmlchange -file env_mach_pes.xml -id NTHRDS_OCN -val 1
./xmlchange -file env_mach_pes.xml -id ROOTPE_OCN -val 0
./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val 16
./xmlchange -file env_mach_pes.xml -id NTHRDS_CPL -val 1
./xmlchange -file env_mach_pes.xml -id ROOTPE_CPL -val 0
./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val 16
./xmlchange -file env_mach_pes.xml -id NTHRDS_GLC -val 1
./xmlchange -file env_mach_pes.xml -id ROOTPE_GLC -val 0
# short term archiving
./xmlchange -file env_run.xml -id DOUT_S -val FALSE
# length of time for model to run for
./xmlchange -file env_run.xml -id STOP_OPTION -val ndays
./xmlchange -file env_run.xml -id STOP_N -val 2
./xmlchange -file env_run.xml -id REST_OPTION -val ndays
./xmlchange -file env_run.xml -id REST_N -val 2
# how many restarts to do
./xmlchange -file env_run.xml -id RESUBMIT -val 0
#./xmlchange -file env_run.xml -id CONTINUE_RUN -val TRUE
