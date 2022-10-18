#!//bin/bash
# Run 20 years, 1 year at a time. Run as-is for initial case, then uncomment CONTINUE_RUN,
# change RESUBMIT to 19 un for 20 years total.
# processor/thread configuration
./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val 288
./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val 1
./xmlchange -file env_mach_pes.xml -id ROOTPE_ATM -val 0
./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val 48
./xmlchange -file env_mach_pes.xml -id NTHRDS_LND -val 1
./xmlchange -file env_mach_pes.xml -id ROOTPE_LND -val 288
./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val 48
./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val 1
./xmlchange -file env_mach_pes.xml -id ROOTPE_ICE -val 336
./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val 48
./xmlchange -file env_mach_pes.xml -id NTHRDS_OCN -val 1
./xmlchange -file env_mach_pes.xml -id ROOTPE_OCN -val 384
./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val 432
./xmlchange -file env_mach_pes.xml -id NTHRDS_CPL -val 1
./xmlchange -file env_mach_pes.xml -id ROOTPE_CPL -val 0
./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val 48
./xmlchange -file env_mach_pes.xml -id NTHRDS_GLC -val 1
./xmlchange -file env_mach_pes.xml -id ROOTPE_GLC -val 0
./xmlchange -file env_mach_pes.xml -id NTASKS_ROF -val 48
./xmlchange -file env_mach_pes.xml -id NTHRDS_ROF -val 1
./xmlchange -file env_mach_pes.xml -id ROOTPE_ROF -val 0
./xmlchange -file env_mach_pes.xml -id NTASKS_WAV -val 48
./xmlchange -file env_mach_pes.xml -id NTHRDS_WAV -val 1
./xmlchange -file env_mach_pes.xml -id ROOTPE_WAV -val 0
./xmlchange -file env_mach_pes.xml -id TOTALPES -val 432
./xmlchange -file env_mach_pes.xml -id MAX_TASKS_PER_NODE -val 48
# short term archiving
./xmlchange -file env_run.xml -id DOUT_S -val TRUE
# Set reference case and date
./xmlchange -file env_run.xml -id RUN_TYPE -val branch
./xmlchange -file env_run.xml -id RUN_REFCASE -val controlrun
./xmlchange -file env_run.xml -id RUN_REFDATE -val 0031-01-01
# length of time for model to run for
./xmlchange -file env_run.xml -id STOP_OPTION -val nyears
./xmlchange -file env_run.xml -id STOP_N -val 5
./xmlchange -file env_run.xml -id REST_OPTION -val nyears
./xmlchange -file env_run.xml -id REST_N -val 5
# how many restarts to do
./xmlchange -file env_run.xml -id RESUBMIT -val 0
# Change CO2 concentration
./xmlchange -file env_run.xml -id CCSM_CO2_PPMV -val 
# Change solar constant
#./change_solar_cont.py 1.0
# Adjust batch submit settings
./xmlchange -file env_run.xml -id BATCHSUBMIT -val "sbatch --account=def-czg"
# SET SOM file
./xmlchange -file env_run.xml -id DOCN_SOM_FILENAME -val pop_frc.b.e11.B1850C5CN.f09_g16.005.082914.nc
# set DEBUG to TRUE
#./xmlchange -file env_build.xml -id DEBUG -val TRUE
