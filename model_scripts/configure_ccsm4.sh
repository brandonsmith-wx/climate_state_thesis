#!//bin/bash
# Run 20 years, 1 year at a time. Run as-is for initial case, then uncomment CONTINUE_RUN,
# change RESUBMIT to 19 run for 20 years total.
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
#./xmlchange -file env_mach_pes.xml -id NTASKS_ROF -val 16
#./xmlchange -file env_mach_pes.xml -id NTHRDS_ROF -val 1
#./xmlchange -file env_mach_pes.xml -id ROOTPE_ROF -val 0
#./xmlchange -file env_mach_pes.xml -id NTASKS_WAV -val 16
#./xmlchange -file env_mach_pes.xml -id NTHRDS_WAV -val 1
#./xmlchange -file env_mach_pes.xml -id ROOTPE_WAV -val 0
# short term archiving
./xmlchange -file env_run.xml -id DOUT_S -val TRUE
# length of time for model to run for
./xmlchange -file env_run.xml -id STOP_OPTION -val nyears
./xmlchange -file env_run.xml -id STOP_N -val 1
./xmlchange -file env_run.xml -id REST_OPTION -val nyears
./xmlchange -file env_run.xml -id REST_N -val 1
# how many restarts to do
./xmlchange -file env_run.xml -id RESUBMIT -val 0
#./xmlchange -file env_run.xml -id CONTINUE_RUN -val TRUE
# SET SOM file - this SOM file is only for testing purposes, see https://bb.cgd.ucar.edu/faq-data-ocean-slab-mode-docn-som
#./xmlchange -file env_run.xml -id DOCN_SOM_FILENAME -val pop_frc.1x1d.090130.nc
