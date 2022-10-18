** SSH into the Compute Canada Cluster:
ssh -X user@cedar.computecanada.ca
 
** Helpful links: 
https://docs.computecanada.ca/wiki/Cedar
http://www.cesm.ucar.edu/models/ccsm4.0/ccsm_doc/book1.html
 
** STEP 1. Create new case **
- Default case is setup to use CAM5 with slab ocean at f19+g16 resolution. Takes a casename as a parameter

- A simple example is setup with: 
source ./baserun_ccsm4.sh baseruntest#

 (The source command is important as it sets variable values in the parent shell of the script.
 This step creates the new $CASEROOT directory as: $CASEROOT=/home/username/cesmruns/baseruntest# 
 where baseruntest# is the case directory. Replace the # with the case number, e.g. baseruntest1)

** STEP 2. Modify configuration for the new case
** IMPORTANT: make configuration changes BEFORE running ./configure -case. 
 Here is where we use xmlchange to set things like DOUT_S, STOP_OPTION, etc. 
 Things that set how long the simulation runs, how often restart files are output, and where model output goes.
 Change the options in env_run.xml by editing them in configure_example.sh and running the script as below. 
 [Note:Use chmod u+x to have execute permissions on configure_example.sh]:

source configure_ccsm4.sh

(Run this script configure_example.sh from $CASEROOT directory. 
Make sure the variables [STOP_N and REST_N] for env_run.xml written in configure_example.sh are reflected as required before the case is run. 
The variables RESUBMIT and CONTINUE_RUN=TRUE need to be set only if re-running the case.)

** STEP 3. Lock in configuration for the case
** IMPORTANT: You cannot change the configuration xml files (except env_run.xml) after running ./configure -case **
 
cd $CASEROOT 
./configure -case

 The baseruntest#.build script should be created in $CASEROOT directory.

 Make sure walltime is set to some value (max of 72 hours) in baseruntest#.run script 
 (This command can be added at the top of the script):
 
#PBS -l walltime=04:00:00

** Changing the solar const and co2:

 Edit the Buildconf/cam.buildnml.csh in the $CASEROOT directory, 
 to set solar_const and co2vmr variables according to values in the table provided.


** STEP 4. Build the case

cd $CASEROOT
./baseruntest#.nestor.build

- Downloading data.
 If the model needs to download data, it should do this mostly automatically. The
 first time, one needs to make the directory that it wants to put the input data in. e.g. 
mkdir /global/scratch/czg/ccsm/inputdata
 Then, then it should download input data automatically. The first time this is
 done (or first few, becuase there are a few sources), one needs to setup access
 to the subversion repository. 
 First, edit 
~/.subversion/servers 
 to set 
store-passwords = yes
 Then, build will ask for username and password
un=guestuser
pw=friendly
 These will then be cached and not need to be entered repeatedly from the same machine. 


** STEP 5. Run the case
**Execute the following command to add the run to the queue:

sbatch --account=def-czg baseruntest#.run

Check run status with: squeue -u $USER
Cedar uses Slurm as the queue manager. More details on how to run cases here: 
https://docs.computecanada.ca/wiki/Running_jobs

**REMINDER: Once you have a successful first run, you must set CONTINUE_RUN to TRUE in env_run.xml before resubmitting, otherwise the job will not progress. 
 You may also need to modify the RESUBMIT, STOP_OPTION, STOP_N, STOP_DATE, REST_OPTION, REST_N and/or REST_DATE variables in env_run.xml before resubmitting. 
