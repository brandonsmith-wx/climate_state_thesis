** SSH into the Compute Canada Cluster:
ssh -X user@cedar.computecanada.ca

** Helpful links: 
https://docs.computecanada.ca/wiki/Cedar
http://www.cesm.ucar.edu/models/cesm1.2/cesm/doc/usersguide/book1.html

**You'll want to clone whatever repo you're using onto your local workspace on cedar before you begin

** Create new case **
 Default case is setup to use CAM5 with slab ocean at f19+g16 resolution. Takes a casename as a parameter

source baserun_cesm.sh testcase

** Copy configuration file to $CASEROOT directory **
 Make any edits to pe layout or input files here after copying

cp configure_cesm.sh ../cesmruns/testcase/
cd ../cesmruns/testcase/
source configure_cesm.sh

** Lock in configuration **
 When happy with configuration, run cesm_setup to create buildfiles and namelists with options chosen.
** IMPORTANT: You cannot change the configuration xml files (except env_run.xml) after running ./cesm_setup **

./cesm_setup

** Change namelist values **
 Here is where you set things like CO2 concentration, solar constant, etc. Edit user_nl_cam:
vi user_nl_cam

and add entries such as:
 co2vmr = 200e-6

to set the CO2 concentration to 200 ppm. All available nameless variables are listed here:
http://www.cesm.ucar.edu/cgi-bin/eaton/namelist/nldef2html-cam5_3
To set the solar constant, run the script change_solar_cont.py with the multiple of 1366 that you want to set the solar constant to, from the CASEROOT directory:

cp gcmexpt/change_solar_cont.py cesmruns/testcase
./change_solar_cont.py 0.9 	


** this creates a solar irradiance file that sets the solar constant to 90% of the default value 1366, ie 1229, and adds an entry to the namelist to use this solar data file instead of the default one. Check user_nl_cam after running this script to make sure it looks right, it's possible that it will have overwritten a previous line of the file such as the CO2 value you set previously.

** Build the case **
./testcase.build

** Change walltime, set up email notifications **
 Edit testcase.run:
vi testcase.run

(for instructions on using Vim as an editor in the terminal see https://vimhelp.appspot.com/vim_faq.txt.html)

** Due to a bug (https://bb.cgd.ucar.edu/cesm122-problem-creating-timing-files-when-doutstrue), 
 edit testcase.run and in the "Perform short term archiving of output" section, replace:

  cd $RUNDIR; $CASETOOLS/st_archive.sh
with
  cd $RUNDIR; $CASETOOLS/st_archive.sh; cd $CASEROOT


** Run the case **
** Execute the following command to add the run to the queue:

sbatch --account=def-czg testcase.run

Check run status with: squeue -u $USER
Cedar uses Slurm as the queue manager. More details on how to run cases here: 
https://docs.computecanada.ca/wiki/Running_jobs

**REMINDER: Once you have a successful first run, you must set CONTINUE_RUN to TRUE in env_run.xml before resubmitting, otherwise the job will not progress. 
 You may also need to modify the RESUBMIT, STOP_OPTION, STOP_N, STOP_DATE, REST_OPTION, REST_N and/or REST_DATE variables in env_run.xml before resubmitting. 

** Details of the run such as whether or not it was successful will be written to /cesmruns/testcase/slurm.out. From here it may refer to you other error logs if necessary for troubleshooting. 

** Actual model output will be written to ~/scratch/ccsm/archive/testcase
** CAM output located at ~/scratch/ccsm/archive/testcase/atm/hist
