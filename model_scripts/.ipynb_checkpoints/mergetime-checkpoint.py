#!/usr/bin/env python3

import os
import sys

#filebase = '/nazko/home/vmcd/projects/modelOutput/climatesensitivity/'+run+'/'
#filebase = '/home/bps5238/scratch/ccsm/archive/'
run = sys.argv[1]
filebase = '/home/brandonsmith/modeloutput/'+run+'/'
outfilebase = 'Merged_'+run+'_'
casenames = {'0.9','0.95','1.0','1.05'}
#merge files together for analysis
for CASENAME in casenames:
	#out = filebase2+outfilebase+CASENAME+'_'+experiment+'.nc'
	out = filebase+CASENAME+'/'+outfilebase+CASENAME+'.nc'
	# check directly if the mergetime file exists
	if not os.path.isfile(out):
#		if os.path.isdir(filebase+run+'/'+CASENAME+'/controlrun/atm/hist/'):
		if os.path.isdir(filebase+run+'_'+CASENAME+'/atm/hist/'):
			filestomerge = filebase+run+'_'+CASENAME+'/atm/hist/'+run+'_*.nc'
			# merge files
			syscall = 'cdo mergetime '+filestomerge+' '+out
			os.system(syscall)
		else:
			print('directory not found: '+filebase+run+'_'+CASENAME+'/atm/hist/')

                        



