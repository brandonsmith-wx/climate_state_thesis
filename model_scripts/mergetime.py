#!/usr/bin/env python3

import os
import sys

run = sys.argv[1]

#filebase = '/nazko/home/vmcd/projects/modelOutput/climatesensitivity/'+run+'/'
infilebase = '/kenes/user/brandonsmith/modeloutput/'+run+'/'
filebase = '/home/brandonsmith/modeloutput/'+run+'/'
outfilebase = 'Merged_'+run+'_'
casenames = {'0.9','0.95','1.0','1.05'}
#merge files together for analysis
for CASENAME in casenames:
	#out = filebase2+outfilebase+CASENAME+'_'+experiment+'.nc'
	out = filebase+CASENAME+'/'+outfilebase+CASENAME+'.nc'
	outshifted = filebase+CASENAME+'/'+outfilebase+CASENAME+'_.nc'
	# check directly if the mergetime file exists
	if not os.path.isfile(out):
		if os.path.isdir(infilebase+run+'_'+CASENAME+'/atm/hist/'):
			filestomerge = infilebase+run+'_'+CASENAME+'/atm/hist/'+run+'_*.nc'
			# merge files												     
			syscall = 'cdo mergetime '+filestomerge+' '+out
			os.system(syscall)
		else:
			print('directory not found: '+infilebase+run+'_'+CASENAME+'/atm/hist/')
	if not os.path.isfile(outshifted):
		syscall = 'cdo shifttime,-1mo '+out+' '+outshifted
		os.system(syscall)
	else:
		print('directory not found: '+filebase+CASENAME)

			



