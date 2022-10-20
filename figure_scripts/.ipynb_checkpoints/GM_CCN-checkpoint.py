#!/usr/bin/env python3
###################################################################################################################################################

#
# Script for plotting figure 3.21. for altercation, make necessary comments and uncomments
# of particular axes settings in the plotting section.
#

import os
import sys
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import time
from matplotlib.colors import Normalize
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import time

figure_path = '/home/brandonsmith/climate_state_thesis/figures/'
casenames = ['0.9','0.95','1.0','1.05']

for CASENAME in casenames:
    run = 'Control'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME # change to location of parent data set
    infile = '/Merged_'+run+'_'+CASENAME+'.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_control = '/CCNs_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_control):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo fldmean -timmean -seltimestep,-360/-1 -select,name=CCN3,T '+ inpath+infile+' '+ outpath+outfile_control
            print(syscall)
    run = '2xCO2'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'_.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_control = '/CCNs_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_control):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo fldmean -timmean -seltimestep,-360/-1 -select,name=CCN3,T '+ inpath+infile+' '+ outpath+outfile_control
            print(syscall)
    run = '4xCO2'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'_.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_control = '/CCNs_levels_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_control):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo fldmean -timmean -seltimestep,-360/-1 -select,name=CCN3,T '+ inpath+infile+' '+ outpath+outfile_control
            print(syscall)
            
# ONCE LINES ARE PRINTED, COPY OVER TO A SHELL SCRIPT ON A SERVER WHERE CLIMATE DATA OPERATORS IS INSTALLED AND RUN TO CREATE THE FILES.
# ALTERNATIVELY, REPLACE print() STATEMENTS WITH os.system() FUNCTION CALL TO CALL DIRECTLY TO COMMAND LINE.

# Load variables and perform calculations from outfiles
i = 0
fig = plt.figure()
outfile_1C = '/home/brandonsmith/modeloutput/Control/1.0/CCNs_Control_1.0.nc'
outfile_1D = '/home/brandonsmith/modeloutput/2xCO2/1.0/CCNs_2xCO2_1.0.nc'
outfile_1Q = '/home/brandonsmith/modeloutput/4xCO2/1.0/CCNs_levels_4xCO2_1.0.nc'
rf = [5.53394852, 4.72283731, 3.88779531, 2.78457941] #defined radiative forcing from doubling CO2 (Byrne and Goldblatt (2014))
rfq = [11.4426504, 9.82042797, 8.15034396, 5.94391217] #defined radiative forcing from quadrupling CO2 (Byrne and Goldblatt (2014))
column_CCN = []
column_CCNd = []
column_CCNq = []
for CASENAME in casenames:
    outfile_2xco2 = '/CCNs_2xCO2_'+CASENAME+'.nc'
    outfile_control = '/CCNs_Control_'+CASENAME+'.nc'
    outfile_4xco2 = '/CCNs_levels_4xCO2_'+CASENAME+'.nc'
    outpath_control = '/home/brandonsmith/modeloutput/Control/'+CASENAME
    outpath_2xco2 = '/home/brandonsmith/modeloutput/2xCO2/'+CASENAME
    outpath_4xco2 = '/home/brandonsmith/modeloutput/4xCO2/'+CASENAME
    dsloc_control = outpath_control+outfile_control
    dsloc_2xco2 = outpath_2xco2+outfile_2xco2
    dsloc_4xco2 = outpath_4xco2+outfile_4xco2
    
    if os.path.isfile(dsloc_control) and os.path.isfile(dsloc_2xco2):
        # Load Variables
        dsc = netCDF4.Dataset(dsloc_control)
        CCN = dsc.variables['CCN3'][:]
        T = dsc.variables['T'][:]
        lev = dsc.variables['lev'][:]
        dsc.close()
        
        dsd = netCDF4.Dataset(dsloc_2xco2)
        CCNd = dsd.variables['CCN3'][:]
        Td = dsd.variables['T'][:]
        dsd.close()
        
        dsq = netCDF4.Dataset(dsloc_4xco2)
        CCNq = dsq.variables['CCN3'][:]
        Tq = dsq.variables['T'][:]
        dsq.close()

        C1 = netCDF4.Dataset(outfile_1C)
        CCN1 = C1.variables['CCN3'][:]
        T1 = C1.variables['T'][:]
        C1.close()
        
        D1 = netCDF4.Dataset(outfile_1D)
        CCN1d = D1.variables['CCN3'][:]
        T1d = D1.variables['T'][:]
        D1.close()
        
        Q1 = netCDF4.Dataset(outfile_1Q)
        CCN1q = Q1.variables['CCN3'][:]
        T1q = Q1.variables['T'][:]
        Q1.close()
        
        CCN = CCN*1e6 #convert to #/m^3
        CCNd = CCNd*1e6
        CCNq = CCNq*1e6
        
        T = np.fliplr(T).squeeze()
        Td = np.fliplr(Td).squeeze()
        Tq = np.fliplr(Tq).squeeze()
        T1 = np.fliplr(T1).squeeze()
        T1d = np.fliplr(T1d).squeeze()
        T1q = np.fliplr(T1q).squeeze()
        
        CCN = np.fliplr(CCN).squeeze()
        CCNd = np.fliplr(CCNd).squeeze()
        CCNq = np.fliplr(CCNq).squeeze()
        CCN1 = np.fliplr(CCN1).squeeze()
        CCN1d = np.fliplr(CCN1d).squeeze()
        CCN1q = np.fliplr(CCN1q).squeeze()
        lev = np.flipud(lev).squeeze()
                
        R = 287
        g = 9.8
        P0 = 100000
        viccn = []
        viccnd = []
        viccnq = []
        viccn1 = []
        viccn1d = []
        viccn1q = []
        
        for j in range(0,len(lev)):
            if j < 29:
                ccn = (-R/g)*CCN[j]*( (T[j+1]*np.log(lev[j+1]/P0)) - (T[j]*np.log(lev[j]/P0)) )
                ccnd = (-R/g)*CCNd[j]*( (Td[j+1]*np.log(lev[j+1]/P0)) - (Td[j]*np.log(lev[j]/P0)) )
                ccnq = (-R/g)*CCNq[j]*( (Tq[j+1]*np.log(lev[j+1]/P0)) - (Tq[j]*np.log(lev[j]/P0)) )
                ccn1 = (-R/g)*CCN1[j]*( (T1[j+1]*np.log(lev[j+1]/P0)) - (T1[j]*np.log(lev[j]/P0)) )
                ccn1d = (-R/g)*CCN1d[j]*( (T1d[j+1]*np.log(lev[j+1]/P0)) - (T1d[j]*np.log(lev[j]/P0)) )
                ccn1q = (-R/g)*CCN1q[j]*( (T1q[j+1]*np.log(lev[j+1]/P0)) - (T1q[j]*np.log(lev[j]/P0)) )
                
                viccn.append(ccn)
                viccnd.append(ccnd)
                viccnq.append(ccnq)
                viccn1.append(ccn1)
                viccn1d.append(ccnd)
                viccn1q.append(ccnq)
        
        column_integrated_CCN = np.sum(viccn)
        column_integrated_CCNd = np.sum(viccnd)
        column_integrated_CCNq = np.sum(viccnq)
        column_integrated_CCN1 = np.sum(viccn1)
        column_integrated_CCN1d = np.sum(viccn1d)
        column_integrated_CCN1q = np.sum(viccn1q)
        
    column_CCN.append(column_integrated_CCN)
    column_CCNd.append(column_integrated_CCNd)
    column_CCNq.append(column_integrated_CCNq)
        
    i+=1
    

ax = fig.add_axes([0,0,1,1])
ax.plot(casenames,column_CCN,'o',color='k',linestyle='solid')
ax.plot(casenames,column_CCNd,'o',color='b',linestyle='solid')
ax.plot(casenames,column_CCNq,'o',color='r',linestyle='solid')
ax.set_title('CCN column number')
ax.set_ylabel('concentration ($m^{-2}$)')
ax.set_xlabel('Solar Case')
ax.set_xticklabels(['90%','95%','100%','105%'])
plt.grid(True)
plt.legend(['Base Climate','2x$CO_2$','4x$CO_2$'])
plt.show()

fig.savefig(figure_path+'vertically_integrated_CCN_response.pdf',bbox_inches='tight')


###################################################################################################################################################