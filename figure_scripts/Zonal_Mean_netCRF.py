#!/usr/bin/env python3

#
# Script for plotting figure 3.24. for altercation, make necessary comments and uncomments
# of particular axes settings in the plotting section.
#

import os
import sys
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from decimal import Decimal

figure_path = '/home/brandonsmith/climate_state_thesis/figures/'
casenames = ['0.9','0.95','1.0','1.05']

field = 'SWCF'
field2 = 'LWCF'
normalize = True
for CASENAME in casenames:
    run = 'Control'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_control = '/ZonalMean_CRF_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_control):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo -timmean -zonmean -seltimestep,-360/-1 -select,name='+field+','+field2+' '+inpath+infile +' '+ outpath+outfile_control
            print(syscall)
    run = '2xCO2'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'_.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_2xCO2 = '/ZonalMean_CRF_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_2xCO2):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo -timmean -zonmean -seltimestep,-360/-1 -select,name='+field+','+field2+' '+inpath+infile +' '+ outpath+outfile_2xCO2
            print(syscall)
    run = '4xCO2'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'_.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_4xCO2 = '/ZonalMean_CRF_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_4xCO2):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo -timmean -zonmean -seltimestep,-360/-1 -select,name='+field+','+field2+' '+inpath+infile +' '+ outpath+outfile_4xCO2
            print(syscall)
# ONCE LINES ARE PRINTED, COPY OVER TO A SHELL SCRIPT ON A SERVER WHERE CLIMATE DATA OPERATORS IS INSTALLED AND RUN TO CREATE THE FILES.
# ALTERNATIVELY, REPLACE print() STATEMENTS WITH os.system() FUNCTION CALL TO CALL DIRECTLY TO COMMAND LINE.

# Load variables and perform calculations from outfiles
i = 0
fig = plt.figure(figsize=(10,10))
outfile_1C = '/home/brandonsmith/modeloutput/Control/1.0/ZonalMean_CRF_Control_1.0.nc'
outfile_1D = '/home/brandonsmith/modeloutput/2xCO2/1.0/ZonalMean_CRF_2xCO2_1.0.nc'
outfile_1Q = '/home/brandonsmith/modeloutput/4xCO2/1.0/ZonalMean_CRF_4xCO2_1.0.nc'
colors = ['b','g','k','r']
rf = [5.53394852, 4.72283731, 3.88779531, 2.78457941] #defined radiative forcing from doubling CO2 (Byrne and Goldblatt (2014))
rfq = [11.4426504, 9.82042797, 8.15034396, 5.94391217] #defined radiative forcing from quadrupling CO2 (Byrne and Goldblatt (2014))
for CASENAME in casenames:
    outpath_control = '/home/brandonsmith/modeloutput/Control/'+CASENAME
    outpath_2xCO2 = '/home/brandonsmith/modeloutput/2xCO2/'+CASENAME
    outpath_4xCO2 = '/home/brandonsmith/modeloutput/4xCO2/'+CASENAME
    outfile_control = '/ZonalMean_CRF_Control_'+CASENAME+'.nc'
    outfile_2xCO2 = '/ZonalMean_CRF_2xCO2_'+CASENAME+'.nc'
    outfile_4xCO2 = '/ZonalMean_CRF_4xCO2_'+CASENAME+'.nc'
    dsloc_control = outpath_control+outfile_control
    dsloc_2xCO2 = outpath_2xCO2+outfile_2xCO2
    dsloc_4xCO2 = outpath_4xCO2+outfile_4xCO2

    if os.path.isfile(dsloc_control) and os.path.isfile(dsloc_2xCO2):
        # Load Variables
        dsc = netCDF4.Dataset(dsloc_control)
        SWCF = dsc.variables[field][:]
        LWCF = dsc.variables[field2][:]
        lat = dsc.variables['lat'][:]
        dsc.close()

        dsd = netCDF4.Dataset(dsloc_2xCO2)
        SWCFd = dsd.variables[field][:]
        LWCFd = dsd.variables[field2][:]
        # no need to reload lat/lon
        dsd.close()

        dsq = netCDF4.Dataset(dsloc_4xCO2)
        SWCFq = dsq.variables[field][:]
        LWCFq = dsq.variables[field2][:]
        # no need to reload lat/lon
        dsq.close()

        C1 = netCDF4.Dataset(outfile_1C)
        SWCF1 = C1.variables[field][:]
        LWCF1 = C1.variables[field2][:]
        C1.close()
        D1 = netCDF4.Dataset(outfile_1D)
        SWCF1d = D1.variables[field][:]
        LWCF1d = D1.variables[field2][:]
        D1.close()
        Q1 = netCDF4.Dataset(outfile_1Q)
        SWCF1q = Q1.variables[field][:]
        LWCF1q = Q1.variables[field2][:]
        Q1.close()
        
        # Because SWCF directly responds to changes in solar luminosity, a correction must be made by dividing each value for SWCF by S/S_0
        SWCF = SWCF.squeeze()/float(CASENAME)
        SWCFd = SWCFd.squeeze()/float(CASENAME)
        SWCFq = SWCFq.squeeze()/float(CASENAME)
        SWCF1 = SWCF1.squeeze()/float(CASENAME)
        SWCF1d = SWCF1d.squeeze()/float(CASENAME)
        SWCF1q = SWCF1q.squeeze()/float(CASENAME)
        LWCF = LWCF.squeeze()
        LWCFd = LWCFd.squeeze()
        LWCFq = LWCFq.squeeze()
        LWCF1 = LWCF1.squeeze()
        LWCF1d = LWCF1d.squeeze()
        LWCF1q = LWCF1q.squeeze()

        net_CRF = SWCF + LWCF
        net_CRFd = SWCFd + LWCFd
        net_CRFq = SWCFq + LWCFq
        net_CRF1 = SWCF1 + LWCF1
        net_CRF1d = SWCF1d + LWCF1d
        net_CRF1q = SWCF1q + LWCF1q
    # Normalize responses with respect to their radiative forcings. If you wish to see raw doubling responses, set "normalize" to False.
    if normalize is True:
        doubling_response = rf[2]*(net_CRFd-net_CRF)/rf[i]
        diff = doubling_response - (net_CRF1d-net_CRF1)
        Quad_response = rfq[2]*(net_CRFq-net_CRF)/rfq[i]
        Qdiff = Quad_response - (net_CRF1q-net_CRF1)
    else:
        doubling_response = net_CRFd - net_CRF
        diff = doubling_response - (net_CRF1d-net_CRF1) 
        Quad_response = net_CRFq - net_CRF
        Qdiff = Quad_response - (net_CRF1q-net_CRF1)

    Casenames = ['90% $S_0$','95% $S_0$','$S_0$','105% $S_0$']

    lat = np.sin(((lat*np.pi)/180))
    ticks = np.linspace(-1,1,9)
    ticklabels = np.around((180*np.arcsin(ticks)/np.pi),2)
    ticklabels = ticklabels.astype(str)
    if i == 0:
        ax = fig.add_subplot(5,1,1)
        ax.axhline(0,-1,1,color='k',linestyle='dashed',linewidth=0.5)
    ax.plot(lat, net_CRF, colors[i], label=Casenames[i], linewidth=2)
    yticks = np.linspace(-50,0,6).astype(int)
    ydiffticks = np.linspace(-20,20,5).astype(int)


    #ax.ticklabel_format(axis='y',style='sci')
    ax.set_xlim(-1,1)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks,fontsize=14)

    if i == 0:
        ax.set_title('Control Climates',fontsize=16)
    ax.set_ylabel('CRF ($W/m^2$)',fontsize=16)
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticklabels,fontsize=14)
    #ax.legend(loc='upper left')
    ax.grid(True)

    if i == 0:
        ax2 = fig.add_subplot(5, 1, 2)
        ax2.axhline(0,-1,1,color='k',linestyle='dashed',linewidth=0.5)
    ax2.plot(lat, doubling_response, colors[i], label=Casenames[i], linewidth=2)

    ax2.set_xlim(-1,1)
    ax2.set_yticks(ydiffticks)
    ax2.set_yticklabels(ydiffticks,fontsize=14)
    ax2.set_ylabel('$\Delta$ CRF ($W/m^2$)',fontsize=16)
    if i == 0:
        if normalize is False:
            plt.title('Response to $CO_2$ Doubling',fontsize=16)
        else:
            plt.title('Response to $F_2$$_x$',fontsize=16)
    ax2.set_xticks(ticks)
    ax2.set_xticklabels(ticklabels,fontsize=14)
    #ax2.legend(loc='upper left')
    ax2.grid(True)

    if i == 0:
        ax3 = fig.add_subplot(5, 1, 3)
        ax3.axhline(0,-1,1,color='k',linestyle='dashed',linewidth=0.5)
    ax3.plot(lat, diff, colors[i], label=Casenames[i], linewidth=2)

    ax3.set_xlim(-1,1)
    ax3.set_yticks(ydiffticks)
    ax3.set_yticklabels(ydiffticks,fontsize=14)
    ax3.set_ylabel('diff ($W/m^2$)',fontsize=16)
    if i == 0:
        ax3.set_title('Difference W.R.T. $S_0$',fontsize=16)
    ax3.set_xticks(ticks)
    ax3.set_xticklabels(ticklabels,fontsize=14)
    #ax3.legend(loc='upper left')
    ax3.grid(True)

    if i == 0:
        ax4 = fig.add_subplot(5, 1, 4)
        ax4.axhline(0,-1,1,color='k',linestyle='dashed',linewidth=0.5)
    ax4.plot(lat, Quad_response, colors[i], label=Casenames[i], linewidth=2)

    ax4.set_xlim(-1,1)
    ax4.set_yticks(ydiffticks)
    ax4.set_yticklabels(ydiffticks,fontsize=14)
    ax4.set_ylabel('$\Delta$ CRF ($W/m^2$)',fontsize=16)
    if i == 0:
        if normalize is False:
            plt.title('Response to $CO_2$ Quadrupling',fontsize=16)
        else:
            plt.title('Response to $F_4$$_x$',fontsize=16)
    ax4.set_xticks(ticks)
    ax4.set_xticklabels(ticklabels,fontsize=14)
    #ax2.legend(loc='upper left')
    ax4.grid(True)

    if i == 0:
        ax5 = fig.add_subplot(5, 1, 5)
        ax5.axhline(0,-1,1,color='k',linestyle='dashed',linewidth=0.5)
    ax5.plot(lat, Qdiff, colors[i], label=Casenames[i], linewidth=2)

    ax5.set_xlim(-1,1)
    ax5.set_yticks(ydiffticks)
    ax5.set_yticklabels(ydiffticks,fontsize=14)
    ax5.set_ylabel('diff ($W/m^2$)',fontsize=16)
    if i == 0:
        ax5.set_title('Difference W.R.T. $S_0$',fontsize=16)
    if i == 3:
        ax5.set_xlabel('Latitude (Deg N)',fontsize=16)
    ax5.set_xticks(ticks)
    ax5.set_xticklabels(ticklabels,fontsize=14)
    #ax3.legend(loc='upper left')
    ax5.grid(True)

    i +=1

if normalize is False:
    fig.suptitle('Zonal Mean Net Cloud Forcing Response to $CO_2$ Doubling',fontsize=24)
else:
    fig.suptitle('Zonal Mean Net Cloud Forcing Response to Equivalent Radiative Forcing',fontsize=18,y=1.05)

handles, labels = ax.get_legend_handles_labels()
fig.legend(handles,labels,loc = (0.3, 0.93), ncol=4 )
fig.tight_layout(pad=0.2)
plt.show()
if normalize is False:
    fig.savefig(figure_path+'ZM_CRF_2x_Response.pdf',bbox_inches='tight')
else:
    fig.savefig(figure_path+'ZM_CRF_n_Response.pdf',bbox_inches='tight')
###################################################################################################################################################


