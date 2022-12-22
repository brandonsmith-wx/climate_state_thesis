#!/usr/bin/env python3

#
# Script for plotting figure 3.15. for altercation, make necessary comments and uncomments
# of particular axes settings in the plotting section.
#

import os
import sys
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.colors as colors

figure_path = '/home/brandonsmith/climate_state_thesis/figures/'
casenames = ['0.9','0.95','1.0','1.05']

field = 'SWCF'
normalize = True
for CASENAME in casenames:
    run = 'Control'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_control = '/ZonalMean_SWCF_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_control):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo -timmean -zonmean -select,name='+field+' '+inpath+infile +' '+ outpath+outfile_control
            print(syscall)
    run = '2xCO2'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'_.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_2xCO2 = '/ZonalMean_SWCF_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_2xCO2):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo -timmean -zonmean -seltimestep,-360/-1 -select,name='+field+' '+inpath+infile +' '+ outpath+outfile_2xCO2
            print(syscall)
    run = '4xCO2'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'_.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_4xCO2 = '/ZonalMean_SWCF_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_4xCO2):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo -timmean -zonmean -seltimestep,-360/-1 -select,name='+field+' '+inpath+infile +' '+ outpath+outfile_4xCO2
            print(syscall)
# ONCE LINES ARE PRINTED, COPY OVER TO A SHELL SCRIPT ON A SERVER WHERE CLIMATE DATA OPERATORS IS INSTALLED AND RUN TO CREATE THE FILES.
# ALTERNATIVELY, REPLACE print() STATEMENTS WITH os.system() FUNCTION CALL TO CALL DIRECTLY TO COMMAND LINE.

# Load variables and perform calculations from outfiles
i = 0
fig = plt.figure(figsize=(5,10))
outfile_1C = '/home/brandonsmith/modeloutput/Control/1.0/ZonalMean_SWCF_Control_1.0.nc'
outfile_1D = '/home/brandonsmith/modeloutput/2xCO2/1.0/ZonalMean_SWCF_2xCO2_1.0.nc'
outfile_1Q = '/home/brandonsmith/modeloutput/4xCO2/1.0/ZonalMean_SWCF_4xCO2_1.0.nc'
colors = ['b','g','k','r']
rf = [5.53394852, 4.72283731, 3.88779531, 2.78457941] #defined radiative forcing from doubling CO2 (Byrne and Goldblatt (2014))
rfq = [11.4426504, 9.82042797, 8.15034396, 5.94391217] #defined radiative forcing from quadrupling CO2 (Byrne and Goldblatt (2014))
for CASENAME in casenames:
    outpath_control = '/home/brandonsmith/modeloutput/Control/'+CASENAME
    outpath_2xCO2 = '/home/brandonsmith/modeloutput/2xCO2/'+CASENAME
    outpath_4xCO2 = '/home/brandonsmith/modeloutput/4xCO2/'+CASENAME
    outfile_control = '/ZonalMean_SWCF_Control_'+CASENAME+'.nc'
    outfile_2xCO2 = '/ZonalMean_SWCF_2xCO2_'+CASENAME+'.nc'
    outfile_4xCO2 = '/ZonalMean_SWCF_4xCO2_'+CASENAME+'.nc'
    dsloc_control = outpath_control+outfile_control
    dsloc_2xCO2 = outpath_2xCO2+outfile_2xCO2
    dsloc_4xCO2 = outpath_4xCO2+outfile_4xCO2

    if os.path.isfile(dsloc_control) and os.path.isfile(dsloc_2xCO2):
        # Load Variables
        dsc = netCDF4.Dataset(dsloc_control)
        Cvar= dsc.variables[field][:]
        lat = dsc.variables['lat'][:]
        dsc.close()

        dsd = netCDF4.Dataset(dsloc_2xCO2)
        Dvar = dsd.variables[field][:]
        # no need to reload lat/lon
        dsd.close()

        dsq = netCDF4.Dataset(dsloc_4xCO2)
        Qvar = dsq.variables[field][:]
        # no need to reload lat/lon
        dsq.close()

        C1 = netCDF4.Dataset(outfile_1C)
        Cvar1 = C1.variables[field][:]
        C1.close()
        D1 = netCDF4.Dataset(outfile_1D)
        Dvar1 = D1.variables[field][:]
        D1.close()
        Q1 = netCDF4.Dataset(outfile_1Q)
        Qvar1 = D1.variables[field][:]
        Q1.close()

        # Because SWCF directly responds to changes in solar luminosity, a correction must be made by dividing each value for SWCF by S/S_0
        Cvar = Cvar/float(CASENAME)
        Cvar = Cvar.squeeze()
        Dvar = Dvar/float(CASENAME)
        Dvar = Dvar.squeeze()
        Qvar = Qvar/float(CASENAME)
        Qvar = Qvar.squeeze()
        Cvar1 = Cvar1/float(CASENAME)
        Cvar1 = Cvar1.squeeze()
        Dvar1 = Dvar1/float(CASENAME)
        Dvar1 = Dvar1.squeeze()
        Qvar1 = Qvar1/float(CASENAME)
        Qvar1 = Qvar1.squeeze()

    # Normalize responses with respect to their radiative forcings. If you wish to see raw doubling responses, set "normalize" to False.
    if normalize is True:
        doubling_response = rf[2]*(Dvar-Cvar)/rf[i]
        diff = doubling_response - (Dvar1-Cvar1)
        Quad_response = rfq[2]*(Qvar-Cvar)/rfq[i]
        Qdiff = Quad_response - (Qvar1-Cvar1)
    else:
        doubling_response = Dvar - Cvar
        diff = doubling_response - (Dvar1-Cvar1) 
        Quad_response = Qvar - Cvar
        Qdiff = Quad_response - (Qvar1-Cvar1)

    Casenames = ['90% $S_0$','95% $S_0$','$S_0$','105% $S_0$']

    lat = np.sin(((lat*np.pi)/180))
    ticks = np.linspace(-90,90,7)
    ticklabels = ticks.astype(int)
    ticks = np.sin(((ticks*np.pi)/180))
    if i == 0:
        ax = fig.add_subplot(5,1,1)
        ax.axhline(0,-1,1,color='k',linestyle='dashed',linewidth=0.5)
    ax.plot(lat, Cvar, colors[i], label=Casenames[i], linewidth=2)
    yticks = np.linspace(-75,0,6).astype(int)
    ydiffticks = np.linspace(-30,30,7).astype(int)
    ax.set_xlim(-1,1)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks,fontsize=14)

    if i == 0:
        ax.set_title('Control Climates',fontsize=16)
    #ax.set_ylabel('SWCF ($W/m^2$)',fontsize=16)
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
    #ax2.set_ylabel('$\Delta$ SWCF ($W/m^2$)',fontsize=16)
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
    #ax3.set_ylabel('diff ($W/m^2$)',fontsize=16)
    if i == 0:
        ax3.set_title('Difference W.R.T. Present',fontsize=16)
    ax3.set_xticks(ticks)
    ax3.set_xticklabels(ticklabels,fontsize=14)
    #ax3.legend(loc='upper left')
    ax3.grid(True)

    if i == 0:
        ax4 = fig.add_subplot(5, 1,4)
        ax4.axhline(0,-1,1,color='k',linestyle='dashed',linewidth=0.5)
    ax4.plot(lat, Quad_response, colors[i], label=Casenames[i], linewidth=2)

    ax4.set_xlim(-1,1)
    ax4.set_yticks(ydiffticks)
    ax4.set_yticklabels(ydiffticks,fontsize=14)
    #ax4.set_ylabel('$\Delta$ SWCF ($W/m^2$)',fontsize=16)
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
    #ax5.set_ylabel('diff ($W/m^2$)',fontsize=16)
    if i == 0:
        ax5.set_title('Difference W.R.T. Present',fontsize=16)
    if i == 3:
        ax5.set_xlabel('Latitude (Deg N)',fontsize=16)
    ax5.set_xticks(ticks)
    ax5.set_xticklabels(ticklabels,fontsize=14)
    #ax3.legend(loc='upper left')
    ax5.grid(True)

    i +=1

if normalize is False:
    fig.suptitle('Shortwave Cloud Forcing',fontsize=24)
else:
    fig.suptitle('Shortwave Cloud Forcing',fontsize=16,y=1.05)
    
'''
if normalize is False:
    fig.suptitle('Zonal Mean Shortwave Cloud Forcing Response to $CO_2$ Doubling',fontsize=24)
else:
    fig.suptitle('Zonal Mean Shortwave Cloud Forcing Response to Equivalent Radiative Forcing',fontsize=12,y=1.05)
'''

handles, labels = ax.get_legend_handles_labels()
#fig.legend(handles,labels,loc = (0.25, 0.93), ncol=4 )
fig.tight_layout(pad=0.25)
plt.show()
if normalize is False:
    fig.savefig(figure_path+'ZM_SWCF_2x_Response.pdf',bbox_inches='tight')
else:
    fig.savefig(figure_path+'ZM_SWCF_n_Response.pdf',bbox_inches='tight')
###################################################################################################################################################


