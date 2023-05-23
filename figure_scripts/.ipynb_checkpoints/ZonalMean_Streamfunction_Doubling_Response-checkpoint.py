#!/usr/bin/env python3
###################################################################################################################################################

#
# Script for plotting figure 3.5. for altercation, make necessary comments and uncomments
# of particular axes settings in the plotting section.
#

import os
import sys
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

figure_path = '/home/brandonsmith/climate_state_thesis/figures/'
casenames = ['0.9','0.95','1.0','1.05']

infield = 'V'
field = 'V'
normalize = False
# For creating smaller file for plotting. Opening parent model output data is too bulky.
for CASENAME in casenames:
    run = 'Control'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_control = '/zmv_'+run+'_'+CASENAME+'.nc'
    outfile_control_stf = '/ZM_Streamfunction_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_control):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo -timmean -zonmean -seltimestep,-360/-1 -select,name='+infield+' '+inpath+infile +' '+ outpath+outfile_control
            print(syscall)
            syscall = '/usr/bin/cdo mastrfu '+outpath+outfile_control+' '+outpath+outfile_control_stf
            print(syscall)
    run = '2xCO2'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'_.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_2xco2 = '/zmv_'+run+'_'+CASENAME+'.nc'
    outfile_2xco2_stf = '/ZM_Streamfunction_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_2xco2):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo -timmean -zonmean -seltimestep,-360/-1 -select,name='+infield+' '+inpath+infile+' '+outpath+outfile_2xco2
            print(syscall)
            syscall = '/usr/bin/cdo mastrfu '+outpath+outfile_2xco2+' '+outpath+outfile_2xco2_stf
            print(syscall)
    run = '4xCO2'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'_.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_4xco2 = '/zmv_'+run+'_'+CASENAME+'.nc'
    outfile_4xco2_stf = '/ZM_Streamfunction_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_4xco2):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo -timmean -zonmean -seltimestep,-360/-1 -select,name='+infield+' '+inpath+infile+' '+outpath+outfile_4xco2
            print(syscall)
            syscall = '/usr/bin/cdo mastrfu '+outpath+outfile_4xco2+' '+outpath+outfile_4xco2_stf
            print(syscall)
            
# ONCE LINES ARE PRINTED, COPY OVER TO A SHELL SCRIPT ON A SERVER WHERE CLIMATE DATA OPERATORS IS INSTALLED AND RUN TO CREATE THE FILES.
# ALTERNATIVELY, REPLACE print() STATEMENTS WITH os.system() FUNCTION CALL TO CALL DIRECTLY TO COMMAND LINE.

# Load variables and perform calculations from outfiles
i = 0
fig = plt.figure(figsize=(20,10))
outfile_1C = '/home/brandonsmith/modeloutput/Control/1.0/ZM_Streamfunction_Control_1.0.nc'
outfile_1D = '/home/brandonsmith/modeloutput/2xCO2/1.0/ZM_Streamfunction_2xCO2_1.0.nc'
outfile_1Q = '/home/brandonsmith/modeloutput/4xCO2/1.0/ZM_Streamfunction_4xCO2_1.0.nc'
rf = [5.53394852, 4.72283731, 3.88779531, 2.78457941] #defined radiative forcing from doubling CO2 (Byrne and Goldblatt (2014))
rfq = [11.4426504, 9.82042797, 8.15034396, 5.94391217] #defined radiative forcing from quadrupling CO2 (Byrne and Goldblatt (2014))
for CASENAME in casenames:
    outfile_2xco2 = '/ZM_Streamfunction_2xCO2_'+CASENAME+'.nc'
    outfile_control = '/ZM_Streamfunction_Control_'+CASENAME+'.nc'
    outfile_4xco2 = '/ZM_Streamfunction_4xCO2_'+CASENAME+'.nc'
    outpath_control = '/home/brandonsmith/modeloutput/Control/'+CASENAME
    outpath_2xco2 = '/home/brandonsmith/modeloutput/2xCO2/'+CASENAME
    outpath_4xco2 = '/home/brandonsmith/modeloutput/4xCO2/'+CASENAME
    dsloc_control = outpath_control+outfile_control
    dsloc_2xco2 = outpath_2xco2+outfile_2xco2
    dsloc_4xco2 = outpath_4xco2+outfile_4xco2

    if os.path.isfile(dsloc_control) and os.path.isfile(dsloc_2xco2):
        # Load Variables
        dsc = netCDF4.Dataset(dsloc_control)
        Cvar= dsc.variables['mastrfu'][:]
        lat = dsc.variables['lat'][:]
        lev = dsc.variables['lev'][:]
        dsc.close()

        dsd = netCDF4.Dataset(dsloc_2xco2)
        Dvar = dsd.variables['mastrfu'][:]
        # no need to reload lat/lon
        dsd.close()

        dsq = netCDF4.Dataset(dsloc_4xco2)
        Qvar = dsq.variables['mastrfu'][:]
        # no need to reload lat/lon
        dsq.close()

        C1 = netCDF4.Dataset(outfile_1C)
        Cvar1 = C1.variables['mastrfu'][:]
        C1.close()
        D1 = netCDF4.Dataset(outfile_1D)
        Dvar1 = D1.variables['mastrfu'][:]
        D1.close()
        Q1 = netCDF4.Dataset(outfile_1Q)
        Qvar1 = Q1.variables['mastrfu'][:]
        Q1.close()
        Cvar1 = Cvar1.squeeze()
        Dvar1 = Dvar1.squeeze()
        Qvar1 = Qvar1.squeeze()
        Cvar = Cvar.squeeze()
        Dvar = Dvar.squeeze()
        Qvar = Qvar.squeeze()
        
        # Normalize responses with respect to their radiative forcings. If you wish to see raw doubling responses, set "normalize" to False.
        if normalize is True:
            doubling_response = rf[2]*(Dvar - Cvar)/rf[i]
            diff = doubling_response - (Dvar1 - Cvar1)
            Quad_response = rfq[2]*(Qvar-Cvar)/rfq[i]
            Qdiff = Quad_response - (Qvar1-Cvar1)
        else:
            doubling_response = Dvar - Cvar
            diff = doubling_response - (Dvar1 - Cvar1)
            Quad_response = (Qvar-Cvar)
            Qdiff = Quad_response - (Qvar1-Cvar1)

        Casenames = ['90% $S_0$','95% $S_0$','$S_0$','105% $S_0$']
        #plot variable
        ax = fig.add_subplot(4,5,5*i+1)
        lat2d, lev2d = np.meshgrid(lat, lev)
        cs = ax.contourf(lat,lev,Cvar,11, cmap='RdBu_r',vmin=-1e9,vmax=1e9)
        #x = plt.contour(lat,lev,np.squeeze(Cvar),colors='black',linewidths=0.75)
        #plt.clabel(x,fontsize='x-small',fmt='%1.0f')
        #plt.clabel(y,fontsize='x-small',fmt='%1.0f')
        #cs.set_edgecolor("face")
        ax.set_yticks(np.linspace(1000,0,5),fontsize=12)
        plt.gca().invert_yaxis()
        ax.set_xticks(np.linspace(-90,90,7),fontsize=12)
        #ax.set_ylabel(Casenames[i],fontsize=14)
        if i == 0:
            ax.set_title('Base Climate Comparison',fontsize=16)
            #cticks=np.linspace(0,50,6)
            cticks=np.around(np.linspace(-1e9,1e9,5),decimals=1)
            cbar_ax = fig.add_axes([0.02,-0.05,0.18,0.02])
            #cb = fig.colorbar(cs, spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
            cb = fig.colorbar(mappable=None, norm=Normalize(vmin=-1e9,vmax=1e9), cmap='RdBu_r',spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
            cb.set_label("Stream Function ($m^2/s$)",fontsize=14)
            #cb.formatter.set_powerlimits((0,0))
            #cb.set_ticklabels(cticks.astype(str),fontsize=12)
            #cticklabels=np.linspace(-10,10,11)
            #cb.set_ticklabels(cticks.astype(int).astype(str),fontsize=12)
        ax.text(0.05, 0.9, Casenames[i], va='center', rotation='horizontal',transform=ax.transAxes,fontsize=16,color='black')
        if i == 3:
            ax.set_xlabel('Latitude (Deg N)',fontsize=12)
        ax2 = fig.add_subplot(4,5,5*i+2)
        lat2d, lev2d = np.meshgrid(lat, lev)
        cs2 = ax2.contourf(lat,lev,doubling_response,101, cmap='RdBu_r',vmin=-2e8,vmax=2e8)
        #cs2 = plt.contourf(lat,lev,doubling_response,60, cmap='BrBG',vmin=-10,vmax=10)

        x = ax2.contour(lat,lev,Dvar,11,colors='black',linewidths=1)
        #plt.clabel(x,fontsize='x-small',fmt='%1.0f')
        #plt.clabel(y,fontsize='x-small',fmt='%1.0f')
        #cs2.set_edgecolor("face")
        ax2.set_yticks(np.linspace(1000,0,5),fontsize=12)
        plt.gca().invert_yaxis()
        ax2.set_xticks(np.linspace(-90,90,7),fontsize=12)
        if i == 0:
            if normalize is False:
                ax2.set_title('Response to $CO_2$ Doubling',fontsize=16)
            else:
                ax2.set_title('Response to $F_2$$_x$',fontsize=16)
            #cticks=np.linspace(-50,50,11)
            cticks=np.around(np.linspace(-2e8,2e8,5),decimals=2)
            cbar_ax = fig.add_axes([0.22,-0.05,0.18,0.02])
            cb2 = fig.colorbar(mappable=None, norm=Normalize(vmin=-2e8,vmax=2e8), cmap='RdBu_r',spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
            #cb2 = fig.colorbar(cs2, spacing='proportional',orientation='horizontal',cax=cbar_ax)
            cb2.set_label("Change in Stream Function",fontsize=14)
            #cb2.formatter.set_powerlimits((0,0))
            #cticklabels=np.linspace(-2,2,5)
            #cb2.set_ticklabels(cticks.astype(int).astype(str),fontsize=12)
            #cb2.set_ticklabels(cticks.astype(str),fontsize=12)
        if i == 3:
            ax2.set_xlabel('Latitude (Deg N)',fontsize=12)
        ax3 = fig.add_subplot(4,5,5*i+3)
        lat2d, lev2d = np.meshgrid(lat, lev)
        cs3 = ax3.contourf(lat,lev,diff,101, cmap='RdBu_r',vmin=-2e8,vmax=2e8)
        #cs3 = plt.contourf(lat,lev,diff,60, cmap='BrBG',vmin=-10,vmax=10)

        x = ax3.contour(lat,lev,Dvar,11,colors='black',linewidths=1)
        #plt.clabel(x,fontsize='x-small',fmt='%1.0f')
        #plt.clabel(y,fontsize='x-small',fmt='%1.0f')
        #cs3.set_edgecolor("face")
        ax3.set_yticks(np.linspace(1000,0,5),fontsize=12)
        plt.gca().invert_yaxis()
        ax3.set_xticks(np.linspace(-90,90,7),fontsize=12)
        if i == 0:
            ax3.set_title('Difference',fontsize=16)
            #cticks=np.linspace(-50,50,11)
            cticks=np.around(np.linspace(-2e8,2e8,5),decimals=2)
            cbar_ax = fig.add_axes([0.42,-0.05,0.18,0.02])
            cb3 = fig.colorbar(mappable=None, norm=Normalize(vmin=-2e8,vmax=2e8), cmap='RdBu_r',spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
            #cb3 = fig.colorbar(cs3, spacing='proportional',orientation='horizontal',cax=cbar_ax)
            cb3.set_label('Difference w.r.t. $S_0$',fontsize=14)
            #cb3.formatter.set_powerlimits((0,0))
            #cticklabels=np.linspace(-2,2,5)
            #cb3.set_ticklabels(cticks.astype(int).astype(str),fontsize=12)
            #cb3.set_ticklabels(cticks.astype(str),fontsize=12)
        if i == 3:
            ax3.set_xlabel('Latitude (Deg N)',fontsize=12)
        ax4 = fig.add_subplot(4,5,5*i+4)
        lat2d, lev2d = np.meshgrid(lat, lev)
        cs4 = ax4.contourf(lat,lev,Quad_response,101, cmap='RdBu_r',vmin=-2e8,vmax=2e8)
        #cs4 = plt.contourf(lat,lev,Quad_response,60, cmap='BrBG',vmin=-10,vmax=10)
        x = ax4.contour(lat,lev,Qvar,11,colors='black',linewidths=1)
        #plt.clabel(x,fontsize='x-small',fmt='%1.0f')
        #plt.clabel(y,fontsize='x-small',fmt='%1.0f')
        #cs3.set_edgecolor("face")
        ax4.set_yticks(np.linspace(1000,0,5),fontsize=12)
        plt.gca().invert_yaxis()
        ax4.set_xticks(np.linspace(-90,90,7),fontsize=12)
        if i == 0:
            if normalize is False:
                ax4.set_title('Response to $CO_2$ Quadrupling',fontsize=16)
            else:
                ax4.set_title('Response to $F_4$$_x$',fontsize=16)
            #cticks=np.linspace(-50,50,11)
            cticks=np.around(np.linspace(-2e8,2e8,5),decimals=2)
            cbar_ax = fig.add_axes([0.62,-0.05,0.18,0.02])
            cb4 = fig.colorbar(mappable=None, norm=Normalize(vmin=-2e8,vmax=2e8), cmap='RdBu_r',spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
            #cb4 = fig.colorbar(cs4, spacing='proportional',orientation='horizontal',cax=cbar_ax)
            cb4.set_label('Change in Stream Function',fontsize=14)
            #cb4.formatter.set_powerlimits((0,0))
            #cticklabels=np.linspace(-2,2,5)
            #cb4.set_ticklabels(cticks.astype(int).astype(str),fontsize=12)
            #cb4.set_ticklabels(cticks.astype(str),fontsize=12)
        if i == 3:
            ax4.set_xlabel('Latitude (Deg N)',fontsize=12)
        ax5 = fig.add_subplot(4,5,5*i+5)
        lat2d, lev2d = np.meshgrid(lat, lev)
        cs5 = ax5.contourf(lat,lev,Qdiff,101, cmap='RdBu_r',vmin=-2e8,vmax=2e8)
        #cs5 = plt.contourf(lat,lev,Qdiff,60, cmap='BrBG',vmin=-10,vmax=10)
        x = ax5.contour(lat,lev,Qvar,11,colors='black',linewidths=1)
        #plt.clabel(x,fontsize='x-small',fmt='%1.0f')
        #plt.clabel(y,fontsize='x-small',fmt='%1.0f')
        #cs3.set_edgecolor("face")
        ax5.set_yticks(np.linspace(1000,0,5),fontsize=12)
        plt.gca().invert_yaxis()
        ax5.set_xticks(np.linspace(-90,90,7),fontsize=12)
        if i == 0:
            ax5.set_title('Difference w.r.t. $S_0$',fontsize=16)
            #cticks=np.linspace(-50,50,11)
            cticks=np.around(np.linspace(-2e8,2e8,5),decimals=2)
            cbar_ax = fig.add_axes([0.82,-0.05,0.18,0.02])
            cb5 = fig.colorbar(mappable=None, norm=Normalize(vmin=-2e8,vmax=2e8), cmap='RdBu_r',spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
            #cb5 = fig.colorbar(cs5, spacing='proportional',orientation='horizontal',cax=cbar_ax)
            cb5.set_label('Difference',fontsize=14)
            #cb5.formatter.set_powerlimits((0,0))
            #cticklabels=np.linspace(-2,2,5)
            #cb5.set_ticklabels(cticks.astype(int).astype(str),fontsize=12)
            #cb5.set_ticklabels(cticks.astype(str),fontsize=12)
        if i == 3:
            ax5.set_xlabel('Latitude (Deg N)',fontsize=12)
#            #cb1 = plt.colorbar(cs1,spacing='proportional',shrink=0.8)
        #cb1.set_label('(K)')

        i = i+1
    else:
        print('No such file or directory:'+dsloc_control)
if normalize is False:        
    fig.suptitle('Zonal Mean Stream Function Response to $CO_2$ Doubling',y=1.01,fontsize=24)
else:
    fig.suptitle('Zonal Mean Stream Function Response to Equivalent Radiative Forcing',y=1.01,fontsize=24)

fig.tight_layout(pad=0.2)
fig.text(-0.02, 0.5, 'Hybrid Sigma-Pressure Level (mb)', va='center', rotation='vertical',fontsize=16)
plt.show()
if normalize is False:
    fig.savefig(figure_path+'Zonal_Mean_streamfunction_2x_response.pdf',bbox_inches='tight')
else:
    fig.savefig(figure_path+'Zonal_Mean_streamfunction_n_response.pdf',bbox_inches='tight')
###################################################################################################################################################
