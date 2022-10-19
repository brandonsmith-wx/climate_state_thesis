#!/usr/bin/env python3

#
# Script for plotting figure 3.2. for altercation, make necessary comments and uncomments
# of particular axes settings in the plotting section.
#

import os
import time
import sys
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import cartopy
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point

figure_path = '/home/brandonsmith/climate-gcm-bps/plots/'
casenames = ['0.9','0.95','1.0','1.05']
normalize = True

field = 'TS'
field2 = 'TREFHT'
for CASENAME in casenames:
    run = 'Control'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME #location of parent data set
    infile = '/Merged_'+run+'_'+CASENAME+'.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME #location of child data set
    outfile_control = '/map_TS_'+run+'_'+CASENAME+'.nc'
    if os.path.isfile(inpath+infile): #check to see if the parent file exists
        if not os.path.isfile(outpath+outfile_control): #create plotting file if it does not exist. Else, continue to plotting
            syscall = '/usr/bin/cdo timmean -seltimestep,-360/-1 -select,name='+field+','+field2+' '+inpath+infile+ ' '+outpath+outfile_control
            print(syscall)
    run = '2xCO2'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_control = '/map_TS_'+run+'_'+CASENAME+'.nc'
    if os.path.isfile(inpath+infile):
        if not os.path.isfile(outpath+outfile_control):
            syscall = '/usr/bin/cdo timmean -seltimestep,-360/-1 -select,name='+field+','+field2+' '+inpath+infile+ ' '+outpath+outfile_control
            print(syscall)
    run = '4xCO2'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_control = '/map_TS_'+run+'_'+CASENAME+'.nc'
    if os.path.isfile(inpath+infile):
        if not os.path.isfile(outpath+outfile_control):
            syscall = '/usr/bin/cdo timmean -seltimestep,-360/-1 -select,name='+field+','+field2+' '+inpath+infile+ ' '+outpath+outfile_control
            print(syscall)
# ONCE LINES ARE PRINTED, COPY OVER TO A SHELL SCRIPT ON A SERVER WHERE CLIMATE DATA OPERATORS IS INSTALLED AND RUN TO CREATE THE FILES.
# ALTERNATIVELY, REPLACE print() STATEMENTS WITH os.system() FUNCTION CALL TO CALL DIRECTLY TO COMMAND LINE.

# Load variables and perform calculations from outfiles
i = 0
fig = plt.figure(figsize=(20,10))
outfile_1C = '/home/brandonsmith/modeloutput/Control/1.0/map_TS_Control_1.0.nc'
outfile_1D = '/home/brandonsmith/modeloutput/2xCO2/1.0/map_TS_2xCO2_1.0.nc'
outfile_1Q = '/home/brandonsmith/modeloutput/4xCO2/1.0/map_TS_4xCO2_1.0.nc'
rf = [5.53394852, 4.72283731, 3.88779531, 2.78457941] #defined radiative forcing from doubling CO2 (Byrne and Goldblatt (2014))
rfq = [11.4426504, 9.82042797, 8.15034396, 5.94391217] #defined radiative forcing from quadrupling CO2 (Byrne and Goldblatt (2014))
for CASENAME in casenames:
    outfile_2xco2 = '/map_TS_2xCO2_'+CASENAME+'.nc'
    outfile_control = '/map_TS_Control_'+CASENAME+'.nc'
    outfile_4xco2 = '/map_TS_4xCO2_'+CASENAME+'.nc'
    outpath_control = '/home/brandonsmith/modeloutput/Control/'+CASENAME
    outpath_2xco2 = '/home/brandonsmith/modeloutput/2xCO2/'+CASENAME
    outpath_4xco2 = '/home/brandonsmith/modeloutput/4xCO2/'+CASENAME
    dsloc_control = outpath_control+outfile_control
    dsloc_2xco2 = outpath_2xco2+outfile_2xco2
    dsloc_4xco2 = outpath_4xco2+outfile_4xco2
    if os.path.isfile(dsloc_control) and os.path.isfile(dsloc_2xco2):
        # Load in variables
        dsc = netCDF4.Dataset(dsloc_control)
        Cvar = dsc.variables[field][:]
        lat = dsc.variables['lat'][:]
        lon = dsc.variables['lon'][:]
        dsc.close()

        dsd = netCDF4.Dataset(dsloc_2xco2)
        Dvar = dsd.variables[field][:]
        # no need to reload lat/lon
        dsd.close()

        C1 = netCDF4.Dataset(outfile_1C)
        Cvar1 = C1.variables[field][:]
        C1.close()
        D1 = netCDF4.Dataset(outfile_1D)
        Dvar1 = D1.variables[field][:]
        D1.close()

        dsq = netCDF4.Dataset(dsloc_4xco2)
        Qvar = dsq.variables[field][:]
        dsq.close()

        Q1 = netCDF4.Dataset(outfile_1Q)
        Qvar1 = Q1.variables[field][:]
        Q1.close()


        Cvar = Cvar.squeeze()
        Dvar = Dvar.squeeze()
        Qvar = Qvar.squeeze()
        Cvar1 = Cvar1.squeeze()
        Dvar1 = Dvar1.squeeze()
        Qvar1 = Qvar1.squeeze()

        # these lines get rid of the irritating white space line through the center of each plot. artifact of cartopy.
        Lon = lon
        Cvar1 = np.squeeze(Cvar1)
        Cvar1, Lon = add_cyclic_point(Cvar1, coord=lon)
        Lon = lon
        Dvar1 = np.squeeze(Dvar1)
        Dvar1, Lon = add_cyclic_point(Dvar1, coord=lon)
        Lon = lon
        Qvar1 = np.squeeze(Qvar1)
        Qvar1, Lon = add_cyclic_point(Qvar1, coord=lon)
        Lon = lon
        Cvar = np.squeeze(Cvar)
        Cvar, Lon = add_cyclic_point(Cvar, coord=lon)
        Lon = lon
        Dvar = np.squeeze(Dvar)
        Dvar, Lon = add_cyclic_point(Dvar, coord=lon)
        Lon = lon
        Qvar = np.squeeze(Qvar)
        Qvar, Lon = add_cyclic_point(Qvar, coord=lon)

        print(np.shape(Cvar))
        
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
        #plot variable
        ax = fig.add_subplot(4,5,5*i+1,projection=ccrs.Robinson())
        cs = ax.contourf(Lon,lat,Cvar,11,transform=ccrs.PlateCarree(),cmap='gist_heat',vmin=210,vmax=310)
        ax.coastlines()
        #ax.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k')
        #ax.set_ylabel('Solar Multiplier: '+CASENAME,fontsize=14)
        if i == 0:
            plt.title('Base Climate Comparison',fontsize=18)
        if i == 0:
            cticks=np.around(np.linspace(210,310,11),decimals=1)
            cbar_ax = fig.add_axes([0.02,-0.05,0.15,0.02])
            cb = fig.colorbar(mappable=None, norm=Normalize(vmin=210,vmax=310), cmap='gist_heat',spacing='proportional',orientation='horizontal',cax=cbar_ax)#,ticks=cticks)
            #cb.formatter.set_powerlimits((0,0))
            #cb.set_ticklabels(cticks.astype(str),fontsize=12)
            cb.set_ticklabels(cticks.astype(int).astype(str),fontsize=12)
            cb.set_label('Surface Temperature (K)',fontsize=14)
        ax.text(0.03, 0.5, Casenames[i], va='center', rotation='horizontal',transform=ax.transAxes,fontsize=14,color='black')
        ax2 = fig.add_subplot(4,5,5*i+2,projection=ccrs.Robinson())
        cs2 = ax2.contourf(Lon,lat,doubling_response,201,transform=ccrs.PlateCarree(),cmap='RdBu_r',vmin=-20,vmax=20)
        ax2.coastlines()
        #ax2.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k')
        if i == 0:
            #plt.title('Response to $CO_2$ Doubling',fontsize=18)
            plt.title('Response to $F_2$$_x$',fontsize=18)
        if i == 1:
            cticks=np.around(np.linspace(-20,20,5),decimals=2)
            cbar_ax = fig.add_axes([0.22,-0.05,0.15,0.02])
            #cb2 = fig.colorbar(cs2, spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
            cb2 = fig.colorbar(mappable=None, norm=Normalize(vmin=-20,vmax=20), cmap='RdBu_r',spacing='proportional',orientation='horizontal',cax=cbar_ax)#,ticks=cticks)
            cb2.set_label('Change in Surface Temperature (K)',fontsize=14)
            #cb2.set_ticklabels(cticks.astype(str),fontsize=12)
            cb2.set_ticklabels(cticks.astype(int).astype(str),fontsize=12)
           # cb2.formatter.set_powerlimits((0,0))
        ax3 = fig.add_subplot(4,5,5*i+3,projection=ccrs.Robinson())
        cs3 = ax3.contourf(Lon,lat,diff,201,transform=ccrs.PlateCarree(),cmap='RdBu_r',vmin=-20,vmax=20)
        ax3.coastlines()
        #ax3.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k')
        if i == 0:
            plt.title('Difference w.r.t. $S_0$',fontsize=18)
        if i == 0:
            cticks=np.around(np.linspace(-20,20,5),decimals=2)
            cbar_ax = fig.add_axes([0.42,-0.05,0.15,0.02])
            cb3 = fig.colorbar(mappable=None, norm=Normalize(vmin=-20,vmax=20), cmap='RdBu_r',spacing='proportional',orientation='horizontal',cax=cbar_ax)#,ticks=cticks)
            #cb3 = fig.colorbar(cs3,spacing='uniform',orientation='horizontal',cax=cbar_ax,ticks=cticks)
            cb3.set_label('Difference (K)',fontsize=14)
            #cb3.formatter.set_powerlimits((0,0))
            #cb3.set_ticklabels(cticks.astype(str),fontsize=12)
            cb3.set_ticklabels(cticks.astype(int).astype(str),fontsize=12)
        ax4 = fig.add_subplot(4,5,5*i+4,projection=ccrs.Robinson())
        cs4 = ax4.contourf(Lon,lat,Quad_response,201,transform=ccrs.PlateCarree(),cmap='RdBu_r',vmin=-20,vmax=20)
        ax4.coastlines()
        #ax4.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k')
        if i == 0:
            #plt.title('Response to $CO_2$ Quadrupling',fontsize=18)
            plt.title('Response to $F_4$$_x$',fontsize=18)
        if i == 0:
            cticks=np.around(np.linspace(-20,20,5),decimals=2)
            cbar_ax = fig.add_axes([0.62,-0.05,0.15,0.02])
            cb4 = fig.colorbar(mappable=None, norm=Normalize(vmin=-20,vmax=20), cmap='RdBu_r',spacing='proportional',orientation='horizontal',cax=cbar_ax)#,ticks=cticks)
            #cb3 = fig.colorbar(cs3,spacing='uniform',orientation='horizontal',cax=cbar_ax,ticks=cticks)
            cb4.set_label('Change in Surface Temperature (K)',fontsize=14)
            #cb3.formatter.set_powerlimits((0,0))
            #cb4.set_ticklabels(cticks.astype(str),fontsize=12)
            cb4.set_ticklabels(cticks.astype(int).astype(str),fontsize=12)
        ax5 = fig.add_subplot(4,5,5*i+5,projection=ccrs.Robinson())
        cs5 = ax5.contourf(Lon,lat,Qdiff,201,transform=ccrs.PlateCarree(),cmap='RdBu_r',vmin=-20,vmax=20)
        ax5.coastlines()
        #ax5.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k')
        if i == 0:
            plt.title('Difference w.r.t. $S_0$',fontsize=18)
        if i == 0:
            cticks=np.around(np.linspace(-20,20,5),decimals=2)
            cbar_ax = fig.add_axes([0.82,-0.05,0.15,0.02])
            cb5 = fig.colorbar(mappable=None, norm=Normalize(vmin=-20,vmax=20), cmap='RdBu_r',spacing='proportional',orientation='horizontal',cax=cbar_ax)#,ticks=cticks)
            #cb3 = fig.colorbar(cs3,spacing='uniform',orientation='horizontal',cax=cbar_ax,ticks=cticks)
            cb5.set_label('Difference (K)',fontsize=14)
            #cb3.formatter.set_powerlimits((0,0))
            #cb5.set_ticklabels(cticks.astype(str),fontsize=12)
            cb5.set_ticklabels(cticks.astype(int).astype(str),fontsize=12)
        #if i == 2:
         #   ax = fig.add_subplot(4,3,7,projection=ccrs.PlateCarree())
          #  cs = ax.contourf(Lon,lat,Cvar1,60,transform=ccrs.PlateCarree(),cmap='gist_heat',vmin=0,vmax=0.2)
           # ax.coastlines()
            #ax.set_ylabel('Solar Multiplier: '+CASENAME)
            #cticks=np.around(np.linspace(0,0.2,6),decimals=2)
            #cbar_ax = fig.add_axes([0.01,0.225,0.01,0.2])
            #cb4 = fig.colorbar(cs, spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=cticks,ticklocation='left')
            #cb4.set_ticklabels(cticks.astype(str),fontsize=6)
                #cb3.formatter.set_powerlimits((0,0))
        i = i+1
    else:
        print('No such file or directory: '+dsloc_control+' or '+dsloc_4xco2)

#    cb.set_label('Temperature (K)')
#plt.suptitle('Surface Temperature Response to $CO_2$ Doubling',fontsize=24)
plt.suptitle('Surface Temperature Response to Equivalent Radiative Forcing',fontsize=24)
fig.tight_layout(pad=0.2)
plt.show()
fig.savefig(figure_path+'map_TS_n_response.pdf',bbox_inches='tight')