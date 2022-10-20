#!/usr/bin/env python3

#
# Script for plotting figure 3.7. for altercation, make necessary comments and uncomments
# of particular axes settings in the plotting section.
#

import os
import time
import sys
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from matplotlib.colors import Normalize
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point

figure_path = '/home/brandonsmith/climate_state_thesis/figures/'

TS = 'TS'
ICEFRAC = 'ICEFRAC'
QFLX = 'QFLX'
CLDLOW = 'CLDLOW'
casenames = ['0.9','0.95','1.0','1.05']
normalize = True

for CASENAME in casenames:
    run = 'Control'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_control = '/map_polar_clouds_and_ice_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_control):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo timmean -seltimestep,-360/-1 -select,name='+TS+','+ICEFRAC+','+QFLX+','+CLDLOW+' '+inpath+infile +' '+ outpath+outfile_control
            print(syscall)
    run = '2xCO2'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'_.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_2xco2 = '/map_polar_clouds_and_ice_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_2xco2):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo timmean -seltimestep,-360/-1 -select,name='+TS+','+ICEFRAC+','+QFLX+','+CLDLOW+' '+inpath+infile +' '+ outpath+outfile_2xco2
            print(syscall)
    run = '4xCO2'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'_.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_4xco2 = '/map_polar_clouds_and_ice_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_4xco2):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo timmean -seltimestep,-360/-1 -select,name='+TS+','+ICEFRAC+','+QFLX+','+CLDLOW+' '+inpath+infile +' '+ outpath+outfile_4xco2
            print(syscall)
# ONCE LINES ARE PRINTED, COPY OVER TO A SHELL SCRIPT ON A SERVER WHERE CLIMATE DATA OPERATORS IS INSTALLED AND RUN TO CREATE THE FILES.
# ALTERNATIVELY, REPLACE print() STATEMENTS WITH os.system() FUNCTION CALL TO CALL DIRECTLY TO COMMAND LINE.    

    # Load variables and perform calculations from outfiles
i = 0
fig = plt.figure(figsize=(10,8))
outfile_1C = '/home/brandonsmith/modeloutput/Control/1.0/map_polar_clouds_and_ice_Control_1.0.nc'
outfile_1D = '/home/brandonsmith/modeloutput/2xCO2/1.0/map_polar_clouds_and_ice_2xCO2_1.0.nc'
outfile_1Q = '/home/brandonsmith/modeloutput/4xCO2/1.0/map_polar_clouds_and_ice_4xCO2_1.0.nc'
rf = [5.53394852, 4.72283731, 3.88779531, 2.78457941] #defined radiative forcing from doubling CO2 (Byrne and Goldblatt (2014))
rfq = [11.4426504, 9.82042797, 8.15034396, 5.94391217] #defined radiative forcing from quadrupling CO2 (Byrne and Goldblatt (2014))
for CASENAME in casenames:
    outfile_2xco2 = '/map_polar_clouds_and_ice_2xCO2_'+CASENAME+'.nc'
    outfile_control = '/map_polar_clouds_and_ice_Control_'+CASENAME+'.nc'
    outfile_4xco2 = '/map_polar_clouds_and_ice_4xCO2_'+CASENAME+'.nc'
    outpath_control = '/home/brandonsmith/modeloutput/Control/'+CASENAME
    outpath_2xco2 = '/home/brandonsmith/modeloutput/2xCO2/'+CASENAME
    outpath_4xco2 = '/home/brandonsmith/modeloutput/4xCO2/'+CASENAME
    dsloc_control = outpath_control+outfile_control
    dsloc_2xco2 = outpath_2xco2+outfile_2xco2
    dsloc_4xco2 = outpath_4xco2+outfile_4xco2
    if os.path.isfile(dsloc_control) and os.path.isfile(dsloc_2xco2):
        dsc = netCDF4.Dataset(dsloc_control)
        icefrac = dsc.variables[ICEFRAC][:]
        lowcloud = dsc.variables[CLDLOW][:]
        lat = dsc.variables['lat'][:]
        lon = dsc.variables['lon'][:]
        #lev = dsc.variables['lev'][:]
        dsc.close()
        
        dsd = netCDF4.Dataset(dsloc_2xco2)
        icefracd = dsd.variables[ICEFRAC][:]
        lowcloudd = dsd.variables[CLDLOW][:]
        #lev = dsc.variables['lev'][:]
        dsd.close()
        
        dsq = netCDF4.Dataset(dsloc_4xco2)
        icefracq = dsq.variables[ICEFRAC][:]
        lowcloudq = dsq.variables[CLDLOW][:]
        dsq.close()
        
        C1 = netCDF4.Dataset(outfile_1C)
        icefrac1 = C1.variables[ICEFRAC][:]
        lowcloud1 = C1.variables[CLDLOW][:]
        C1.close()
        D1 = netCDF4.Dataset(outfile_1D)
        icefrac1d = D1.variables[ICEFRAC][:]
        lowcloud1d = D1.variables[CLDLOW][:]
        D1.close()
        
        Q1 = netCDF4.Dataset(outfile_1Q)
        icefrac1q = Q1.variables[ICEFRAC][:]
        lowcloud1q = Q1.variables[CLDLOW][:]
        Q1.close()
        
        Lon = lon
        icefrac = np.squeeze(icefrac)
        icefrac, Lon = add_cyclic_point(icefrac, coord=lon)
        Lon = lon
        lowcloud = np.squeeze(lowcloud)
        lowcloud, Lon = add_cyclic_point(lowcloud, coord=lon)
        Lon = lon
        icefracd = np.squeeze(icefracd)
        icefracd, Lon = add_cyclic_point(icefracd, coord=lon)
        Lon = lon
        icefracq = np.squeeze(icefracq)
        icefracq, Lon = add_cyclic_point(icefracq, coord=lon)
        Lon = lon
        lowcloudd = np.squeeze(lowcloudd)
        lowcloudd, Lon = add_cyclic_point(lowcloudd, coord=lon)
        Lon = lon
        lowcloudq = np.squeeze(lowcloudq)
        lowcloudq, Lon = add_cyclic_point(lowcloudq, coord=lon)
        Lon = lon
        icefrac1 = np.squeeze(icefrac1)
        icefrac1, Lon = add_cyclic_point(icefrac1, coord=lon)
        Lon = lon
        lowcloud1 = np.squeeze(lowcloud1)
        lowcloud1, Lon = add_cyclic_point(lowcloud1, coord=lon)
        Lon = lon
        icefrac1d = np.squeeze(icefrac1d)
        icefrac1d, Lon = add_cyclic_point(icefrac1d, coord=lon)
        Lon = lon
        icefrac1q = np.squeeze(icefrac1q)
        icefrac1q, Lon = add_cyclic_point(icefrac1q, coord=lon)
        Lon = lon
        lowcloud1d = np.squeeze(lowcloud1d)
        lowcloud1d, Lon = add_cyclic_point(lowcloud1d, coord=lon)
        Lon = lon
        lowcloud1q = np.squeeze(lowcloud1q)
        lowcloud1q, Lon = add_cyclic_point(lowcloud1q, coord=lon)
        
        
        if normalize is True:
            icefracDR = rf[2]*(icefracd - icefrac)/rf[i]
            lowcloudDR = rf[2]*(lowcloudd - lowcloud)/rf[i]
            icefracQR = rfq[2]*(icefracq - icefrac)/rfq[i]
            lowcloudQR = rfq[2]*(lowcloudq - lowcloud)/rfq[i]
        else:
            icefracDR = icefracd - icefrac
            lowcloudDR = lowcloudd - lowcloud
            icefracQR = icefracq - icefrac
            lowcloudQR = lowcloudq - lowcloud
        
        #plot variable
        ax = fig.add_subplot(4,5,5*i+3,projection=ccrs.SouthPolarStereo())
        ax.set_extent([-180,180,-90,-60],ccrs.PlateCarree())
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax.set_boundary(circle, transform=ax.transAxes)
        
        cs = ax.contourf(Lon,lat,lowcloudDR,11,transform=ccrs.PlateCarree(),cmap='BrBG',vmin=-0.2,vmax=0.2)
        ax.coastlines()
        ax.set_ylabel('Solar Multiplier: '+CASENAME)
        ax.gridlines(draw_labels=True,color='gray',linewidth=0.25,linestyle='dotted')
        
        if i == 0:
            if normalize is True:
                ax.set_title('$F_2$$_x$ Response',fontsize=11)
            else:
                ax.set_title('$2xCO_2$ Response',fontsize=11)
            cticks=np.around(np.linspace(-0.2,0.2,5),decimals=1)
            cbar_ax = fig.add_axes([0.41,-0.05,0.18,0.02])
            #cbar_ax = fig.add_axes([0,-0.05,0.23,0.02])
            #cb = fig.colorbar(mappable=None, norm=Normalize(vmin=-0.2,vmax=0.2), cmap='BrBG',spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
            cb = fig.colorbar(cs, spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
            #cb.formatter.set_powerlimits((0,0))
            cb.set_ticklabels(cticks.astype(str),fontsize=10)
            cb.set_label('Cloud Fraction Change',fontsize=12)
                
        ax2 = fig.add_subplot(4,5,5*i+2,projection=ccrs.SouthPolarStereo())
        ax2.set_extent([-180,180,-90,-60],ccrs.PlateCarree())
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax2.set_boundary(circle, transform=ax2.transAxes)
        
        cs2 = ax2.contourf(Lon,lat,icefracDR,11,transform=ccrs.PlateCarree(),cmap='RdBu_r',vmin=-1,vmax=1)
        ax2.coastlines()
        ax2.gridlines(draw_labels=True,color='gray',linewidth=0.25,linestyle='dotted')
        
        if i == 0:
            if normalize is True:
                ax2.set_title('$F_2$$_x$ Response',fontsize=11)
            else:
                ax2.set_title('$2xCO_2$ Response',fontsize=11)
            cticks=np.around(np.linspace(-1,1,11),decimals=1)
            cbar_ax = fig.add_axes([0.21,-0.05,0.18,0.02])
            cb2 = fig.colorbar(cs2, spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
            #cb2 = fig.colorbar(mappable=None, norm=Normalize(vmin=-1,vmax=1), cmap='RdBu_r',spacing='uniform',orientation='horizontal',cax=cbar_ax,ticks=cticks)
            cb2.set_label('Ice Fraction Change',fontsize=12)
            cb2.set_ticklabels(cticks.astype(str),fontsize=10)
       # cb2.formatter.set_powerlimits((0,0))
    
        ax3 = fig.add_subplot(4,5,5*i+1,projection=ccrs.SouthPolarStereo())
        ax3.set_extent([-180,180,-90,-60],ccrs.PlateCarree())
        
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax3.set_boundary(circle, transform=ax3.transAxes)
        cs3 = ax3.contourf(Lon,lat,icefrac,10,transform=ccrs.PlateCarree(),cmap='Blues_r',vmin=0,vmax=1)
        ax3.coastlines()
        ax3.gridlines(draw_labels=True,color='gray',linewidth=0.25,linestyle='dotted')
        
        if i == 0:
            ax3.set_title('Initial Ice Fraction',fontsize=12)
            cticks=np.around(np.linspace(0,1,6),decimals=1)
            cbar_ax = fig.add_axes([0,-0.05,0.18,0.02])
            #cbar_ax = fig.add_axes([0.51,-0.05,0.23,0.02])
            #cb3 = fig.colorbar(mappable=None, norm=Normalize(vmin=0,vmax=1), cmap='Blues_r',spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
            cb3 = plt.colorbar(cs3,spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
            cb3.set_label('Ice Fraction',fontsize=12)
            cb3.set_ticklabels(cticks.astype(str),fontsize=10)
        
        ax4 = fig.add_subplot(4,5,5*i+4,projection=ccrs.SouthPolarStereo())
        ax4.set_extent([-180,180,-90,-60],ccrs.PlateCarree())
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax4.set_boundary(circle, transform=ax4.transAxes)
        cs4 = ax4.contourf(Lon,lat,icefracQR,11,transform=ccrs.PlateCarree(),cmap='RdBu_r',vmin=-1,vmax=1)
        ax4.coastlines()
        ax4.gridlines(draw_labels=True,color='gray',linewidth=0.25,linestyle='dotted')
        
        if i == 0:
            if normalize is True:
                ax4.set_title('$F_4$$_x$ Response',fontsize=11)
            else:
                ax4.set_title('$4xCO_2$ Response',fontsize=11)
        if i == 2:
            cticks=np.around(np.linspace(-1,1,11),decimals=1)
            cbar_ax = fig.add_axes([0.61,-0.05,0.18,0.02])
            #cb4 = fig.colorbar(mappable=None, norm=Normalize(vmin=-1,vmax=1), cmap='RdBu_r',spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
            cb4 = fig.colorbar(cs4,spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
            cb4.set_label('Ice Fraction Change',fontsize=12)
            cb4.set_ticklabels(cticks.astype(str),fontsize=10)
            
        ax5 = fig.add_subplot(4,5,5*i+5,projection=ccrs.SouthPolarStereo())
        ax5.set_extent([-180,180,-90,-60],ccrs.PlateCarree())
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax5.set_boundary(circle, transform=ax5.transAxes)
        
        cs5 = ax5.contourf(Lon,lat,lowcloudQR,11,transform=ccrs.PlateCarree(),cmap='BrBG',vmin=-0.2,vmax=0.2)
        ax5.coastlines()
        ax5.set_ylabel('Solar Multiplier: '+CASENAME)
        ax5.gridlines(draw_labels=True,color='gray',linewidth=0.25,linestyle='dotted')
        
        if i == 0:
            if normalize is True:
                ax5.set_title('$F_4$$_x$ Response',fontsize=11)
            else:
                ax5.set_title('$4xCO_2$ Response',fontsize=11)
        if i == 2:
            cticks=np.around(np.linspace(-0.2,0.2,5),decimals=1)
            cbar_ax = fig.add_axes([0.81,-0.05,0.18,0.02])
            #cbar_ax = fig.add_axes([0,-0.05,0.23,0.02])
            #cb5 = fig.colorbar(mappable=None, norm=Normalize(vmin=-0.2,vmax=0.2), cmap='BrBG',spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
            cb5 = fig.colorbar(cs5, spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
            #cb.formatter.set_powerlimits((0,0))
            cb5.set_ticklabels(cticks.astype(str),fontsize=10)
            cb5.set_label('Cloud Fraction Change',fontsize=12)
            
        i = i+1
    else:
        print('No such file or directory: '+dsloc_control+' or '+dsloc_4xco2)
        
if normalize is True:
    plt.suptitle('Antarctic Low Cloud Fraction and Sea Ice Response to Equivalent Radiative Forcing',y=1.01,fontsize=20)
else:
    plt.suptitle('Antarctic Low Cloud Fraction and Sea Ice Response to $CO_2$ Doublings',y=1.01,fontsize=20)
fig.text(-0.02, 0.85, 'Solar Constant: 90%', va='center', ha='left', rotation='vertical',fontsize=10)
fig.text(-0.02, 0.61, 'Solar Constant: 95%', va='center', ha='left', rotation='vertical',fontsize=10)
fig.text(-0.02, 0.37, 'Solar Constant: 100%', va='center', ha='left', rotation='vertical',fontsize=10)
fig.text(-0.02, 0.13, 'Solar Constant: 105%', va='center', ha='left', rotation='vertical',fontsize=10)
fig.tight_layout(pad=0.2)
fig.subplots_adjust(hspace=0.05,wspace=0.05)
plt.show()
if normalize is True:
    fig.savefig(figure_path+'map_Antarctic_clouds_and_ice_n.pdf',bbox_inches='tight')
else:
    fig.savefig(figure_path+'map_Antarctic_clouds_and_ice.pdf',bbox_inches='tight')