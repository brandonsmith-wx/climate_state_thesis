#!/usr/bin/env python3
import os
import time
import sys
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
import matplotlib.path as mpath

TS = 'TGCLDCWP'
ICEFRAC = 'ICEFRAC'
QFLX = 'FSUTOA'
CLDLOW = 'CLDLOW'
FLDS = 'LWCF'

figure_path = '/home/brandonsmith/climate-gcm-bps/plots/'
casenames = ['0.9','0.95','1.0','1.05']

run = 'Control'
inpath = '/home/brandonsmith/modeloutput/'+run+'/1.0'
infile = '/Merged_'+run+'_1.0_.nc'
outpath = '/home/brandonsmith/modeloutput/'+run+'/1.0'
outfile_control = '/map_polar_clouds_and_ice_seasonal_'+run+'_1.0.nc'
if not os.path.isfile(outpath+outfile_control):
    if os.path.isfile(inpath+infile):
        syscall = '/usr/bin/cdo yseasmean -seltimestep,-360/-1 -select,name='+TS+','+ICEFRAC+','+QFLX+','+CLDLOW+','+FLDS+' '+inpath+infile +' '+ outpath+outfile_control
        os.system(syscall)
run = '2xCO2'
inpath = '/home/brandonsmith/modeloutput/'+run+'/1.0'
infile = '/Merged_'+run+'_1.0.nc'
outpath = '/home/brandonsmith/modeloutput/'+run+'/1.0'
outfile_2xco2 = '/map_polar_clouds_and_ice_seasonal_'+run+'_1.0.nc'
if not os.path.isfile(outpath+outfile_2xco2):
    if os.path.isfile(inpath+infile):
        syscall = '/usr/bin/cdo yseasmean -seltimestep,-360/-1 -select,name='+TS+','+ICEFRAC+','+QFLX+','+CLDLOW+','+FLDS+' '+inpath+infile +' '+ outpath+outfile_2xco2
        os.system(syscall)
run = '4xCO2'
inpath = '/home/brandonsmith/modeloutput/'+run+'/1.0'
infile = '/Merged_'+run+'_1.0.nc'
outpath = '/home/brandonsmith/modeloutput/'+run+'/1.0'
outfile_4xco2 = '/map_polar_clouds_and_ice_seasonal_'+run+'_1.0.nc'
if not os.path.isfile(outpath+outfile_4xco2):
    if os.path.isfile(inpath+infile):
        syscall = '/usr/bin/cdo yseasmean -seltimestep,-360/-1 -select,name='+TS+','+ICEFRAC+','+QFLX+','+CLDLOW+','+FLDS+' '+inpath+infile +' '+ outpath+outfile_4xco2
        os.system(syscall)
    
    # Load variables and perform calculations from outfiles
i = 0
fig = plt.figure(figsize=(8,8))
outfile_1C = '/home/brandonsmith/modeloutput/Control/1.0/map_polar_clouds_and_ice_seasonal_Control_1.0.nc'
outfile_1D = '/home/brandonsmith/modeloutput/2xCO2/1.0/map_polar_clouds_and_ice_seasonal_2xCO2_1.0.nc'
outfile_1Q = '/home/brandonsmith/modeloutput/4xCO2/1.0/map_polar_clouds_and_ice_seasonal_4xCO2_1.0.nc'
outfile_2xco2 = '/map_polar_clouds_and_ice_seasonal_2xCO2_1.0.nc'
outfile_control = '/map_polar_clouds_and_ice_seasonal_Control_1.0.nc'
outfile_4xco2 = '/map_polar_clouds_and_ice_seasonal_4xCO2_1.0.nc'
outpath_control = '/home/brandonsmith/modeloutput/Control/1.0'
outpath_2xco2 = '/home/brandonsmith/modeloutput/2xCO2/1.0'
outpath_4xco2 = '/home/brandonsmith/modeloutput/4xCO2/1.0'
dsloc_control = outpath_control+outfile_control
dsloc_2xco2 = outpath_2xco2+outfile_2xco2
dsloc_4xco2 = outpath_4xco2+outfile_4xco2
if os.path.isfile(dsloc_control) and os.path.isfile(dsloc_2xco2):
    dsc = netCDF4.Dataset(dsloc_control)
    ts = dsc.variables[TS][:,:,:]
    icefrac = dsc.variables[ICEFRAC][:,:,:]
    qflx = dsc.variables[QFLX][:,:,:]
    lowcloud = dsc.variables[CLDLOW][:,:,:]
    flds = dsc.variables[FLDS][:,:,:]
    lat = dsc.variables['lat'][:]
    lon = dsc.variables['lon'][:]
    #lev = dsc.variables['lev'][:]
    dsc.close()

    dsd = netCDF4.Dataset(dsloc_2xco2)
    tsd = dsd.variables[TS][:,:,:]
    icefracd = dsd.variables[ICEFRAC][:,:,:]
    qflxd = dsd.variables[QFLX][:,:,:]
    lowcloudd = dsd.variables[CLDLOW][:,:,:]
    fldsd = dsd.variables[FLDS][:,:,:]
    #lev = dsc.variables['lev'][:]
    dsd.close()

    dsq = netCDF4.Dataset(dsloc_4xco2)
    tsq = dsq.variables[TS][:,:,:]
    icefracq = dsq.variables[ICEFRAC][:,:,:]
    qflxq = dsq.variables[QFLX][:,:,:]
    lowcloudq = dsq.variables[CLDLOW][:,:,:]
    fldsq = dsq.variables[FLDS][:,:,:]
    #lev = dsc.variables['lev'][:]
    dsq.close()
    
    C1 = netCDF4.Dataset(outfile_1C)
    ts1 = C1.variables[TS][:,:,:]
    icefrac1 = C1.variables[ICEFRAC][:,:,:]
    qflx1 = C1.variables[QFLX][:,:,:]
    lowcloud1 = C1.variables[CLDLOW][:,:,:]
    C1.close()
    D1 = netCDF4.Dataset(outfile_1D)
    ts1d = D1.variables[TS][:,:,:]
    icefrac1d = D1.variables[ICEFRAC][:,:,:]
    qflx1d = D1.variables[QFLX][:,:,:]
    lowcloud1d = D1.variables[CLDLOW][:,:,:]
    D1.close()

    #Cvar = Cvar - Cvar2
    #Cvar1 = Cvar1 - Cvar1_2
    #Dvar = Dvar - Dvar2
    #Qvar = Qvar - Qvar2
    icefrac = np.squeeze(icefrac)
    icefracd = np.squeeze(icefracd)
    icefracDR = icefracd - icefrac
    Lon = lon
    ts = np.squeeze(ts)
    ts, Lon = add_cyclic_point(ts, coord=lon)
    Lon = lon
    tsd = np.squeeze(tsd)
    tsd, Lon = add_cyclic_point(tsd, coord=lon)
    Lon = lon
    tsq = np.squeeze(tsq)
    tsq, Lon = add_cyclic_point(tsq, coord=lon)
    Lon = lon
    ts1 = np.squeeze(ts1)
    ts1, Lon = add_cyclic_point(ts1, coord=lon)
    Lon = lon
    ts1d = np.squeeze(ts1d)
    ts1d, Lon = add_cyclic_point(ts1d, coord=lon)
    Lon = lon
    flds, Lon = add_cyclic_point(flds, coord=lon)
    Lon = lon
    fldsd, Lon = add_cyclic_point(fldsd, coord=lon)
    Lon = lon
    fldsq, Lon = add_cyclic_point(fldsq, coord=lon)
    Lon = lon
    icefrac, Lon = add_cyclic_point(icefrac, coord=lon)
    Lon = lon
    qflx = np.squeeze(qflx)
    qflx, Lon = add_cyclic_point(qflx, coord=lon)
    Lon = lon
    lowcloud = np.squeeze(lowcloud)
    lowcloud, Lon = add_cyclic_point(lowcloud, coord=lon)
    Lon = lon
    icefracd = np.squeeze(icefracd)
    icefracd, Lon = add_cyclic_point(icefracd, coord=lon)
    #Lon = lon
    #icefracq = np.squeeze(icefracq)
    #icefracq, Lon = add_cyclic_point(icefracq, coord=lon)
    Lon = lon
    qflxd = np.squeeze(qflxd)
    qflxd, Lon = add_cyclic_point(qflxd, coord=lon)
    Lon = lon
    qflxq = np.squeeze(qflxq)
    qflxq, Lon = add_cyclic_point(qflxq, coord=lon)
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
    qflx1 = np.squeeze(qflx1)
    qflx1, Lon = add_cyclic_point(qflx1, coord=lon)
    Lon = lon
    lowcloud1 = np.squeeze(lowcloud1)
    lowcloud1, Lon = add_cyclic_point(lowcloud1, coord=lon)
    Lon = lon
    icefrac1d = np.squeeze(icefrac1d)
    icefrac1d, Lon = add_cyclic_point(icefrac1d, coord=lon)
    Lon = lon
    qflx1d = np.squeeze(qflx1d)
    qflx1d, Lon = add_cyclic_point(qflx1d, coord=lon)
    Lon = lon
    lowcloud1d = np.squeeze(lowcloud1d)
    lowcloud1d, Lon = add_cyclic_point(lowcloud1d, coord=lon)
    
    tsDR = tsd - ts
    
    qflxDR = qflxd - qflx
    lowcloudDR = lowcloudd - lowcloud
    fldsDR = fldsd - flds
    
    print(np.average(icefracDR[1,:,:]))
    
Seasons = ['Winter, Spring, Summer, Autumn']
for i in range(0,4):
    #plot variable
    print(i)
    ax = fig.add_subplot(4,5,5*i+3,projection=ccrs.NorthPolarStereo())
    ax.set_extent([-180,180,60,90],ccrs.PlateCarree())
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)

    cs = ax.contourf(Lon,lat,lowcloudDR[i],12,transform=ccrs.PlateCarree(),cmap='BrBG',vmin=-0.4,vmax=0.4)
    ax.coastlines()
    #ax.gridlines(draw_labels=True,color='gray',linewidth=0.25,linestyle='dotted')

    if i == 0:
        cticks=np.around(np.linspace(-0.2,0.3,6),decimals=1)
        cbar_ax = fig.add_axes([0.41,-0.05,0.18,0.02])
        #cbar_ax = fig.add_axes([0,-0.05,0.23,0.02])
        #cb = fig.colorbar(mappable=None, norm=Normalize(vmin=-0.25,vmax=0.25), cmap='BrBG',spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
        cb = fig.colorbar(cs, spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
        #cb.formatter.set_powerlimits((0,0))
        cb.set_ticklabels(cticks.astype(str),fontsize=7)
        cb.set_label('Low Cloud Fraction Response',fontsize=7)

    ax2 = fig.add_subplot(4,5,5*i+2,projection=ccrs.NorthPolarStereo())
    ax2.set_extent([-180,180,60,90],ccrs.PlateCarree())
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax2.set_boundary(circle, transform=ax2.transAxes)

    cs2 = ax2.contourf(lon,lat,icefracDR[i],10,transform=ccrs.PlateCarree(),cmap='RdBu_r',vmin=-1,vmax=1)
    ax2.coastlines()
    #ax2.gridlines(draw_labels=True,color='gray',linewidth=0.25,linestyle='dotted')

    if i == 0:
        cticks=np.around(np.linspace(-1,0,6),decimals=1)
        cbar_ax = fig.add_axes([0.21,-0.05,0.18,0.02])
        cb2 = fig.colorbar(cs2, spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
        #cb2 = fig.colorbar(mappable=None, norm=Normalize(vmin=-0.8,vmax=0.8), cmap='RdBu_r',spacing='uniform',orientation='horizontal',cax=cbar_ax,ticks=cticks)
        cb2.set_label('Sea Ice Response ($W/m^2$)',fontsize=7)
        cb2.set_ticklabels(cticks.astype(str),fontsize=7)
   # cb2.formatter.set_powerlimits((0,0))

    ax3 = fig.add_subplot(4,5,5*i+1,projection=ccrs.NorthPolarStereo())
    ax3.set_extent([-180,180,60,90],ccrs.PlateCarree())

    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax3.set_boundary(circle, transform=ax3.transAxes)
    cs3 = ax3.contourf(Lon,lat,icefrac[i],10,transform=ccrs.PlateCarree(),cmap='Blues_r',vmin=0,vmax=1)
    ax3.coastlines()
    #ax3.gridlines(draw_labels=True,color='gray',linewidth=0.25,linestyle='dotted')


    if i == 0:
        cticks=np.around(np.linspace(0,1,6),decimals=1)
        cbar_ax = fig.add_axes([0,-0.05,0.18,0.02])
        #cbar_ax = fig.add_axes([0.51,-0.05,0.23,0.02])
        #cb3 = fig.colorbar(mappable=None, norm=Normalize(vmin=0,vmax=1), cmap='Blues_r',spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
        cb3 = plt.colorbar(cs3,spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
        cb3.set_label('Sea Ice Fraction (Reference Climates)',fontsize=7)
        cb3.set_ticklabels(cticks.astype(str),fontsize=7)

    ax4 = fig.add_subplot(4,5,5*i+4,projection=ccrs.NorthPolarStereo())
    ax4.set_extent([-180,180,60,90],ccrs.PlateCarree())
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax4.set_boundary(circle, transform=ax4.transAxes)
    cs4 = ax4.contourf(Lon,lat,qflxDR[i],100,transform=ccrs.PlateCarree(),cmap='RdBu_r',vmin=-50,vmax=50)
    ax4.coastlines()
    #ax4.gridlines(draw_labels=True,color='gray',linewidth=0.25,linestyle='dotted')

    if i == 0:
        cticks=np.around(np.linspace(-50,50,11),decimals=0)
        cbar_ax = fig.add_axes([0.61,-0.05,0.18,0.02])
        cb4 = fig.colorbar(mappable=None, norm=Normalize(vmin=-50,vmax=50), cmap='RdBu_r',spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
        #cb4 = fig.colorbar(cs4,spacing='proportional',orientation='horizontal',cax=cbar_ax)
        cb4.set_label('$\Delta$ Upwelling TOA Solar Flux ($W/m^2$)',fontsize=6)
        cb4.set_ticklabels(cticks.astype(int).astype(str),fontsize=7)
        #cb4.formatter.set_useOffset(True)
        
    ax5 = fig.add_subplot(4,5,5*i+5,projection=ccrs.NorthPolarStereo())
    ax5.set_extent([-180,180,60,90],ccrs.PlateCarree())
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax5.set_boundary(circle, transform=ax5.transAxes)
    cs5 = ax5.contourf(Lon,lat,fldsDR[i],100,transform=ccrs.PlateCarree(),cmap='RdBu_r',vmin=-50,vmax=50)
    ax5.coastlines()
    #ax4.gridlines(draw_labels=True,color='gray',linewidth=0.25,linestyle='dotted')

    if i == 0:
        cticks=np.around(np.linspace(-50,50,11),decimals=0)
        cbar_ax = fig.add_axes([0.81,-0.05,0.18,0.02])
        cb4 = fig.colorbar(mappable=None, norm=Normalize(vmin=-50,vmax=50), cmap='RdBu_r',spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
        #cb4 = fig.colorbar(cs4,spacing='proportional',orientation='horizontal',cax=cbar_ax)
        cb4.set_label('$\Delta$ Longwave Cloud Forcing ($W/m^2$)',fontsize=6)
        cb4.set_ticklabels(cticks.astype(int).astype(str),fontsize=7)
        #cb4.formatter.set_useOffset(True)
    
#    cb.set_label('Temperature (K)')
plt.suptitle('Seasonal Arctic Cloud and Sea Ice Response to $CO_2$ Doubling',y=1.01,fontsize=16)
fig.text(-0.02, 0.85, 'Winter', va='center', ha='left', rotation='vertical',fontsize=8)
fig.text(-0.02, 0.61, 'Spring', va='center', ha='left', rotation='vertical',fontsize=8)
fig.text(-0.02, 0.37, 'Summer', va='center', ha='left', rotation='vertical',fontsize=8)
fig.text(-0.02, 0.13, 'Autumn', va='center', ha='left', rotation='vertical',fontsize=8)
fig.tight_layout(pad=0.2)
fig.subplots_adjust(hspace=0.05,wspace=0.05)
plt.show()
fig.savefig(figure_path+'Map_seasonal_Polar_Clouds_and_Ice_2x.pdf',bbox_inches='tight')