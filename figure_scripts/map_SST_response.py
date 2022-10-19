#!/usr/bin/env python3

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

field = 'SST'
field2 = 'TS'
Zonal_Mean = False
Global_Mean = False
normalize = False
stf = False
if field == 'SST' and Zonal_Mean is False and Global_Mean is False:
    for CASENAME in casenames:
        run = 'Control'
        inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        infile = '/Merged_'+run+'_'+CASENAME+'.nc'
        outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        outfile_control = '/map_SST_'+run+'_'+CASENAME+'.nc'
        if os.path.isfile(inpath+infile):
            if not os.path.isfile(outpath+outfile_control):
                syscall = '/usr/bin/cdo timmean -seltimestep,-360/-1 -select,name='+field+','+field2+' '+inpath+infile+ ' '+outpath+outfile_control
                print(syscall)
        run = '2xCO2'
        inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        infile = '/Merged_'+run+'_'+CASENAME+'_.nc'
        outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        outfile_control = '/map_SST_'+run+'_'+CASENAME+'.nc'
        if os.path.isfile(inpath+infile):
            if not os.path.isfile(outpath+outfile_control):
                syscall = '/usr/bin/cdo timmean -seltimestep,-360/-1 -select,name='+field+','+field2+' '+inpath+infile+ ' '+outpath+outfile_control
                print(syscall)
        run = '4xCO2'
        inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        infile = '/Merged_'+run+'_'+CASENAME+'_.nc'
        outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        outfile_control = '/map_SST_'+run+'_'+CASENAME+'.nc'
        if os.path.isfile(inpath+infile):
            if not os.path.isfile(outpath+outfile_control):
                syscall = '/usr/bin/cdo timmean -seltimestep,-360/-1 -select,name='+field+','+field2+' '+inpath+infile+ ' '+outpath+outfile_control
                print(syscall)
    # Load variables and perform calculations from outfiles
    i = 0
    fig = plt.figure(figsize=(20,10))
    outfile_1C = '/home/brandonsmith/modeloutput/Control/1.0/map_SST_Control_1.0.nc'
    outfile_1D = '/home/brandonsmith/modeloutput/2xCO2/1.0/map_SST_2xCO2_1.0.nc'
    outfile_1Q = '/home/brandonsmith/modeloutput/4xCO2/1.0/map_SST_4xCO2_1.0.nc'
    rf = [5.53394852, 4.72283731, 3.88779531, 2.78457941]
    rfq = [11.4426504, 9.82042797, 8.15034396, 5.94391217]
    for CASENAME in casenames:
        outfile_2xco2 = '/map_SST_2xCO2_'+CASENAME+'.nc'
        outfile_control = '/map_SST_Control_'+CASENAME+'.nc'
        outfile_4xco2 = '/map_SST_4xCO2_'+CASENAME+'.nc'
        outpath_control = '/home/brandonsmith/modeloutput/Control/'+CASENAME
        outpath_2xco2 = '/home/brandonsmith/modeloutput/2xCO2/'+CASENAME
        outpath_4xco2 = '/home/brandonsmith/modeloutput/4xCO2/'+CASENAME
        dsloc_control = outpath_control+outfile_control
        dsloc_2xco2 = outpath_2xco2+outfile_2xco2
        dsloc_4xco2 = outpath_4xco2+outfile_4xco2
        if os.path.isfile(dsloc_control) and os.path.isfile(dsloc_2xco2):
            
            dsc = netCDF4.Dataset(dsloc_control)
            Cvar = dsc.variables[field][:]
            lat = dsc.variables['lat'][:]
            lon = dsc.variables['lon'][:]
            #lev = dsc.variables['lev'][:]
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
            

            '''
            print("Loading Variables")
            dsc = netCDF4.Dataset(dsloc_control)
            U = dsc.variables['U'][:]
            V = dsc.variables['V'][:]
            UU = dsc.variables['UU'][:]
            VV = dsc.variables['VV'][:]
            lat = dsc.variables['lat'][:]
            lon = dsc.variables['lon'][:]
            #lev = dsc.variables['lev'][:]
            dsc.close()
            
            dsd = netCDF4.Dataset(dsloc_2xco2)
            Ud = dsd.variables['U'][:]
            Vd = dsd.variables['V'][:]
            UUd = dsd.variables['UU'][:]
            VVd = dsd.variables['VV'][:]
            # no need to reload lat/lon
            dsd.close()
            
            C1 = netCDF4.Dataset(outfile_1C)
            U1 = C1.variables['U'][:]
            V1 = C1.variables['V'][:]
            UU1 = C1.variables['UU'][:]
            VV1 = C1.variables['VV'][:]
            C1.close()
            D1 = netCDF4.Dataset(outfile_1D)
            U1d = D1.variables['U'][:]
            V1d = D1.variables['V'][:]
            UU1d = D1.variables['UU'][:]
            VV1d = D1.variables['VV'][:]
            D1.close()
            
            dsq = netCDF4.Dataset(dsloc_4xco2)
            Uq = dsq.variables['U'][:]
            Vq = dsq.variables['V'][:]
            UUq = dsq.variables['UU'][:]
            VVq = dsq.variables['VV'][:]
            dsq.close()
            
            Q1 = netCDF4.Dataset(outfile_1Q)
            U1q = Q1.variables['U'][:]
            V1q = Q1.variables['V'][:]
            UU1q = Q1.variables['UU'][:]
            VV1q = Q1.variables['VV'][:]
            Q1.close()
            
            print('computing eddy components')
            eddy_VV1 = VV1 - V1**2
            eddy_UU1 = UU1 - U1**2
            eddy_VV1d = VV1d - V1d**2
            eddy_UU1d = UU1d - U1d**2
            eddy_VV1q = VV1q - V1q**2
            eddy_UU1q = UU1q - U1q**2
            
            eddy_VV = VV - V**2
            eddy_UU = UU - U**2
            eddy_VVd = VVd - Vd**2
            eddy_UUd = UUd - Ud**2
            eddy_VVq = VVq - Vq**2
            eddy_UUq = UUq - Uq**2
            
            print('computing EKE')
            EKE = 0.5*(eddy_VV+eddy_UU)
            EKEd = 0.5*(eddy_VVd+eddy_UUd)
            EKEq = 0.5*(eddy_VVq+eddy_UUq)
            EKE1 = 0.5*(eddy_VV1+eddy_UU1)
            EKE1d = 0.5*(eddy_VV1d+eddy_UU1d)
            EKE1q = 0.5*(eddy_VV1q+eddy_UU1q)
            
            print('computing averages')
            EKE1 = np.average(EKE1,0)
            EKE1 = np.average(EKE1,0)
            EKE1d = np.average(EKE1d,0)
            EKE1d = np.average(EKE1d,0)
            EKE1q = np.average(EKE1q,0)
            EKE1q = np.average(EKE1q,0)
            EKE = np.average(EKE,0)
            EKE = np.average(EKE,0)
            EKEd = np.average(EKEd,0)
            EKEd = np.average(EKEd,0)
            EKEq = np.average(EKEq,0)
            EKEq = np.average(EKEq,0)
            '''
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
            cs = ax.contourf(Lon,lat,Cvar,181,transform=ccrs.PlateCarree(),cmap='gist_heat',vmin=273,vmax=303)
            #ax.coastlines()
            ax.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k')
            #ax.set_ylabel('Solar Multiplier: '+CASENAME,fontsize=14)
            if i == 0:
                plt.title('Base Climate Comparison',fontsize=18)
            if i == 0:
                cticks=np.around(np.linspace(273,303,5),decimals=1)
                cbar_ax = fig.add_axes([0.02,-0.05,0.15,0.02])
                cb = fig.colorbar(mappable=None, norm=Normalize(vmin=273,vmax=303), cmap='gist_heat',spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
                #cb.formatter.set_powerlimits((0,0))
                #cb.set_ticklabels(cticks.astype(str),fontsize=12)
                cb.set_ticklabels(cticks.astype(int).astype(str),fontsize=12)
                cb.set_label('Sea Surface Temperature (K)',fontsize=14)
            ax.text(0.03, 0.5, Casenames[i], va='center', rotation='horizontal',transform=ax.transAxes,fontsize=14,color='black')
            ax2 = fig.add_subplot(4,5,5*i+2,projection=ccrs.Robinson())
            cs2 = ax2.contourf(Lon,lat,doubling_response,201,transform=ccrs.PlateCarree(),cmap='Reds',vmin=0,vmax=15)
            #ax2.coastlines()
            ax2.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k')
            if i == 0:
                plt.title('Response to $CO_2$ Doubling',fontsize=18)
                #plt.title('Response to $F_2$$_x$',fontsize=18)
            if i == 1:
                cticks=np.around(np.linspace(0,15,6),decimals=2)
                cbar_ax = fig.add_axes([0.22,-0.05,0.15,0.02])
                #cb2 = fig.colorbar(cs2, spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
                cb2 = fig.colorbar(mappable=None, norm=Normalize(vmin=0,vmax=15), cmap='Reds',spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
                cb2.set_label('Change in SST (K)',fontsize=14)
                #cb2.set_ticklabels(cticks.astype(str),fontsize=12)
                cb2.set_ticklabels(cticks.astype(int).astype(str),fontsize=12)
               # cb2.formatter.set_powerlimits((0,0))
            ax3 = fig.add_subplot(4,5,5*i+3,projection=ccrs.Robinson())
            cs3 = ax3.contourf(Lon,lat,diff,201,transform=ccrs.PlateCarree(),cmap='RdBu_r',vmin=-15,vmax=15)
            #ax3.coastlines()
            ax3.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k')
            if i == 0:
                plt.title('Difference w.r.t. $S_0$',fontsize=18)
            if i == 0:
                cticks=np.around(np.linspace(-15,15,11),decimals=2)
                cbar_ax = fig.add_axes([0.42,-0.05,0.15,0.02])
                cb3 = fig.colorbar(mappable=None, norm=Normalize(vmin=-15,vmax=15), cmap='RdBu_r',spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
                #cb3 = fig.colorbar(cs3,spacing='uniform',orientation='horizontal',cax=cbar_ax,ticks=cticks)
                cb3.set_label('Difference',fontsize=14)
                #cb3.formatter.set_powerlimits((0,0))
                #cb3.set_ticklabels(cticks.astype(str),fontsize=12)
                cb3.set_ticklabels(cticks.astype(int).astype(str),fontsize=12)
            ax4 = fig.add_subplot(4,5,5*i+4,projection=ccrs.Robinson())
            cs4 = ax4.contourf(Lon,lat,Quad_response,201,transform=ccrs.PlateCarree(),cmap='Reds',vmin=0,vmax=15)
            #ax4.coastlines()
            ax4.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k')
            if i == 0:
                plt.title('Response to $CO_2$ Quadrupling',fontsize=18)
                #plt.title('Response to $F_4$$_x$',fontsize=18)
            if i == 0:
                cticks=np.around(np.linspace(0,15,6),decimals=2)
                cbar_ax = fig.add_axes([0.62,-0.05,0.15,0.02])
                cb4 = fig.colorbar(mappable=None, norm=Normalize(vmin=0,vmax=15), cmap='Reds',spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
                #cb3 = fig.colorbar(cs3,spacing='uniform',orientation='horizontal',cax=cbar_ax,ticks=cticks)
                cb4.set_label('Change in SST (K)',fontsize=14)
                #cb3.formatter.set_powerlimits((0,0))
                #cb4.set_ticklabels(cticks.astype(str),fontsize=12)
                cb4.set_ticklabels(cticks.astype(int).astype(str),fontsize=12)
            ax5 = fig.add_subplot(4,5,5*i+5,projection=ccrs.Robinson())
            cs5 = ax5.contourf(Lon,lat,Qdiff,201,transform=ccrs.PlateCarree(),cmap='RdBu_r',vmin=-15,vmax=15)
            #ax5.coastlines()
            ax5.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k')
            if i == 0:
                plt.title('Difference w.r.t. $S_0$',fontsize=18)
            if i == 0:
                cticks=np.around(np.linspace(-15,15,11),decimals=2)
                cbar_ax = fig.add_axes([0.82,-0.05,0.15,0.02])
                cb5 = fig.colorbar(mappable=None, norm=Normalize(vmin=-15,vmax=15), cmap='RdBu_r',spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
                #cb3 = fig.colorbar(cs3,spacing='uniform',orientation='horizontal',cax=cbar_ax,ticks=cticks)
                cb5.set_label('Difference',fontsize=14)
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
    #plt.suptitle('Sea Surface Temperature Response to Equivalent Radiative Forcing',fontsize=24)
    plt.suptitle('Sea Surface Temperature Response to $CO_2$ Doublings',fontsize=24)
    fig.tight_layout(pad=0.2)
    plt.show()
    fig.savefig(figure_path+'map_SST_corrected_2x_response.pdf',bbox_inches='tight')
    '''
    fig2 = plt.figure(2)
    outfile_1C = '/home/brandonsmith/modeloutput/Control/1.0/map_CLDLOW_Control_1.0.nc'
    if os.path.isfile(outfile_1C):
        ds = netCDF4.Dataset(outfile_1C)
        var = ds.variables['CLDLOW'][:]
        lon = ds.variables['lon'][:]
        lat = ds.variables['lat'][:]
        ds.close()
        var = var.squeeze()
    var, lon = add_cyclic_point(var, coord=lon)
    print(type(var))
    ax = plt.axes(projection=ccrs.Robinson())
    cs = ax.contourf(lon,lat,var,60,transform=ccrs.PlateCarree(),cmap='bone',vmin=-100,vmax=0)
    ax.coastlines()
    #ax.set_ylabel('Solar Multiplier: '+CASENAME)
    ax.set_title('Shortwave Cloud Forcing, Present Day Climate',fontsize=14)
    cticks=np.linspace(220,320,11)
    cbar_ax = fig2.add_axes([0.95,0.2,0.02,0.6])
    cb = fig2.colorbar(cs, spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=cticks)
    cb.set_label('Forcing ($W/m^2$)',fontsize=8)
    fig2.savefig(figure_path+'Control_CLDLOW.pdf',bbox_inches='tight')
    '''
    