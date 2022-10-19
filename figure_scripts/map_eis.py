##### !/usr/bin/env python3

import os
import time
import sys
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point

figure_path = '/home/brandonsmith/climate-gcm-bps/plots/'
casenames = ['0.9','0.95','1.0','1.05']

field = 'T'
field2 = 'T'
Zonal_Mean = False
Global_Mean = False
normalize = False
stf = False
if field == 'T' and Zonal_Mean is False and Global_Mean is False:
    for CASENAME in casenames:
        run = 'Control'
        inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        infile = '/Merged_'+run+'_'+CASENAME+'.nc'
        outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        T700 = '/map_lts_700Temp_'+run+'_'+CASENAME+'.nc'
        T1000 = '/map_lts_1000Temp_'+run+'_'+CASENAME+'.nc'
        outfilesub = '/map_lts_sub_'+run+'_'+CASENAME+'.nc'
        outfile700 = '/map_lts_700_'+run+'_'+CASENAME+'.nc'
        outfile1000 = '/map_lts_1000_'+run+'_'+CASENAME+'.nc'
        q850 = '/map_q850_'+run+'_'+CASENAME+'.nc'
        qs850 = '/map_qs850_'+run+'_'+CASENAME+'.nc'
        tempsum = '/map_tempsum_'+run+'_'+CASENAME+'.nc'
        z700 = '/map_z700_'+run+'_'+CASENAME+'.nc'
        LCL = '/map_lcl_'+run+'_'+CASENAME+'.nc'
        outfile_LCL = '/map_LCL_'+run+'_'+CASENAME+'.nc'
        outfile_eis = '/map_eis_'+run+'_'+CASENAME+'.nc'
        if os.path.isfile(inpath+infile):
            if not os.path.isfile(outpath+T700):
                syscall = '/usr/bin/cdo seltimestep,-360/-1 -sellevidx,21 -select,name=T '+inpath+infile+ ' '+outpath+T700
                print(syscall)
            if not os.path.isfile(outpath+outfile700):
                syscall = '/usr/bin/cdo timmean -select,name=lts -expr,\'lts=(T*(1000/691.389430314302)^0.286)\' '+outpath+T700+' '+outpath+outfile700
                print(syscall)
            if not os.path.isfile(outpath+T1000):
                syscall = '/usr/bin/cdo seltimestep,-360/-1 -select,name=TREFHT,PS '+inpath+infile+ ' '+outpath+T1000
                print(syscall)
            if not os.path.isfile(outpath+outfile1000):
                syscall = '/usr/bin/cdo timmean -select,name=lts -expr,\'lts=(TREFHT*(1000/(PS*0.01))^0.286)\' '+inpath+infile+' '+outpath+outfile1000
                print(syscall)
            if not os.path.isfile(outpath+outfilesub):
                syscall = '/usr/bin/cdo sub '+outpath+outfile700+' '+outpath+outfile1000 + ' '+outpath+outfilesub
                print(syscall)
            if not os.path.isfile(outpath+q850):
                syscall = '/usr/bin/cdo -seltimestep,-360/-1 -sellevidx,24 -select,name=Q,RELHUM '+inpath+infile+' '+outpath+q850
                print(syscall)
            if not os.path.isfile(outpath+qs850):
                syscall = '/usr/bin/cdo timmean -select,name=smr -expr,\'smr=(Q/RELHUM)\'  '+outpath+q850+' '+outpath+qs850
                print(syscall)
                syscall = 'ncrename -d lev,level -v lev,lev850 '+outpath+qs850
                print(syscall)
            if not os.path.isfile(outpath+z700):
                syscall = '/usr/bin/cdo timmean -seltimestep,-360/-1 -sellevidx,21 -select,name=Z3 '+inpath+infile+' '+outpath+z700
                print(syscall)
                syscall = 'ncrename -d lev,level -v lev,lev700 '+outpath+z700
                print(syscall)
            if not os.path.isfile(outpath+tempsum):
                syscall = '/usr/bin/cdo timmean -add '+outpath+T700+' '+outpath+T1000+' '+outpath+tempsum
                print(syscall)
            if not os.path.isfile(outpath+LCL):
                syscall = '/usr/bin/cdo -seltimestep,-360/-1 -select,name=TREFHT,RELHUM '+inpath+infile+' '+outpath+LCL
                print(syscall)
            if not os.path.isfile(outpath+outfile_LCL):
                syscall = '/usr/bin/cdo timmean -select,name=LCL -sellevidx,30 -expr,\'LCL=((20+(TREFHT-273.15)/5)*(100-RELHUM))\' '+outpath+LCL+' '+outpath+outfile_LCL
                print(syscall)
                syscall = 'ncrename -d lev,level -v lev,lev1000 '+outpath+outfile_LCL
                print(syscall)
            if not os.path.isfile(outpath+outfile_eis):
                syscall = '/usr/bin/cdo merge '+outpath+qs850+' '+outpath+z700+' '+outpath+tempsum+' '+outpath+outfile_LCL+' '+outpath+outfilesub+' '+outpath+outfile_eis
                print(syscall)
        
                
                
        run = '2xCO2'
        inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        infile = '/Merged_'+run+'_'+CASENAME+'.nc'
        outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        T700 = '/map_lts_700Temp_'+run+'_'+CASENAME+'.nc'
        T1000 = '/map_lts_1000Temp_'+run+'_'+CASENAME+'.nc'
        outfilesub = '/map_lts_sub_'+run+'_'+CASENAME+'.nc'
        outfile700 = '/map_lts_700_'+run+'_'+CASENAME+'.nc'
        outfile1000 = '/map_lts_1000_'+run+'_'+CASENAME+'.nc'
        q850 = '/map_q850_'+run+'_'+CASENAME+'.nc'
        qs850 = '/map_qs850_'+run+'_'+CASENAME+'.nc'
        tempsum = '/map_tempsum_'+run+'_'+CASENAME+'.nc'
        z700 = '/map_z700_'+run+'_'+CASENAME+'.nc'
        LCL = '/map_lcl_'+run+'_'+CASENAME+'.nc'
        outfile_LCL = '/map_LCL_'+run+'_'+CASENAME+'.nc'
        outfile_eis = '/map_eis_'+run+'_'+CASENAME+'.nc'
        if os.path.isfile(inpath+infile):
            if not os.path.isfile(outpath+T700):
                syscall = '/usr/bin/cdo seltimestep,-360/-1 -sellevidx,21 -select,name=T '+inpath+infile+ ' '+outpath+T700
                print(syscall)
            if not os.path.isfile(outpath+outfile700):
                syscall = '/usr/bin/cdo timmean -select,name=lts -expr,\'lts=(T*(1000/691.389430314302)^0.286)\' '+outpath+T700+' '+outpath+outfile700
                print(syscall)
            if not os.path.isfile(outpath+T1000):
                syscall = '/usr/bin/cdo seltimestep,-360/-1 -select,name=TREFHT,PS '+inpath+infile+ ' '+outpath+T1000
                print(syscall)
            if not os.path.isfile(outpath+outfile1000):
                syscall = '/usr/bin/cdo timmean -select,name=lts -expr,\'lts=(TREFHT*(1000/(PS*0.01))^0.286)\' '+inpath+infile+' '+outpath+outfile1000
                print(syscall)
            if not os.path.isfile(outpath+outfilesub):
                syscall = '/usr/bin/cdo sub '+outpath+outfile700+' '+outpath+outfile1000 + ' '+outpath+outfilesub
                print(syscall)
            if not os.path.isfile(outpath+q850):
                syscall = '/usr/bin/cdo -seltimestep,-360/-1 -sellevidx,24 -select,name=Q,RELHUM '+inpath+infile+' '+outpath+q850
                print(syscall)
            if not os.path.isfile(outpath+qs850):
                syscall = '/usr/bin/cdo timmean -select,name=smr -expr,\'smr=(Q/RELHUM)\'  '+outpath+q850+' '+outpath+qs850
                print(syscall)
                syscall = 'ncrename -d lev,level -v lev,lev850 '+outpath+qs850
                print(syscall)
            if not os.path.isfile(outpath+z700):
                syscall = '/usr/bin/cdo timmean -seltimestep,-360/-1 -sellevidx,21 -select,name=Z3 '+inpath+infile+' '+outpath+z700
                print(syscall)
                syscall = 'ncrename -d lev,level -v lev,lev700 '+outpath+z700
                print(syscall)
            if not os.path.isfile(outpath+tempsum):
                syscall = '/usr/bin/cdo timmean -add '+outpath+T700+' '+outpath+T1000+' '+outpath+tempsum
                print(syscall)
            if not os.path.isfile(outpath+LCL):
                syscall = '/usr/bin/cdo -seltimestep,-360/-1 -select,name=TREFHT,RELHUM '+inpath+infile+' '+outpath+LCL
                print(syscall)
            if not os.path.isfile(outpath+outfile_LCL):
                syscall = '/usr/bin/cdo timmean -select,name=LCL -sellevidx,30 -expr,\'LCL=((20+(TREFHT-273.15)/5)*(100-RELHUM))\' '+outpath+LCL+' '+outpath+outfile_LCL
                print(syscall)
                syscall = 'ncrename -d lev,level -v lev,lev1000 '+outpath+outfile_LCL
                print(syscall)
            if not os.path.isfile(outpath+outfile_eis):
                syscall = '/usr/bin/cdo merge '+outpath+qs850+' '+outpath+z700+' '+outpath+tempsum+' '+outpath+outfile_LCL+' '+outpath+outfilesub+' '+outpath+outfile_eis
                print(syscall)
                
        run = '4xCO2'
        inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        infile = '/Merged_'+run+'_'+CASENAME+'.nc'
        outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        T700 = '/map_lts_700Temp_'+run+'_'+CASENAME+'.nc'
        T1000 = '/map_lts_1000Temp_'+run+'_'+CASENAME+'.nc'
        outfilesub = '/map_lts_sub_'+run+'_'+CASENAME+'.nc'
        outfile700 = '/map_lts_700_'+run+'_'+CASENAME+'.nc'
        outfile1000 = '/map_lts_1000_'+run+'_'+CASENAME+'.nc'
        q850 = '/map_q850_'+run+'_'+CASENAME+'.nc'
        qs850 = '/map_qs850_'+run+'_'+CASENAME+'.nc'
        tempsum = '/map_tempsum_'+run+'_'+CASENAME+'.nc'
        z700 = '/map_z700_'+run+'_'+CASENAME+'.nc'
        LCL = '/map_lcl_'+run+'_'+CASENAME+'.nc'
        outfile_LCL = '/map_LCL_'+run+'_'+CASENAME+'.nc'
        outfile_eis = '/map_eis_'+run+'_'+CASENAME+'.nc'
        if os.path.isfile(inpath+infile):
            if not os.path.isfile(outpath+T700):
                syscall = '/usr/bin/cdo seltimestep,-360/-1 -sellevidx,21 -select,name=T '+inpath+infile+ ' '+outpath+T700
                print(syscall)
            if not os.path.isfile(outpath+outfile700):
                syscall = '/usr/bin/cdo timmean -select,name=lts -expr,\'lts=(T*(1000/691.389430314302)^0.286)\' '+outpath+T700+' '+outpath+outfile700
                print(syscall)
            if not os.path.isfile(outpath+T1000):
                syscall = '/usr/bin/cdo seltimestep,-360/-1 -select,name=TREFHT,PS '+inpath+infile+ ' '+outpath+T1000
                print(syscall)
            if not os.path.isfile(outpath+outfile1000):
                syscall = '/usr/bin/cdo timmean -select,name=lts -expr,\'lts=(TREFHT*(1000/(PS*0.01))^0.286)\' '+inpath+infile+' '+outpath+outfile1000
                print(syscall)
            if not os.path.isfile(outpath+outfilesub):
                syscall = '/usr/bin/cdo sub '+outpath+outfile700+' '+outpath+outfile1000 + ' '+outpath+outfilesub
                print(syscall)
            if not os.path.isfile(outpath+q850):
                syscall = '/usr/bin/cdo -seltimestep,-360/-1 -sellevidx,24 -select,name=Q,RELHUM '+inpath+infile+' '+outpath+q850
                print(syscall)
            if not os.path.isfile(outpath+qs850):
                syscall = '/usr/bin/cdo timmean -select,name=smr -expr,\'smr=(Q/RELHUM)\'  '+outpath+q850+' '+outpath+qs850
                print(syscall)
                syscall = 'ncrename -d lev,level -v lev,lev850 '+outpath+qs850
                print(syscall)
            if not os.path.isfile(outpath+z700):
                syscall = '/usr/bin/cdo timmean -seltimestep,-360/-1 -sellevidx,21 -select,name=Z3 '+inpath+infile+' '+outpath+z700
                print(syscall)
                syscall = 'ncrename -d lev,level -v lev,lev700 '+outpath+z700
                print(syscall)
            if not os.path.isfile(outpath+tempsum):
                syscall = '/usr/bin/cdo timmean -add '+outpath+T700+' '+outpath+T1000+' '+outpath+tempsum
                print(syscall)
            if not os.path.isfile(outpath+LCL):
                syscall = '/usr/bin/cdo -seltimestep,-360/-1 -select,name=TREFHT,RELHUM '+inpath+infile+' '+outpath+LCL
                print(syscall)
            if not os.path.isfile(outpath+outfile_LCL):
                syscall = '/usr/bin/cdo timmean -select,name=LCL -sellevidx,30 -expr,\'LCL=((20+(TREFHT-273.15)/5)*(100-RELHUM))\' '+outpath+LCL+' '+outpath+outfile_LCL
                print(syscall)
                syscall = 'ncrename -d lev,level -v lev,lev1000 '+outpath+outfile_LCL
                print(syscall)
            if not os.path.isfile(outpath+outfile_eis):
                syscall = '/usr/bin/cdo merge '+outpath+qs850+' '+outpath+z700+' '+outpath+tempsum+' '+outpath+outfile_LCL+' '+outpath+outfilesub+' '+outpath+outfile_eis
                print(syscall)
    
    # Load variables and perform calculations from outfiles
    i = 0
    fig = plt.figure(figsize=(20,10))
    outfile_1C = '/home/brandonsmith/modeloutput/Control/1.0/map_eis_Control_1.0.nc'
    outfile_1D = '/home/brandonsmith/modeloutput/2xCO2/1.0/map_eis_2xCO2_1.0.nc'
    outfile_1Q = '/home/brandonsmith/modeloutput/4xCO2/1.0/map_eis_4xCO2_1.0.nc'
    rf = [5.53394852, 4.72283731, 3.88779531, 2.78457941]
    rfq = [11.4426504, 9.82042797, 8.15034396, 5.94391217]
    for CASENAME in casenames:
        outfile_2xco2 = '/map_eis_2xCO2_'+CASENAME+'.nc'
        outfile_control = '/map_eis_Control_'+CASENAME+'.nc'
        outfile_4xco2 = '/map_eis_4xCO2_'+CASENAME+'.nc'
        outpath_control = '/home/brandonsmith/modeloutput/Control/'+CASENAME
        outpath_2xco2 = '/home/brandonsmith/modeloutput/2xCO2/'+CASENAME
        outpath_4xco2 = '/home/brandonsmith/modeloutput/4xCO2/'+CASENAME
        dsloc_control = outpath_control+outfile_control
        dsloc_2xco2 = outpath_2xco2+outfile_2xco2
        dsloc_4xco2 = outpath_4xco2+outfile_4xco2
        if os.path.isfile(dsloc_control) and os.path.isfile(dsloc_2xco2):
            
            dsc = netCDF4.Dataset(dsloc_control)
            lts= dsc.variables['lts'][:]
            tempsum = dsc.variables['TREFHT'][:]
            qs850 = dsc.variables['smr'][:]
            z700 = dsc.variables['Z3'][:]
            LCL = dsc.variables['LCL'][:]
            lat = dsc.variables['lat'][:]
            lon = dsc.variables['lon'][:]
            #lev = dsc.variables['lev'][:]
            dsc.close()
            
            dsd = netCDF4.Dataset(dsloc_2xco2)
            ltsd= dsd.variables['lts'][:]
            tempsumd = dsd.variables['TREFHT'][:]
            qs850d = dsd.variables['smr'][:]
            z700d = dsd.variables['Z3'][:]
            LCLd = dsd.variables['LCL'][:]
            # no need to reload lat/lon
            dsd.close()
            
            C1 = netCDF4.Dataset(outfile_1C)
            lts1= C1.variables['lts'][:]
            tempsum1 = C1.variables['TREFHT'][:]
            qs8501 = C1.variables['smr'][:]
            z7001 = C1.variables['Z3'][:]
            LCL1 = C1.variables['LCL'][:]
            C1.close()
            D1 = netCDF4.Dataset(outfile_1D)
            ltsd1= D1.variables['lts'][:]
            tempsumd1 = D1.variables['TREFHT'][:]
            qs850d1 = D1.variables['smr'][:]
            z700d1 = D1.variables['Z3'][:]
            LCLd1 = D1.variables['LCL'][:]
            D1.close()
            
            dsq = netCDF4.Dataset(dsloc_4xco2)
            ltsq= dsq.variables['lts'][:]
            tempsumq = dsq.variables['TREFHT'][:]
            qs850q = dsq.variables['smr'][:]
            z700q = dsq.variables['Z3'][:]
            LCLq = dsq.variables['LCL'][:]
            dsq.close()
            
            Q1 = netCDF4.Dataset(outfile_1Q)
            ltsq1= Q1.variables['lts'][:]
            tempsumq1 = Q1.variables['TREFHT'][:]
            qs850q1 = Q1.variables['smr'][:]
            z700q1 = Q1.variables['Z3'][:]
            LCLq1 = Q1.variables['LCL'][:]
            Q1.close()
            
            pal = (0.00989)*(1-(1+2450000*qs850/(287.058*(tempsum/2)))/(1+(2450000**2)*qs850/(993*461.4*((tempsum/2)**2))))	
            pald = (0.00989)*(1-(1+2450000*qs850d/(287.058*(tempsumd/2)))/(1+(2450000**2)*qs850d/(993*461.4*((tempsumd/2)**2))))	
            palq = (0.00989)*(1-(1+2450000*qs850q/(287.058*(tempsumq/2)))/(1+(2450000**2)*qs850q/(993*461.4*((tempsumq/2)**2))))
            pal1 = (0.00989)*(1-(1+2450000*qs8501/(287.058*(tempsum1/2)))/(1+(2450000**2)*qs8501/(993*461.4*((tempsum1/2)**2))))	
            pald1 = (0.00989)*(1-(1+2450000*qs850d1/(287.058*(tempsumd1/2)))/(1+(2450000**2)*qs850d1/(993*461.4*((tempsumd1/2)**2))))
            palq1 = (0.00989)*(1-(1+2450000*qs850q1/(287.058*(tempsumq1/2)))/(1+(2450000**2)*qs850q1/(993*461.4*((tempsumq1/2)**2))))	
            
            Cvar = lts - pal*(z700-LCL)
            Dvar = ltsd - pald*(z700d-LCLd)
            Qvar = ltsq - palq*(z700q-LCLq)
            Cvar1 = lts1 - pal1*(z7001-LCL1)
            Dvar1 = ltsd1 - pald1*(z700d1-LCLd1)
            Qvar1 = ltsq1 - palq1*(z700q1-LCLq1)
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
            print(np.shape(Cvar))
            ax = fig.add_subplot(4,5,5*i+1,projection=ccrs.Robinson())
            cs = ax.contourf(Lon,lat,Cvar,21,transform=ccrs.PlateCarree(),cmap='PRGn',vmin=-40,vmax=40)
            ax.coastlines()
            #ax.set_ylabel('Solar Multiplier: '+CASENAME,fontsize=14)
            if i == 0:
                plt.title('Base Climate Comparison',fontsize=14)
            if i == 0:
                cticks=np.around(np.linspace(-40,40,9),decimals=1)
                cbar_ax = fig.add_axes([0.02,-0.05,0.15,0.02])
                cb = fig.colorbar(cs, spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
                #cb.formatter.set_powerlimits((0,0))
                #cb.set_ticklabels(cticks.astype(str),fontsize=12)
                cb.set_ticklabels(cticks.astype(int).astype(str),fontsize=12)
                cb.set_label('EIS)',fontsize=14)
            ax.text(0.03, 0.5, Casenames[i], va='center', rotation='horizontal',transform=ax.transAxes,fontsize=14,color='black')
            ax2 = fig.add_subplot(4,5,5*i+2,projection=ccrs.Robinson())
            cs2 = ax2.contourf(Lon,lat,doubling_response,30,transform=ccrs.PlateCarree(),cmap='PRGn',vmin=-10,vmax=10)
            ax2.coastlines()
            if i == 0:
                plt.title('Response to $CO_2$ Doubling',fontsize=14)
                #plt.title('Response to $F_2$$_x$',fontsize=14)
                
            if i == 1:
                cticks=np.around(np.linspace(-10,10,11),decimals=1)
                cbar_ax = fig.add_axes([0.22,-0.05,0.15,0.02])
                #cb2 = fig.colorbar(cs2, spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
                cb2 = fig.colorbar(mappable=None, norm=Normalize(vmin=-10,vmax=10), cmap='PRGn',spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
                cb2.set_label('Change in EIS (K)',fontsize=14)
                #cb2.set_ticklabels(cticks.astype(str),fontsize=12)
                cb2.set_ticklabels(cticks.astype(int).astype(str),fontsize=12)
               # cb2.formatter.set_powerlimits((0,0))
            ax3 = fig.add_subplot(4,5,5*i+3,projection=ccrs.Robinson())
            cs3 = ax3.contourf(Lon,lat,diff,30,transform=ccrs.PlateCarree(),cmap='PRGn',vmin=-10,vmax=10)
            ax3.coastlines()
            if i == 0:
                plt.title('Difference w.r.t. $S_0$',fontsize=14)
            if i == 0:
                cticks=np.around(np.linspace(-10,10,11),decimals=1)
                cbar_ax = fig.add_axes([0.42,-0.05,0.15,0.02])
                cb3 = fig.colorbar(mappable=None, norm=Normalize(vmin=-10,vmax=10), cmap='PRGn',spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
                #cb3 = fig.colorbar(cs3,spacing='uniform',orientation='horizontal',cax=cbar_ax,ticks=cticks)
                cb3.set_label('Difference (K)',fontsize=14)
                #cb3.formatter.set_powerlimits((0,0))
                #cb3.set_ticklabels(cticks.astype(str),fontsize=12)
                cb3.set_ticklabels(cticks.astype(int).astype(str),fontsize=12)
            ax4 = fig.add_subplot(4,5,5*i+4,projection=ccrs.Robinson())
            cs4 = ax4.contourf(Lon,lat,Quad_response,30,transform=ccrs.PlateCarree(),cmap='PRGn',vmin=-10,vmax=10)
            ax4.coastlines()
            if i == 0:
                plt.title('Response to $CO_2$ Quadrupling',fontsize=14)
                #plt.title('Response to $F_4$$_x$',fontsize=14)
            if i == 0:
                cticks=np.around(np.linspace(-10,10,11),decimals=1)
                cbar_ax = fig.add_axes([0.62,-0.05,0.15,0.02])
                cb4 = fig.colorbar(mappable=None, norm=Normalize(vmin=-10,vmax=10), cmap='PRGn',spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
                #cb3 = fig.colorbar(cs3,spacing='uniform',orientation='horizontal',cax=cbar_ax,ticks=cticks)
                cb4.set_label('Change in EIS (K)',fontsize=14)
                #cb3.formatter.set_powerlimits((0,0))
                #cb4.set_ticklabels(cticks.astype(str),fontsize=12)
                cb4.set_ticklabels(cticks.astype(int).astype(str),fontsize=12)
            ax5 = fig.add_subplot(4,5,5*i+5,projection=ccrs.Robinson())
            cs5 = ax5.contourf(Lon,lat,Qdiff,30,transform=ccrs.PlateCarree(),cmap='PRGn',vmin=-10,vmax=10)
            ax5.coastlines()
            if i == 0:
                plt.title('Difference w.r.t. $S_0$',fontsize=14)
            if i == 0:
                cticks=np.around(np.linspace(-10,10,11),decimals=1)
                cbar_ax = fig.add_axes([0.82,-0.05,0.15,0.02])
                cb5 = fig.colorbar(mappable=None, norm=Normalize(vmin=-10,vmax=10), cmap='PRGn',spacing='proportional',orientation='horizontal',cax=cbar_ax,ticks=cticks)
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
    plt.suptitle('Estimated Inversion Strength Response to $CO_2$ Doubling',fontsize=24,y=1.02)
    #plt.suptitle('Estimated Inversion Strength Response to Equivalent Radiative Forcing',fontsize=24,y=1.02)
    
    fig.tight_layout(pad=0.2)
    plt.show()
    fig.savefig(figure_path+'Map_eis_2x_response.pdf',bbox_inches='tight')
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
    