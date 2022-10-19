#!/usr/bin/env python3

#
# Script for plotting figure 3.2. for altercation, make necessary comments and uncomments
# of particular axes settings in the plotting section.
#

import os
import sys
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

figure_path = '/home/brandonsmith/climate-gcm-bps/plots/' # path where figures are to be saved
casenames = ['1.0','0.9','0.95','1.05']

for CASENAME in casenames:
    run = 'Control'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME #location of parent data set
    infile = '/Merged_'+run+'_'+CASENAME+'.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME #location of child data set
    outfile_control = '/ZM_Background_mappable_'+run+'_'+CASENAME+'.nc'
    T700 = '/map_lts_700Temp_'+run+'_'+CASENAME+'.nc' #####
    T1000 = '/map_lts_1000Temp_'+run+'_'+CASENAME+'.nc'   #
    outfilesub = '/map_lts_sub_'+run+'_'+CASENAME+'.nc'   #
    outfile700 = '/map_lts_700_'+run+'_'+CASENAME+'.nc'   #
    outfile1000 = '/map_lts_1000_'+run+'_'+CASENAME+'.nc' #
    q850 = '/map_q850_'+run+'_'+CASENAME+'.nc'            ### THESE FILES ARE FOR THE CALCULATION OF ESTIMATED INVERSION STRENGTH.
    qs850 = '/map_qs850_'+run+'_'+CASENAME+'.nc'          #
    tempsum = '/map_tempsum_'+run+'_'+CASENAME+'.nc'      #
    z700 = '/map_z700_'+run+'_'+CASENAME+'.nc'            #
    LCL = '/map_lcl_'+run+'_'+CASENAME+'.nc'              #
    outfile_LCL = '/map_LCL_'+run+'_'+CASENAME+'.nc'      #
    outfile_eis = '/map_eis_'+run+'_'+CASENAME+'.nc' ######
    if os.path.isfile(inpath+infile): # check that input file exists
        if not os.path.isfile(outpath+outfile_control): # create plotting file if it does not exist. Else, continue to plotting.
            syscall = '/usr/bin/cdo timmean -seltimestep,-360/-1 -select,name=TS,CLDLOW,CLDHGH,TGCLDLWP,TGCLDIWP,LHFLX '+inpath+infile+' '+outpath+outfile_control
            print(syscall)
        if not os.path.isfile(outpath+T700):
            syscall = '/usr/bin/cdo seltimestep,-360/-1 -sellevel,691.389430314302 -select,name=T '+inpath+infile+ ' '+outpath+T700
            print(syscall)
        if not os.path.isfile(outpath+outfile700):
            syscall = '/usr/bin/cdo timmean -select,name=lts -expr,\'lts=(T*(1000/691.389430314302)^0.286)\' '+outpath+T700+' '+outpath+outfile700
            print(syscall)
        if not os.path.isfile(outpath+T1000):
            syscall = '/usr/bin/cdo seltimestep,-360/-1 -select,name=TS,PS '+inpath+infile+ ' '+outpath+T1000
            print(syscall)
        if not os.path.isfile(outpath+outfile1000):
            syscall = '/usr/bin/cdo timmean -select,name=lts -expr,\'lts=(TS*(1000/(PS*0.01))^0.286)\' '+inpath+infile+' '+outpath+outfile1000
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
# ONCE LINES ARE PRINTED, COPY OVER TO A SHELL SCRIPT ON A SERVER WHERE CLIMATE DATA OPERATORS IS INSTALLED AND RUN TO CREATE THE FILES.
# ALTERNATIVELY, REPLACE print() STATEMENTS WITH os.system() FUNCTION CALL TO CALL DIRECTLY TO COMMAND LINE.
    
# Load variables and perform calculations from outfiles
i = 0
fig = plt.figure(figsize=(20,20))
outfile_1C = '/home/brandonsmith/modeloutput/Control/1.0/ZM_Background_mappable_Control_1.0.nc' # file for other variables
outfile_eis1 = '/home/brandonsmith/modeloutput/Control/1.0/map_eis_Control_1.0.nc' # Estimated Inversion Strength (EIS) file
for CASENAME in casenames:
    outfile_control = '/ZM_Background_mappable_Control_'+CASENAME+'.nc'
    outfilesub = '/map_eis_Control_'+CASENAME+'.nc'
    outpath_control = '/home/brandonsmith/modeloutput/Control/'+CASENAME
    dsloc = outpath_control+outfilesub
    dsloc_control = outpath_control+outfile_control
    if os.path.isfile(dsloc_control) and os.path.isfile(dsloc):
        # Load in variables
        dsc = netCDF4.Dataset(dsloc_control)
        T= dsc.variables['TS'][:]
        cldlow = dsc.variables['CLDLOW'][:]
        cldhgh = dsc.variables['CLDHGH'][:]
        LWP = dsc.variables['TGCLDLWP'][:]
        IWP = dsc.variables['TGCLDIWP'][:]
        LHFLX = dsc.variables['LHFLX'][:]
        lat = dsc.variables['lat'][:]
        lon = dsc.variables['lon'][:]
        dsc.close()
        
        # Load in variables needed for calculating EIS
        ds = netCDF4.Dataset(dsloc)
        lts = ds.variables['lts'][:]
        tempsum = ds.variables['TREFHT'][:]
        qs850 = ds.variables['smr'][:]
        z700 = ds.variables['Z3'][:]
        LCL = ds.variables['LCL'][:]
        ds.close()

        C1 = netCDF4.Dataset(outfile_1C)
        T1= C1.variables['TS'][:]
        cldlow1 = C1.variables['CLDLOW'][:]
        cldhgh1 = C1.variables['CLDHGH'][:]
        LWP1 = C1.variables['TGCLDLWP'][:]
        IWP1 = C1.variables['TGCLDIWP'][:]
        LHFLX1 = C1.variables['LHFLX'][:]
        C1.close()

        C2 = netCDF4.Dataset(outfile_eis1)
        lts1 = C2.variables['lts'][:]
        tempsum1 = C2.variables['TREFHT'][:]
        qs8501 = C2.variables['smr'][:]
        z7001 = C2.variables['Z3'][:]
        LCL1 = C2.variables['LCL'][:]
        C2.close()
        
        # calculate EIS (defined as "eis" and "eis1")
        pal = (0.00989)*(1-(1+2450000*qs850/(287.058*(tempsum/2)))/(1+(2450000**2)*qs850/(993*461.4*((tempsum/2)**2))))
        pal1 = (0.00989)*(1-(1+2450000*qs8501/(287.058*(tempsum1/2)))/(1+(2450000**2)*qs8501/(993*461.4*((tempsum1/2)**2))))

        T = T.squeeze()
        T1 = T1.squeeze()
        LWP = LWP*1000
        IWP = IWP*1000
        LWP1 = LWP1*1000
        IWP1 = IWP1*1000

        eis = lts - pal*(z700-LCL)
        eis1 = lts1 - pal1*(z7001-LCL1)
        
        # these lines get rid of the irritating white space line through the center of each plot. artifact of cartopy.
        Lon = lon
        T = np.squeeze(T)
        T, Lon = add_cyclic_point(T, coord=lon)
        Lon = lon
        T1 = np.squeeze(T1)
        T1, Lon = add_cyclic_point(T1, coord=lon)
        Lon = lon
        cldlow = np.squeeze(cldlow)
        cldlow, Lon = add_cyclic_point(cldlow, coord=lon)
        Lon = lon
        cldlow1 = np.squeeze(cldlow1)
        cldlow1, Lon = add_cyclic_point(cldlow1, coord=lon)
        Lon = lon
        cldhgh = np.squeeze(cldhgh)
        cldhgh, Lon = add_cyclic_point(cldhgh, coord=lon)
        Lon = lon
        cldhgh1 = np.squeeze(cldhgh1)
        cldhgh1, Lon = add_cyclic_point(cldhgh1, coord=lon)
        Lon = lon
        LWP = np.squeeze(LWP)
        LWP, Lon = add_cyclic_point(LWP, coord=lon)
        Lon = lon
        LWP1 = np.squeeze(LWP1)
        LWP1, Lon = add_cyclic_point(LWP1, coord=lon)
        Lon = lon
        IWP = np.squeeze(IWP)
        IWP, Lon = add_cyclic_point(IWP, coord=lon)
        Lon = lon
        IWP1 = np.squeeze(IWP1)
        IWP1, Lon = add_cyclic_point(IWP1, coord=lon)
        Lon = lon
        LHFLX = np.squeeze(LHFLX)
        LHFLX, Lon = add_cyclic_point(LHFLX, coord=lon)
        Lon = lon
        LHFLX1 = np.squeeze(LHFLX1)
        LHFLX1, Lon = add_cyclic_point(LHFLX1, coord=lon)
        Lon = lon
        eis = np.squeeze(eis)
        eis, Lon = add_cyclic_point(eis, coord=lon)
        Lon = lon
        eis1 = np.squeeze(eis1)
        eis1, Lon = add_cyclic_point(eis1, coord=lon)
        
        # Calculate differences between each background climate and the control run.
        deltaT = T.squeeze() - T1.squeeze()
        delta_cldlow = cldlow.squeeze() - cldlow1.squeeze()
        delta_cldhgh = cldhgh.squeeze() - cldhgh1.squeeze()
        deltaLWP = LWP.squeeze() - LWP1.squeeze()
        deltaIWP = IWP.squeeze() - IWP1.squeeze()
        delta_LHFLX = LHFLX.squeeze() - LHFLX1.squeeze()
        delta_eis = eis.squeeze() - eis1.squeeze()
        
        #plot variable
        #casenames = ['1.0','0.9','0.95','1.05']
        gs = gridspec.GridSpec(21, 13)
        if i >= 1:
            ax = plt.subplot(gs[:3, 3*i+1:3*i+4],projection=ccrs.Robinson())
        else:
            ax = plt.subplot(gs[:3, 3*i:3*i+3],projection=ccrs.Robinson())
        if i == 0:
            cs = ax.contourf(Lon,lat,T.squeeze(),13,transform=ccrs.PlateCarree(),cmap='gist_heat',vmin=180,vmax=300)
            ax.set_title('Solar Multiplier: $S_0$',fontsize=16)
            ax.coastlines()
            Cticks=np.around(np.linspace(180,300,13),decimals=1)
            Cbar_ax = inset_axes(ax,width="5%",  height="100%",loc='center right',borderpad=-5)
            cb = fig.colorbar(cs, spacing='proportional',orientation='vertical',cax=Cbar_ax,ticks=Cticks)
            #cb.formatter.set_powerlimits((0,0))
            cb.set_ticklabels(Cticks.astype(int).astype(str),fontsize=12)
            #cb.set_ticklabels(cticks.astype(str),fontsize=10)
            cb.set_label('Surface Temperature (K)',fontsize=10)
            Cbar_ax.yaxis.set_label_position('left')
        else:
            cs = ax.contourf(Lon,lat,deltaT,61,transform=ccrs.PlateCarree(),cmap='RdBu_r',vmin=-5,vmax=5)
            c = ax.contour(Lon,lat,T.squeeze(),13,vmin=180,vmax=300,colors='black',linewidths=0.5)
            ax.set_title('Solar Multiplier: '+CASENAME+'*$S_0$',fontsize=16)
            ax.coastlines()
            if i == 3:
                cticks=np.around(np.linspace(-5,5,11),decimals=2)
                cbar_ax = inset_axes(ax,width="5%",  height="100%",loc='center right',borderpad=-5)
                ax.set_ylabel('Hybrid Sigma-Pressure level (mb)',fontsize=12)
                #cb = fig.colorbar(cs, spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=cticks)
                cb = fig.colorbar(mappable=None, norm=Normalize(vmin=-5,vmax=5), cmap='RdBu_r',spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=cticks)
            #cb.formatter.set_powerlimits((0,0))
            #cb.set_ticklabels(cticks.astype(str),fontsize=7)
                cb.set_ticklabels(cticks.astype(int).astype(str),fontsize=12)
                cbar_ax.yaxis.set_ticks_position('left')
                cb.set_label('Change in Temperature (K)',fontsize=10)
        if i >= 1:
            ax2 = plt.subplot(gs[3:6, 3*i+1:3*i+4],projection=ccrs.Robinson())
        else:
            ax2 = plt.subplot(gs[3:6, 3*i:3*i+3],projection=ccrs.Robinson())
        if i == 0:
            cs2 = ax2.contourf(Lon,lat,cldlow.squeeze(),11,transform=ccrs.PlateCarree(),cmap='gist_gray',vmin=0,vmax=1)
            ax2.coastlines()
            Cticks=np.around(np.linspace(0,1,6),decimals=1)
            Cbar_ax2 = inset_axes(ax2,width="5%",  height="100%",loc='center right',borderpad=-5)
            cb2 = fig.colorbar(cs2, spacing='proportional',orientation='vertical',cax=Cbar_ax2,ticks=Cticks)
            #cb.formatter.set_powerlimits((0,0))
            #cb2.set_ticklabels(Cticks.astype(int).astype(str),fontsize=7)
            cb2.set_ticklabels(Cticks.astype(str),fontsize=12)
            cb2.set_label('Low Cloud Fraction',fontsize=10)
            Cbar_ax2.yaxis.set_label_position('left')
        else:
            cs2 = ax2.contourf(Lon,lat,delta_cldlow,61,transform=ccrs.PlateCarree(),cmap='BrBG',vmin=-0.2,vmax=0.2)
            c2 = ax2.contour(Lon,lat,cldlow.squeeze(),13,vmin=0,vmax=1,colors='black',linewidths=0.5)
            ax2.coastlines()
            if i == 3:
                cticks=np.around(np.linspace(-0.2,0.2,5),decimals=2)
                cbar_ax = inset_axes(ax2,width="5%",  height="100%",loc='center right',borderpad=-5)
                ax2.set_ylabel('Hybrid Sigma-Pressure level (mb)',fontsize=12)
                #cb = fig.colorbar(cs, spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=cticks)
                cb2 = fig.colorbar(mappable=None, norm=Normalize(vmin=-0.2,vmax=0.2), cmap='BrBG',spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=cticks)
            #cb.formatter.set_powerlimits((0,0))
                cb2.set_ticklabels(cticks.astype(str),fontsize=12)
                cbar_ax.yaxis.set_ticks_position('left')
                #cb2.set_ticklabels(cticks.astype(int).astype(str),fontsize=10)
                cb2.set_label('Cloud Fraction Difference',fontsize=10)

        if i >= 1:
            ax3 = plt.subplot(gs[6:9, 3*i+1:3*i+4],projection=ccrs.Robinson())
        else:
            ax3 = plt.subplot(gs[6:9, 3*i:3*i+3],projection=ccrs.Robinson())
        if i == 0:
            cs3 = ax3.contourf(Lon,lat,cldhgh.squeeze(),12,transform=ccrs.PlateCarree(),cmap='gist_gray',vmin=0,vmax=1)
            ax3.coastlines()
            Cticks=np.around(np.linspace(0,1,6),decimals=1)
            Cbar_ax3 = inset_axes(ax3,width="5%",  height="100%",loc='center right',borderpad=-5)
            cb3 = fig.colorbar(cs3, spacing='proportional',orientation='vertical',cax=Cbar_ax3,ticks=Cticks)
            #cb2.set_ticklabels(Cticks.astype(int).astype(str),fontsize=7)
            cb3.set_ticklabels(Cticks.astype(str),fontsize=12)
            cb3.set_label('High Cloud Fraction',fontsize=10)
            Cbar_ax3.yaxis.set_label_position('left')
        else:
            cs3 = ax3.contourf(Lon,lat,delta_cldhgh,61,transform=ccrs.PlateCarree(),cmap='BrBG',vmin=-0.2,vmax=0.2)
            c3 = ax3.contour(Lon,lat,cldhgh.squeeze(),12,vmin=0,vmax=1,colors='black',linewidths=0.5)
            ax3.coastlines()
            if i == 3:
                cticks=np.around(np.linspace(-0.2,0.2,5),decimals=1)
                cbar_ax = inset_axes(ax3,width="5%",  height="100%",loc='center right',borderpad=-5)
                ax3.set_ylabel('Hybrid Sigma-Pressure level (mb)',fontsize=12)
                #cb = fig.colorbar(cs, spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=cticks)
                cb3 = fig.colorbar(mappable=None, norm=Normalize(vmin=-0.2,vmax=0.2), cmap='BrBG',spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=cticks)
                cb3.set_ticklabels(cticks.astype(str),fontsize=12)
                cbar_ax.yaxis.set_ticks_position('left')
                #cb2.set_ticklabels(cticks.astype(int).astype(str),fontsize=10)
                cb3.set_label('Cloud Fraction Difference',fontsize=10)
        if i >= 1:
            ax4 = plt.subplot(gs[9:12, 3*i+1:3*i+4],projection=ccrs.Robinson())
        else:
            ax4 = plt.subplot(gs[9:12, 3*i:3*i+3],projection=ccrs.Robinson())
        if i == 0:
            cs4 = ax4.contourf(Lon,lat,LWP.squeeze(),11,transform=ccrs.PlateCarree(),cmap='Greens',vmin=0,vmax=200)
            ax4.coastlines()
            ax4.set_xlabel('Latitude (Deg N)',fontsize=14)
            Cticks=np.around(np.linspace(0,200,11),decimals=1)
            Cbar_ax4 = inset_axes(ax4,width="5%",  height="100%",loc='center right',borderpad=-5)
            cb4 = fig.colorbar(cs4, spacing='proportional',orientation='vertical',cax=Cbar_ax4,ticks=Cticks)
            #cb.formatter.set_powerlimits((0,0))
            cb4.set_ticklabels(Cticks.astype(int).astype(str),fontsize=12)
            #cb.set_ticklabels(cticks.astype(str),fontsize=10)
            cb4.set_label('LWP ($gm^2$)',fontsize=10)
            Cbar_ax4.yaxis.set_label_position('left')
        else:
            cs4 = ax4.contourf(Lon,lat,deltaLWP,61,transform=ccrs.PlateCarree(),cmap='BrBG',vmin=-20,vmax=20)
            c4 = ax4.contour(Lon,lat,LWP.squeeze(),12,vmin=0,vmax=100,colors='black',linewidths=0.5)
            ax4.coastlines()
            ax4.set_xlabel('Latitude (Deg N)',fontsize=14)
            if i == 3:
                cticks=np.around(np.linspace(-20,20,11),decimals=2)
                cbar_ax = inset_axes(ax4,width="5%",  height="100%",loc='center right',borderpad=-5)
                ax4.set_ylabel('Hybrid Sigma-Pressure level (mb)',fontsize=12)
                #cb = fig.colorbar(cs, spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=cticks)
                cb4 = fig.colorbar(mappable=None, norm=Normalize(vmin=-20,vmax=20), cmap='BrBG',spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=cticks)
            #cb.formatter.set_powerlimits((0,0))
            #cb.set_ticklabels(cticks.astype(str),fontsize=7)
                cb4.set_ticklabels(cticks.astype(int).astype(str),fontsize=12)
                cbar_ax.yaxis.set_ticks_position('left')
                cb4.set_label('LWP Difference($gm^2$)',fontsize=10)
        if i >= 1:
            ax5 = plt.subplot(gs[12:15, 3*i+1:3*i+4],projection=ccrs.Robinson())
        else:
            ax5 = plt.subplot(gs[12:15, 3*i:3*i+3],projection=ccrs.Robinson())
        if i == 0:
            cs5 = ax5.contourf(Lon,lat,IWP.squeeze(),12,transform=ccrs.PlateCarree(),cmap='Blues',vmin=0,vmax=100)
            ax5.coastlines()
            ax5.set_xlabel('Latitude (Deg N)',fontsize=14)
            Cticks=np.around(np.linspace(0,100,6),decimals=1)
            Cbar_ax5 = inset_axes(ax5,width="5%",  height="100%",loc='center right',borderpad=-5)
            cb5 = fig.colorbar(cs5, spacing='proportional',orientation='vertical',cax=Cbar_ax5,ticks=Cticks)
            #cb.formatter.set_powerlimits((0,0))
            cb5.set_ticklabels(Cticks.astype(int).astype(str),fontsize=12)
            #cb.set_ticklabels(cticks.astype(str),fontsize=10)
            cb5.set_label('IWP ($gm^2$)',fontsize=10)
            Cbar_ax5.yaxis.set_label_position('left')
        else:
            cs5 = ax5.contourf(Lon,lat,deltaIWP,61,transform=ccrs.PlateCarree(),cmap='BrBG',vmin=-20,vmax=20)
            c5 = ax5.contour(Lon,lat,IWP.squeeze(),12,vmin=0,vmax=240,colors='black',linewidths=0.5)
            ax5.coastlines()
            ax5.set_xlabel('Latitude (Deg N)',fontsize=14)
            if i == 3:
                cticks=np.around(np.linspace(-20,20,11),decimals=2)
                cbar_ax = inset_axes(ax5,width="5%",  height="100%",loc='center right',borderpad=-5)
                ax5.set_ylabel('Hybrid Sigma-Pressure level (mb)',fontsize=12)
                #cb = fig.colorbar(cs, spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=cticks)
                cb5 = fig.colorbar(mappable=None, norm=Normalize(vmin=-20,vmax=20), cmap='BrBG',spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=cticks)
            #cb.formatter.set_powerlimits((0,0))
            #cb.set_ticklabels(cticks.astype(str),fontsize=7)
                cb5.set_ticklabels(cticks.astype(int).astype(str),fontsize=12)
                cbar_ax.yaxis.set_ticks_position('left')
                cb5.set_label('IWP Difference ($gm^2$)',fontsize=10)
        if i >= 1:
            ax6 = plt.subplot(gs[15:18, 3*i+1:3*i+4],projection=ccrs.Robinson())
        else:
            ax6 = plt.subplot(gs[15:18, 3*i:3*i+3],projection=ccrs.Robinson())
        if i == 0:
            cs6 = ax6.contourf(Lon,lat,LHFLX.squeeze(),61,transform=ccrs.PlateCarree(),cmap='gist_heat',vmin=0,vmax=240)
            ax6.coastlines()
            ax6.set_xlabel('Latitude (Deg N)',fontsize=14)
            Cticks=np.around(np.linspace(0,240,6),decimals=1)
            #Cbar_ax6 = fig.add_axes([0.31,0.16,0.01,0.10])
            Cbar_ax6 = inset_axes(ax6,width="5%",  height="100%",loc='center right',borderpad=-5)
            cb6 = fig.colorbar(cs6, spacing='proportional',orientation='vertical',cax=Cbar_ax6,ticks=Cticks)
            #cb.formatter.set_powerlimits((0,0))
            cb6.set_ticklabels(Cticks.astype(int).astype(str),fontsize=12)
            #cb6.set_ticklabels(Cticks.astype(str),fontsize=12)
            cb6.set_label('Latent Heat Flux ($W/m^2$)',fontsize=10)
            Cbar_ax6.yaxis.set_label_position('left')
        else:
            cs6= ax6.contourf(Lon,lat,delta_LHFLX,61,transform=ccrs.PlateCarree(),cmap='PuOr_r',vmin=-50,vmax=50)
            c6 = ax6.contour(Lon,lat,LHFLX.squeeze(),12,vmin=0,vmax=240,colors='black',linewidths=0.5)
            ax6.coastlines()
            ax6.set_xlabel('Latitude (Deg N)',fontsize=14)
            if i == 3:
                cticks=np.around(np.linspace(-50,50,11),decimals=2)
                cbar_ax = inset_axes(ax6,width="5%",  height="100%",loc='center right',borderpad=-5)
                ax6.set_ylabel('Hybrid Sigma-Pressure level (mb)',fontsize=12)
                #cb = fig.colorbar(cs, spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=cticks)
                cb6 = fig.colorbar(mappable=None, norm=Normalize(vmin=-25,vmax=25), cmap='PuOr_r',spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=cticks)
            #cb.formatter.set_powerlimits((0,0))
                #cb6.set_ticklabels(cticks.astype(str),fontsize=12)
                cb6.set_ticklabels(cticks.astype(int).astype(str),fontsize=12)
                cbar_ax.yaxis.set_ticks_position('left')
                cb6.set_label('LHF difference ($W/m^2$)',fontsize=10)
        if i >= 1:
            ax7 = plt.subplot(gs[18:, 3*i+1:3*i+4],projection=ccrs.Robinson())
        else:
            ax7 = plt.subplot(gs[18:, 3*i:3*i+3],projection=ccrs.Robinson())
        if i == 0:
            cs7 = ax7.contourf(Lon,lat,eis.squeeze(),36,transform=ccrs.PlateCarree(),cmap='PRGn',vmin=-40,vmax=40)
            ax7.coastlines()
            Cticks=np.around(np.linspace(-40,40,9),decimals=1)
            Cbar_ax7 = inset_axes(ax7,width="5%",  height="100%",loc='center right',borderpad=-5)
            cb7 = fig.colorbar(cs7, spacing='proportional',orientation='vertical',cax=Cbar_ax7,ticks=Cticks)
            #cb.formatter.set_powerlimits((0,0))
            cb7.set_ticklabels(Cticks.astype(int).astype(str),fontsize=12)
            #cb.set_ticklabels(cticks.astype(str),fontsize=10)
            cb7.set_label('EIS (K)',fontsize=10)
            Cbar_ax7.yaxis.set_label_position('left')            
        else:
            cs7 = ax7.contourf(Lon,lat,delta_eis,61,transform=ccrs.PlateCarree(),cmap='PRGn',vmin=-5,vmax=5)
            c7 = ax7.contour(Lon,lat,eis.squeeze(),13,vmin=0,vmax=40,colors='black',linewidths=0.5)
            ax7.coastlines()
            if i == 3:
                cticks=np.around(np.linspace(-5,5,11),decimals=2)
                cbar_ax = inset_axes(ax7,width="5%",  height="100%",loc='center right',borderpad=-5)
                ax7.set_ylabel('Hybrid Sigma-Pressure level (mb)',fontsize=12)
                #cb = fig.colorbar(cs, spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=cticks)
                cb7 = fig.colorbar(mappable=None, norm=Normalize(vmin=-5,vmax=5), cmap='PRGn',spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=cticks)
            #cb.formatter.set_powerlimits((0,0))
            #cb.set_ticklabels(cticks.astype(str),fontsize=7)
                cb7.set_ticklabels(cticks.astype(int).astype(str),fontsize=12)
                cbar_ax.yaxis.set_ticks_position('left')
                cb7.set_label('EIS difference (K)',fontsize=10)
        i = i+1
    else:
        print('No such file or directory: '+dsloc_control+' or '+dsloc)

#    cb.set_label('Temperature (K)')
plt.suptitle('Comparison of Geospatially-Resolved Fields Between Background Climate States',fontsize=24,y=0.93)
#fig.text(-0.04, 0.5, 'Sigma Pressure Level (mb)', va='center', rotation='vertical')
#fig.tight_layout()
plt.show()
fig.savefig(figure_path+'Reference_Climate_Mappable_Variables.pdf',bbox_inches='tight')

    