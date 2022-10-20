#!/usr/bin/env python3

#
# Script for plotting figure 3.1. for altercation, make necessary comments and uncomments
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

figure_path = '/home/brandonsmith/climate-gcm-bps/plots/' # Change to the path where you want your figures saved
casenames = ['1.0','0.9','0.95','1.05']

# Create the files that will be read for plotting purposes
for CASENAME in casenames:
    run = 'Control'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME #location of parent data set
    infile = '/Merged_'+run+'_'+CASENAME+'.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME #location of child data set
    outfile_vc = '/ZMV_'+run+'_'+CASENAME+'.nc' #file for computing mass stream function from V profiles
    outfile_control_stf = '/MASTRFU_'+run+'_'+CASENAME+'.nc' #file for mass stream function
    outfile_control_vars = '/ZM_background_vars_'+run+'_'+CASENAME+'.nc' #remaining variables
    outfile_control = '/ZM_Background_'+run+'_'+CASENAME+'.nc' #file for plotting
    outfile_winds = '/winds_and_variances_'+run+'_'+CASENAME+'.nc' #file for calculating Eddy Kinetic Energy
    if os.path.isfile(inpath+infile): #check that the parent file exists
        if not os.path.isfile(outpath+outfile_vc): #perform file creation if the plotting file does not exist. Else, continue to plotting.
            syscall = '/usr/bin/cdo timmean -zonmean -seltimestep,-360/-1 -select,name='+field+' '+inpath+infile+' '+outpath+outfile_vc
            print(syscall)
        if not os.path.isfile(outpath+outfile_control_stf):
            syscall = '/usr/bin/cdo mastrfu '+outpath+outfile_vc+' '+outpath+outfile_control_stf
            print(syscall)
        if not os.path.isfile(outpath+outfile_control_vars):
            syscall = '/usr/bin/cdo timmean -zonmean -seltimestep,-360/-1 -select,name=RELHUM,T,OMEGA,Q '+inpath+infile +' '+ outpath+outfile_control_vars
            print(syscall)
        if not os.path.isfile(outpath+outfile_control):
            syscall = '/usr/bin/cdo merge '+outpath+outfile_control_stf+' '+outpath+outfile_control_vars+' '+outpath+outfile_control
            print(syscall)
        if not os.path.isfile(outpath+outfile_winds):
            syscall = '/usr/bin/cdo -seltimestep,-360/-1 -select,name=V,VV,U,UU '+inpath+infile+' '+outpath+outfile_winds
            print(syscall)
# ONCE LINES ARE PRINTED, COPY OVER TO A SHELL SCRIPT ON A SERVER WHERE CLIMATE DATA OPERATORS IS INSTALLED AND RUN TO CREATE THE FILES.
# ALTERNATIVELY, REPLACE print() STATEMENTS WITH os.system() FUNCTION CALL TO CALL DIRECTLY TO COMMAND LINE.
    

# Load variables and perform calculations from outfiles
i = 0
fig = plt.figure(figsize=(30,20))
outfile_1C = '/home/brandonsmith/modeloutput/Control/1.0/ZM_Background_Control_1.0.nc' #file for calculating most variables
outfile_1C_2 = '/home/brandonsmith/modeloutput/Control/1.0/winds_and_variances_Control_1.0.nc' #file for calculating Eddy Kinetic Energy Profiles
for CASENAME in casenames:
    outfile_control = '/ZM_Background_Control_'+CASENAME+'.nc'
    outfile_winds = '/winds_and_variances_Control_'+CASENAME+'.nc'
    outpath_control = '/home/brandonsmith/modeloutput/Control/'+CASENAME
    dsloc = outpath_control+outfile_winds
    dsloc_control = outpath_control+outfile_control
    if os.path.isfile(dsloc_control) and os.path.isfile(dsloc):
        # Load in variables for plotting and calculations
        print('loading variables')
        dsc = netCDF4.Dataset(dsloc_control)
        T= dsc.variables['T'][:]
        mastrfu = dsc.variables['OMEGA'][:]
        CLOUD = dsc.variables['RELHUM'][:]
        Q = dsc.variables['Q'][:]
        lat = dsc.variables['lat'][:]
        lev = dsc.variables['lev'][:]
        dsc.close()

        ds = netCDF4.Dataset(dsloc)
        U = ds.variables['U'][:]
        V= ds.variables['V'][:]
        UU = ds.variables['UU'][:]
        VV = ds.variables['VV'][:]
        ds.close()

        C1 = netCDF4.Dataset(outfile_1C)
        T1 = C1.variables['T'][:]
        CLOUD1 = C1.variables['RELHUM'][:]
        mastrfu1 = C1.variables['OMEGA'][:]
        Q1 = C1.variables['Q'][:]
        C1.close()

        C2 = netCDF4.Dataset(outfile_1C_2)
        U1 = C2.variables['U'][:]
        V1= C2.variables['V'][:]
        UU1 = C2.variables['UU'][:]
        VV1 = C2.variables['VV'][:]
        C2.close()
        print('computing eddy components')
        # Calculate Eddy Components of wind
        eddy_VV = VV - V**2
        eddy_UU = UU - U**2
        eddy_VV1 = VV1 - V1**2
        eddy_UU1 = UU1 - U1**2
        print('Computing EKE')
        # Calculate Eddy Kinetic Energy from Eddy Components
        EKE = 0.5*(eddy_VV+eddy_UU)
        EKE1 = 0.5*(eddy_VV1+eddy_UU1)
        # Compute zonal and temporal averages
        EKE = np.average(EKE,0)
        EKE = np.average(EKE,2)
        EKE1 = np.average(EKE1,0)
        EKE1 = np.average(EKE1,2)
        T = T.squeeze()
        T1 = T1.squeeze()
        CLOUD = CLOUD.squeeze()
        CLOUD1 = CLOUD1.squeeze()
        mastrfu = mastrfu.squeeze()
        mastrfu1 = mastrfu1.squeeze()
        Q = Q.squeeze() * 1000
        Q1 = Q1.squeeze() * 1000
        # Calculate difference between each background climate state and the control run
        deltaT = T - T1
        deltaCloud = CLOUD - CLOUD1
        deltaChi = mastrfu - mastrfu1
        EKE_diff = EKE.squeeze() - EKE1.squeeze()
        Q_diff = Q - Q1


        #plot variable
        gs = gridspec.GridSpec(15, 13)
        if i >= 1:
            ax = plt.subplot(gs[:3, 3*i+1:3*i+4])
        else:
            ax = plt.subplot(gs[:3, 3*i:3*i+3])
        if i == 0:
            cs = ax.contourf(lat,lev,T,13,cmap='gist_heat',vmin=180,vmax=300)
            ax.set_title('Solar Constant: $S_0$',fontsize=18)
            plt.yticks(np.linspace(1000,0,5),fontsize=14)
            plt.gca().invert_yaxis()
            plt.xticks(np.linspace(-90,90,7),fontsize=14)
            Cticks=np.around(np.linspace(180,300,7),decimals=1)
            cbar_ax = inset_axes(ax,width="5%",  height="100%",loc='center right',borderpad=-5)
            cb = fig.colorbar(cs, spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=Cticks)
            #cb.formatter.set_powerlimits((0,0))
            cb.set_ticklabels(Cticks.astype(int).astype(str),fontsize=12)
            #cb.set_ticklabels(cticks.astype(str),fontsize=10)
            cb.set_label('Temperature (K)',fontsize=12)
            cbar_ax.yaxis.set_ticks_position('left')
            #ax.set_ylabel('Hybrid Sigma-Pressure level (mb)',fontsize=10)
            ax.text(-0.3, 0.5, 'a', va='center', rotation='horizontal',transform=ax.transAxes,fontsize=16,color='black')
        else:
            cs = ax.contourf(lat,lev,deltaT,61,cmap='RdBu_r',vmin=-6,vmax=6)
            c = ax.contour(lat,lev,T,13,vmin=180,vmax=300,colors='black',linewidths=0.5)
            ax.set_title('Solar Constant: '+CASENAME+'*$S_0$',fontsize=18)
            plt.yticks(np.linspace(1000,0,5),fontsize=14)
            plt.gca().invert_yaxis()
            plt.xticks(np.linspace(-90,90,7),fontsize=14)
            if i == 3:
                cticks=np.around(np.linspace(-6,6,7),decimals=2)
                cbar_ax = inset_axes(ax,width="5%",  height="100%",loc='center right',borderpad=-5)
                #cb = fig.colorbar(cs, spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=cticks)
                cb = fig.colorbar(mappable=None, norm=Normalize(vmin=-6,vmax=6), cmap='RdBu_r',spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=cticks)
            #cb.formatter.set_powerlimits((0,0))
            #cb.set_ticklabels(cticks.astype(str),fontsize=7)
                cb.set_ticklabels(cticks.astype(int).astype(str),fontsize=12)
                cbar_ax.yaxis.set_ticks_position('left')
                cb.set_label('Difference (K)',fontsize=12)
        if i >= 1:
            ax2 = plt.subplot(gs[3:6, 3*i+1:3*i+4])
        else:
            ax2 = plt.subplot(gs[3:6, 3*i:3*i+3])
        if i == 0:
            cs2 = ax2.contourf(lat,lev,CLOUD,11,cmap='Greens',vmin=0,vmax=100)
            plt.yticks(np.linspace(1000,0,5),fontsize=14)
            plt.gca().invert_yaxis()
            plt.xticks(np.linspace(-90,90,7),fontsize=14)
            Cticks=np.around(np.linspace(0,100,6),decimals=1)
            cbar_ax = inset_axes(ax2,width="5%",  height="100%",loc='center right',borderpad=-5)
            cb2 = fig.colorbar(cs2, spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=Cticks)
            #cb.formatter.set_powerlimits((0,0))
            cb2.set_ticklabels(Cticks.astype(int).astype(str),fontsize=12)
            #cb2.set_ticklabels(Cticks.astype(str),fontsize=12)
            cb2.set_label('RH %',fontsize=12)
            #cbar_ax.yaxis.set_label_position('left')
            cbar_ax.yaxis.set_ticks_position('left')
            #ax2.set_ylabel('Hybrid Sigma-Pressure level (mb)',fontsize=10)
            ax2.text(-0.3, 0.5, 'b', va='center', rotation='horizontal',transform=ax2.transAxes,fontsize=16,color='black')
        else:
            cs2 = ax2.contourf(lat,lev,deltaCloud,61,cmap='BrBG',vmin=-10,vmax=10)
            c2 = ax2.contour(lat,lev,CLOUD,13,vmin=0,vmax=1,colors='black',linewidths=0.5)
            plt.yticks(np.linspace(1000,0,5),fontsize=14)
            plt.gca().invert_yaxis()
            plt.xticks(np.linspace(-90,90,7),fontsize=14)
            if i == 3:
                cticks=np.around(np.linspace(-10,10,5),decimals=2)
                cbar_ax = inset_axes(ax2,width="5%",  height="100%",loc='center right',borderpad=-5)
                #cb = fig.colorbar(cs, spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=cticks)
                cb2 = fig.colorbar(mappable=None, norm=Normalize(vmin=-10,vmax=10), cmap='BrBG',spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=cticks)
            #cb.formatter.set_powerlimits((0,0))
                #cb2.set_ticklabels(cticks.astype(str),fontsize=12)
                cb2.set_ticklabels(cticks.astype(int).astype(str),fontsize=12)
                cbar_ax.yaxis.set_ticks_position('left')
                cb2.set_label('Difference',fontsize=12)

        if i >= 1:
            ax3 = plt.subplot(gs[6:9, 3*i+1:3*i+4])
        else:
            ax3 = plt.subplot(gs[6:9, 3*i:3*i+3])
        if i == 0:
            cs3 = ax3.contourf(lat,lev,Q,181,cmap='Greens',vmin=0,vmax=25)
            plt.yticks(np.linspace(1000,0,5),fontsize=14)
            plt.gca().invert_yaxis()
            plt.xticks(np.linspace(-90,90,7),fontsize=14)
            Cticks=np.around(np.linspace(0,25,6),decimals=1)
            cbar_ax = inset_axes(ax3,width="5%",  height="100%",loc='center right',borderpad=-5)
            #cb3 = fig.colorbar(cs3, spacing='proportional',orientation='vertical',cax=cbar_ax)
            cb3 = fig.colorbar(mappable=None, norm=Normalize(vmin=0,vmax=25), cmap='Greens',spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=Cticks)
            #cb3.formatter.set_powerlimits((0,0))
            #cb2.set_ticklabels(Cticks.astype(int).astype(str),fontsize=7)
            #cb3.set_ticklabels(Cticks.astype(str),fontsize=12)
            cb3.set_label('Water MMV (g/Kg)',fontsize=12)
            #cbar_ax.yaxis.set_label_position('left')
            cbar_ax.yaxis.set_ticks_position('left')
            #ax3.set_ylabel('Hybrid Sigma-Pressure level (mb)',fontsize=10)
            ax3.text(-0.3, 0.5, 'c', va='center', rotation='horizontal',transform=ax3.transAxes,fontsize=16,color='black')
        else:
            cs3 = ax3.contourf(lat,lev,Q_diff,181,cmap='BrBG',vmin=-2,vmax=2)
            c3 = ax3.contour(lat,lev,Q,6,vmin=-0.01,vmax=0.01,colors='black',linewidths=0.5)
            plt.yticks(np.linspace(1000,0,5),fontsize=14)
            plt.gca().invert_yaxis()
            plt.xticks(np.linspace(-90,90,7),fontsize=14)
            if i == 3:
                cticks=np.around(np.linspace(-2,2,5),decimals=2)
                cbar_ax = inset_axes(ax3,width="5%",  height="100%",loc='center right',borderpad=-5)
                #cb = fig.colorbar(cs, spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=cticks)
                cb3 = fig.colorbar(mappable=None, norm=Normalize(vmin=-2,vmax=2), cmap='BrBG',spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=cticks)
                #cb3.formatter.set_powerlimits((0,0))
                #cb3.set_ticklabels(cticks.astype(str),fontsize=12)
                cb3.set_ticklabels(cticks.astype(int).astype(str),fontsize=10)
                cbar_ax.yaxis.set_ticks_position('left')
                cb3.set_label('Difference (g/Kg)',fontsize=10)

        if i >= 1:
            ax4 = plt.subplot(gs[9:12, 3*i+1:3*i+4])
        else:
            ax4 = plt.subplot(gs[9:12, 3*i:3*i+3])
        if i == 0:
            cs4 = ax4.contourf(lat,lev,mastrfu,181,cmap='RdBu',vmin=-0.05,vmax=0.05)
            plt.yticks(np.linspace(1000,0,5),fontsize=14)
            plt.gca().invert_yaxis()
            plt.xticks(np.linspace(-90,90,7),fontsize=14)
            #Cticks=np.around(np.linspace(-0.05,0.05,5),decimals=1)
            cbar_ax = inset_axes(ax4,width="5%",  height="100%",loc='center right',borderpad=-5)
            #cb3 = fig.colorbar(cs3, spacing='proportional',orientation='vertical',cax=cbar_ax)
            cb4 = fig.colorbar(mappable=None, norm=Normalize(vmin=-0.05,vmax=0.05), cmap='RdBu',spacing='proportional',orientation='vertical',cax=cbar_ax)
            #cb4.formatter.set_powerlimits((0,0))
            #cb2.set_ticklabels(Cticks.astype(int).astype(str),fontsize=7)
            #cb3.set_ticklabels(Cticks.astype(str),fontsize=12)
            cb4.set_label('Vertical Velocity (Pa/s)',fontsize=12)
            #cbar_ax.yaxis.set_label_position('left')
            cbar_ax.yaxis.set_ticks_position('left')
            #ax3.set_ylabel('Hybrid Sigma-Pressure level (mb)',fontsize=10)
            ax4.text(-0.3, 0.5, 'd', va='center', rotation='horizontal',transform=ax4.transAxes,fontsize=16,color='black')
        else:
            cs4 = ax4.contourf(lat,lev,deltaChi,181,cmap='RdBu',vmin=-0.01,vmax=0.01)
            c4 = ax4.contour(lat,lev,mastrfu,51,vmin=-0.01,vmax=0.01,colors='black',linewidths=0.5)
            plt.yticks(np.linspace(1000,0,5),fontsize=14)
            plt.gca().invert_yaxis()
            plt.xticks(np.linspace(-90,90,7),fontsize=14)
            if i == 3:
                #cticks=np.around(np.linspace(-5e8,5e8,11),decimals=2)
                cbar_ax = inset_axes(ax4,width="5%",  height="100%",loc='center right',borderpad=-5)
                #cb = fig.colorbar(cs, spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=cticks)
                cb4 = fig.colorbar(mappable=None, norm=Normalize(vmin=-0.05,vmax=0.05), cmap='RdBu_r',spacing='proportional',orientation='vertical',cax=cbar_ax)
                #cb4.formatter.set_powerlimits((0,0))
                #cb3.set_ticklabels(cticks.astype(str),fontsize=12)
                #cb2.set_ticklabels(cticks.astype(int).astype(str),fontsize=10)
                cbar_ax.yaxis.set_ticks_position('left')
                cb4.set_label('Difference (Pa/s)',fontsize=10)
        if i >= 1:
            ax5 = plt.subplot(gs[12:, 3*i+1:3*i+4])
        else:
            ax5 = plt.subplot(gs[12:, 3*i:3*i+3])
        if i == 0:
            cs5 = ax5.contourf(lat,lev,EKE,12,cmap='plasma',vmin=0,vmax=240)
            plt.yticks(np.linspace(1000,0,5),fontsize=14)
            plt.gca().invert_yaxis()
            plt.xticks(np.linspace(-90,90,7),fontsize=14)
            ax5.set_xlabel('Latitude (Deg N)',fontsize=14)
            Cticks=np.around(np.linspace(0,240,9),decimals=1)
            cbar_ax = inset_axes(ax5,width="5%",  height="100%",loc='center right',borderpad=-5)
            cb5 = fig.colorbar(cs5, spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=Cticks)
            #cb.formatter.set_powerlimits((0,0))
            cb5.set_ticklabels(Cticks.astype(int).astype(str),fontsize=12)
            #cb.set_ticklabels(cticks.astype(str),fontsize=10)
            cb5.set_label('EKE ($m^2/s^2$)',fontsize=12)
            #cbar_ax.yaxis.set_label_position('left')
            cbar_ax.yaxis.set_ticks_position('left')
            #ax4.set_ylabel('Hybrid Sigma-Pressure level (mb)',fontsize=10)
            ax5.text(-0.3, 0.5, 'e', va='center', rotation='horizontal',transform=ax5.transAxes,fontsize=16,color='black')
        else:
            cs5 = ax5.contourf(lat,lev,EKE_diff,61,cmap='RdBu_r',vmin=-10,vmax=10)
            c5 = ax5.contour(lat,lev,EKE,12,vmin=0,vmax=240,colors='black',linewidths=0.5)
            plt.yticks(np.linspace(1000,0,5),fontsize=14)
            plt.gca().invert_yaxis()
            plt.xticks(np.linspace(-90,90,7),fontsize=14)
            ax5.set_xlabel('Latitude (Deg N)',fontsize=14)
            if i == 3:
                cticks=np.around(np.linspace(-10,10,5),decimals=2)
                cbar_ax = inset_axes(ax5,width="5%",  height="100%",loc='center right',borderpad=-5)
                #cb = fig.colorbar(cs, spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=cticks)
                cb5 = fig.colorbar(mappable=None, norm=Normalize(vmin=-10,vmax=10), cmap='RdBu_r',spacing='proportional',orientation='vertical',cax=cbar_ax,ticks=cticks)
            #cb.formatter.set_powerlimits((0,0))
            #cb.set_ticklabels(cticks.astype(str),fontsize=7)
                cb5.set_ticklabels(cticks.astype(int).astype(str),fontsize=12)
                cbar_ax.yaxis.set_ticks_position('left')
                cb5.set_label('Difference ($m^2/s^2$)',fontsize=12)
        i = i+1
    else:
        print('No such file or directory: '+dsloc_control+' or '+dsloc)

#    cb.set_label('Temperature (K)')
plt.suptitle('Comparison of Vertical Profiles of Background Climate States',fontsize=24,y=1.01)
fig.text(-0.04, 0.5, 'Hybrid Sigma-Pressure Level (mb)', va='center', rotation='vertical',fontsize=16)
#fig.text(-0.04, 0.5, 'Sigma Pressure Level (mb)', va='center', rotation='vertical')
fig.tight_layout(pad=5)
plt.show()
fig.savefig(figure_path+'Reference_Climate_Comparison_of_Variables.pdf',bbox_inches='tight')

    