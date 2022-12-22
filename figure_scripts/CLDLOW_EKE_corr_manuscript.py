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
from scipy import stats
from pdf_contours import *
from fastkde import fastKDE


figure_path = '/home/brandonsmith/climate_state_thesis/figures/'
casenames = ['0.9','0.95','1.0','1.05']

field = 'CLDLOW'
Zonal_Mean = False
Global_Mean = False
normalize = True
stf = False
if field == 'CLDLOW' and Zonal_Mean is False and Global_Mean is False:
    for CASENAME in casenames:
        run = 'Control'
        inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        infile = '/Merged_'+run+'_'+CASENAME+'.nc'
        outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        outfile_CLDLOW = '/map_CLDLOW_'+run+'_'+CASENAME+'.nc'
        outfile_EKE = '/map_winds_and_variances_L_'+run+'_'+CASENAME+'.nc'
        
        if os.path.isfile(inpath+infile):
            if not os.path.isfile(outpath+outfile_EKE):
                syscall = '/usr/bin/cdo -sellevidx,21/30 -seltimestep,-360/-1 -select,name=U,V,UU,VV '+inpath+infile+ ' '+outpath+outfile_EKE
                print(syscall)
            if not os.path.isfile(outpath+outfile_CLDLOW):
                syscall = '/usr/bin/cdo timmean -seltimestep,-360/-1 -select,name='+field+' '+inpath+infile+' '+outpath+outfile_CLDLOW
                print(syscall)
        run = '2xCO2'
        inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        infile = '/Merged_'+run+'_'+CASENAME+'_.nc'
        outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        outfile_CLDLOW = '/map_CLDLOW_'+run+'_'+CASENAME+'.nc'
        outfile_EKE = '/map_winds_and_variances_L_'+run+'_'+CASENAME+'.nc'
        
        if os.path.isfile(inpath+infile):
            if not os.path.isfile(outpath+outfile_EKE):
                syscall = '/usr/bin/cdo -sellevidx,21/30 -seltimestep,-360/-1 -select,name=U,V,UU,VV '+inpath+infile+ ' '+outpath+outfile_EKE
                print(syscall)
            if not os.path.isfile(outpath+outfile_CLDLOW):
                syscall = '/usr/bin/cdo timmean -seltimestep,-360/-1 -select,name='+field+' '+inpath+infile+' '+outpath+outfile_CLDLOW
                print(syscall)
        run = '4xCO2'
        inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        infile = '/Merged_'+run+'_'+CASENAME+'_.nc'
        outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        outfile_CLDLOW = '/map_CLDLOW_'+run+'_'+CASENAME+'.nc'
        outfile_EKE = '/map_winds_and_variances_L_'+run+'_'+CASENAME+'.nc'
        
        if os.path.isfile(inpath+infile):
            if not os.path.isfile(outpath+outfile_EKE):
                syscall = '/usr/bin/cdo -sellevidx,21/30 -seltimestep,-360/-1 -select,name=U,V,UU,VV '+inpath+infile+ ' '+outpath+outfile_EKE
                print(syscall)
            if not os.path.isfile(outpath+outfile_CLDLOW):
                syscall = '/usr/bin/cdo timmean -seltimestep,-360/-1 -select,name='+field+' '+inpath+infile+' '+outpath+outfile_CLDLOW
                print(syscall)

    i = 0
    fig = plt.figure(figsize=(10,10),rasterized=True)
    outfile_1C = '/home/brandonsmith/modeloutput/Control/1.0/map_CLDLOW_Control_1.0.nc'
    outfile_1D = '/home/brandonsmith/modeloutput/2xCO2/1.0/map_CLDLOW_2xCO2_1.0.nc'
    outfile_1Q = '/home/brandonsmith/modeloutput/4xCO2/1.0/map_CLDLOW_4xCO2_1.0.nc'
    outfile_EKE_1C = '/home/brandonsmith/modeloutput/Control/1.0/map_winds_and_variances_L_Control_1.0.nc'
    outfile_EKE_1D = '/home/brandonsmith/modeloutput/2xCO2/1.0/map_winds_and_variances_L_2xCO2_1.0.nc'
    outfile_EKE_1Q = '/home/brandonsmith/modeloutput/4xCO2/1.0/map_winds_and_variances_L_4xCO2_1.0.nc'
    rf = [5.53394852, 4.72283731, 3.88779531, 2.78457941]
    rfq = [11.4426504, 9.82042797, 8.15034396, 5.94391217]
    for CASENAME in casenames:
        outfile_2xco2 = '/map_CLDLOW_2xCO2_'+CASENAME+'.nc'
        outfile_control = '/map_CLDLOW_Control_'+CASENAME+'.nc'
        outfile_4xco2 = '/map_CLDLOW_4xCO2_'+CASENAME+'.nc'
        outfile_EKE_2xco2 = '/map_winds_and_variances_L_2xCO2_'+CASENAME+'.nc'
        outfile_EKE_control = '/map_winds_and_variances_L_Control_'+CASENAME+'.nc'
        outfile_EKE_4xco2 = '/map_winds_and_variances_L_4xCO2_'+CASENAME+'.nc'
        outpath_control = '/home/brandonsmith/modeloutput/Control/'+CASENAME
        outpath_2xco2 = '/home/brandonsmith/modeloutput/2xCO2/'+CASENAME
        outpath_4xco2 = '/home/brandonsmith/modeloutput/4xCO2/'+CASENAME
        dsloc_control = outpath_control+outfile_control
        dsloc_2xco2 = outpath_2xco2+outfile_2xco2
        dsloc_4xco2 = outpath_4xco2+outfile_4xco2
        dsloc_EKE_control = outpath_control+outfile_EKE_control
        dsloc_EKE_2xco2 = outpath_2xco2+outfile_EKE_2xco2
        dsloc_EKE_4xco2 = outpath_4xco2+outfile_EKE_4xco2
        if os.path.isfile(dsloc_control) and os.path.isfile(dsloc_EKE_control):
            
            dsc = netCDF4.Dataset(dsloc_control)
            print("Loading Variables")
            Cvar = dsc.variables[field][:]
            lat = dsc.variables['lat'][:]
            lon = dsc.variables['lon'][:]
            dsc.close()
            
            dsd = netCDF4.Dataset(dsloc_2xco2)
            Dvar = dsd.variables[field][:]
            dsd.close()
            
            dsq = netCDF4.Dataset(dsloc_4xco2)
            Qvar = dsq.variables[field][:]
            dsq.close()
            
            C1 = netCDF4.Dataset(outfile_1C)
            Cvar1 = C1.variables[field][:]
            C1.close()
            
            D1 = netCDF4.Dataset(outfile_1D)
            Dvar1 = D1.variables[field][:]
            D1.close()
            
            Q1 = netCDF4.Dataset(outfile_1Q)
            Qvar1 = Q1.variables[field][:]
            Q1.close()
            
            dseke = netCDF4.Dataset(dsloc_EKE_control)
            U = dseke.variables['U'][:]
            V = dseke.variables['V'][:]
            UU = dseke.variables['UU'][:]
            VV = dseke.variables['VV'][:]
            dseke.close()
            
            dsd_eke = netCDF4.Dataset(dsloc_EKE_2xco2)
            Ud = dsd_eke.variables['U'][:]
            Vd = dsd_eke.variables['V'][:]
            UUd = dsd_eke.variables['UU'][:]
            VVd = dsd_eke.variables['VV'][:]
            dsd_eke.close()
            
            dsq_eke = netCDF4.Dataset(dsloc_EKE_4xco2)
            Uq = dsq_eke.variables['U'][:]
            Vq = dsq_eke.variables['V'][:]
            UUq = dsq_eke.variables['UU'][:]
            VVq = dsq_eke.variables['VV'][:]
            dsq_eke.close()
            
            dseke1 = netCDF4.Dataset(outfile_EKE_1C)
            U1 = dseke1.variables['U'][:]
            V1 = dseke1.variables['V'][:]
            UU1 = dseke1.variables['UU'][:]
            VV1 = dseke1.variables['VV'][:]
            dseke1.close()
            
            dsd_eke1 = netCDF4.Dataset(outfile_EKE_1D)
            U1d = dsd_eke1.variables['U'][:]
            V1d = dsd_eke1.variables['V'][:]
            UU1d = dsd_eke1.variables['UU'][:]
            VV1d = dsd_eke1.variables['VV'][:]
            dsd_eke1.close()
            
            dsq_eke1 = netCDF4.Dataset(outfile_EKE_1Q)
            U1q = dsq_eke1.variables['U'][:]
            V1q = dsq_eke1.variables['V'][:]
            UU1q = dsq_eke1.variables['UU'][:]
            VV1q = dsq_eke1.variables['VV'][:]
            dsq_eke1.close()
            
            Cvar = Cvar.squeeze()
            Dvar = Dvar.squeeze()
            Qvar = Qvar.squeeze()
            Cvar1 = Cvar1.squeeze()
            Dvar1 = Dvar1.squeeze()
            Qvar1 = Qvar1.squeeze()
            
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
            
            EKE1 = np.squeeze(EKE1)
            EKE1d = np.squeeze(EKE1d)
            EKE1q = np.squeeze(EKE1q)
            EKE = np.squeeze(EKE)
            EKEd = np.squeeze(EKEd)
            EKEq = np.squeeze(EKEq)
            
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
            if normalize is True:
                EKE_doubling_response = rf[2]*(EKEd-EKE)/rf[i]
                EKE_diff = EKE_doubling_response - (EKE1d-EKE1)
                EKE_Quad_response = rfq[2]*(EKEq-EKE)/rfq[i]
                EKE_Qdiff = EKE_Quad_response - (EKE1q-EKE1)
            else:
                EKE_doubling_response = EKEd - EKE
                EKE_diff = EKE_doubling_response - (EKE1d-EKE1)
                EKE_Quad_response = EKEq - EKE
                EKE_Qdiff = EKE_Quad_response - (EKE1q-EKE1)
            
            doubling_response = doubling_response.flatten()
            EKE_doubling_response = EKE_doubling_response.flatten()
            Quad_response = Quad_response.flatten()
            EKE_Quad_response = EKE_Quad_response.flatten()
            
        
        #Estimate the bivariate pdf 
        mypdf,axis = fastKDE.pdf(EKE_doubling_response,doubling_response)
        v1,v2 = axis
        #Estimate the two marginal pdf
        #mypdf_x,axis_x=fastKDE.pdf(x)
        #mypdf_y,axis_y=fastKDE.pdf(y)
        
        ax = fig.add_subplot(4,2,2*i+1)
        regr = stats.linregress(EKE_doubling_response,doubling_response)
        ax.scatter(EKE_doubling_response,doubling_response,s=1,alpha=0.1)
        ax.contour(v1,v2,mypdf,colors='black',levels=contour_levels_2d(mypdf,v1,v2))
        ax.plot(EKE_doubling_response,regr.intercept + regr.slope*EKE_doubling_response,'r')
        ax.set_xlim([-10,5])
        ax.set_ylim([-0.2,0.2])
        if i == 0:
            ax.set_title('$2xCO_2$ Response',fontsize=14)
        if i == 3:
            ax.set_xlabel('Eddy Kinetic Energy Response',fontsize=10)
        ax.set_ylabel('Cloud Fraction Response',fontsize=8)
        ax.text(0.03, 0.8, 'R = '+str(np.around(regr.rvalue,2)), va='center', rotation='horizontal',transform=ax.transAxes,fontsize=10,color='black')
        plt.grid(True)
        
        mypdf,axis = fastKDE.pdf(EKE_Quad_response,Quad_response)
        v1,v2 = axis
        ax2 = fig.add_subplot(4,2,2*i+2)
        regr = stats.linregress(EKE_Quad_response,Quad_response)
        ax2.scatter(EKE_Quad_response,Quad_response,s=1,alpha=0.1)
        ax2.contour(v1,v2,mypdf,colors='black',levels=contour_levels_2d(mypdf,v1,v2))
        ax2.plot(EKE_Quad_response,regr.intercept + regr.slope*EKE_Quad_response,'r')
        ax2.set_xlim([-20,10])
        ax2.set_ylim([-0.3,0.3])
        if i == 0:
            ax2.set_title('$4xCO_2$ Response',fontsize=14)
        if i == 3:
            ax2.set_xlabel('Eddy Kinetic Energy Response',fontsize=8)
        ax2.text(0.03, 0.8,'R = '+str(np.around(regr.rvalue,2)), va='center', rotation='horizontal',transform=ax2.transAxes,fontsize=10,color='black')
        plt.grid(True)
        
        i+=1
        
    fig.suptitle('Low Cloud Fraction Response and Low-Level Eddy Kinetic Energy Response Regression',fontsize=16)
    plt.show()
    fig.savefig(figure_path+'CLDLOW_EKE_Regression.pdf',bbox_inches='tight')
###################################################################################################################################################