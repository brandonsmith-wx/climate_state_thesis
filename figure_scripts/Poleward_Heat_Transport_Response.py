#!/usr/bin/env python3
###################################################################################################################################################

#
# Script for plotting figure 3.7. for altercation, make necessary comments and uncomments
# of particular axes settings in the plotting section. Code from Brian Rose, University of Albany, Lecture 13 (2015)
# 
# https://www.atmos.albany.edu/facstaff/brose/classes/ATM623_Spring2015/Notes/Lectures/Lecture13%20--%20Heat%20transport.html#section7
#
# REQUIRES PYTHON PACKAGES: climlab, pooch, xarray, future


import os
import sys
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import climlab
from climlab import constants as const

casenames = ['0.9','0.95','1.0','1.05']
figure_path = '/home/brandonsmith/climate-gcm-bps/plots/'
normalize = True

# For creating smaller file for plotting. Opening parent model output data is too bulky.
for CASENAME in casenames:
        run = 'Control'
        inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        infile = '/Merged_'+run+'_'+CASENAME+'.nc'
        outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        outfile_control = '/ZonalMean_Heat_Transport_'+run+'_'+CASENAME+'.nc'
        outfile_control_precip = '/Precipitation_'+run+'_'+CASENAME+'.nc'
        if not os.path.isfile(outpath+outfile_control) and os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo -timmean -zonmean -seltimestep,-360/-1 -select,name=FSNT,FLNT,LHFLX,SHFLX,FLNS,FSNS,QFLX, '+inpath+infile +' '+ outpath+outfile_control
            print(syscall)
        if not os.path.isfile(outpath+outfile_control_precip) and os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo -timmean -seltimestep,-360/-1 -select,name=PRECSC,PRECSL,PRECC,PRECL '+inpath+infile +' '+ outpath+outfile_control_precip
            print(syscall)
        run = '2xCO2'
        inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        infile = '/Merged_'+run+'_'+CASENAME+'.nc'
        outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        outfile_2xCO2 = '/ZonalMean_Heat_Transport_'+run+'_'+CASENAME+'.nc'
        outfile_2xCO2_precip = '/Precipitation_'+run+'_'+CASENAME+'.nc'
        if not os.path.isfile(outpath+outfile_2xCO2) and os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo -timmean -zonmean -seltimestep,-360/-1 -select,name=FSNT,FLNT,LHFLX,SHFLX,FLNS,FSNS,QFLX, '+inpath+infile +' '+ outpath+outfile_2xCO2
            print(syscall)
        if not os.path.isfile(outpath+outfile_2xCO2_precip) and os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo -timmean -seltimestep,-360/-1 -select,name=PRECSC,PRECSL,PRECC,PRECL '+inpath+infile +' '+ outpath+outfile_2xCO2_precip
            print(syscall)
        run = '4xCO2'
        inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        infile = '/Merged_'+run+'_'+CASENAME+'.nc'
        outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        outfile_4xCO2 = '/ZonalMean_Heat_Transport_'+run+'_'+CASENAME+'.nc'
        outfile_4xCO2_precip = '/Precipitation_'+run+'_'+CASENAME+'.nc'
        if not os.path.isfile(outpath+outfile_4xCO2) and os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo -timmean -zonmean -seltimestep,-360/-1 -select,name=FSNT,FLNT,LHFLX,SHFLX,FLNS,FSNS,QFLX, '+inpath+infile +' '+ outpath+outfile_4xCO2
            print(syscall)
        if not os.path.isfile(outpath+outfile_4xCO2_precip) and os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo -timmean -seltimestep,-360/-1 -select,name=PRECSC,PRECSL,PRECC,PRECL '+inpath+infile +' '+ outpath+outfile_4xCO2_precip
            print(syscall)
# ONCE LINES ARE PRINTED, COPY OVER TO A SHELL SCRIPT ON A SERVER WHERE CLIMATE DATA OPERATORS IS INSTALLED AND RUN TO CREATE THE FILES.
# ALTERNATIVELY, REPLACE print() STATEMENTS WITH os.system() FUNCTION CALL TO CALL DIRECTLY TO COMMAND LINE.

# Read Plotting Files
i = 0
outfile_1C = '/home/brandonsmith/modeloutput/Control/1.0/ZonalMean_Heat_Transport_Control_1.0.nc'
outfile_1D = '/home/brandonsmith/modeloutput/2xCO2/1.0/ZonalMean_Heat_Transport_2xCO2_1.0.nc'
outfile_1Q = '/home/brandonsmith/modeloutput/4xCO2/1.0/ZonalMean_Heat_Transport_4xCO2_1.0.nc'
outfile_1PC = '/home/brandonsmith/modeloutput/Control/1.0/Precipitation_Control_1.0.nc'
outfile_1PD = '/home/brandonsmith/modeloutput/2xCO2/1.0/Precipitation_2xCO2_1.0.nc'
outfile_1PQ = '/home/brandonsmith/modeloutput/4xCO2/1.0/Precipitation_4xCO2_1.0.nc'
colors = ['b','g','k','r']
rf = [5.53394852, 4.72283731, 3.88779531, 2.78457941] #defined radiative forcing from doubling CO2 (Byrne and Goldblatt (2014))
rfq = [11.4426504, 9.82042797, 8.15034396, 5.94391217] #defined radiative forcing from quadrupling CO2 (Byrne and Goldblatt (2014))
fig = plt.figure(figsize=(20,16))
for CASENAME in casenames:
    outpath_control = '/home/brandonsmith/modeloutput/Control/'+CASENAME
    outpath_2xCO2 = '/home/brandonsmith/modeloutput/2xCO2/'+CASENAME
    outpath_4xCO2 = '/home/brandonsmith/modeloutput/4xCO2/'+CASENAME
    outfile_control = '/ZonalMean_Heat_Transport_Control_'+CASENAME+'.nc'
    outfile_2xCO2 = '/ZonalMean_Heat_Transport_2xCO2_'+CASENAME+'.nc'
    outfile_4xCO2 = '/ZonalMean_Heat_Transport_4xCO2_'+CASENAME+'.nc'
    outfile_control_precip = '/Precipitation_Control_'+CASENAME+'.nc'
    outfile_2xCO2_precip = '/Precipitation_2xCO2_'+CASENAME+'.nc'
    outfile_4xCO2_precip = '/Precipitation_4xCO2_'+CASENAME+'.nc'
    dsloc_control = outpath_control+outfile_control
    dsloc_2xCO2 = outpath_2xCO2+outfile_2xCO2
    dsloc_4xCO2 = outpath_4xCO2+outfile_4xCO2
    precip_control = outpath_control+outfile_control_precip
    precip_2xCO2 = outpath_2xCO2+outfile_2xCO2_precip
    precip_4xCO2 = outpath_4xCO2+outfile_4xCO2_precip
    
    if os.path.isfile(dsloc_control) and os.path.isfile(dsloc_2xCO2):
        # Load Variables
        dsc = netCDF4.Dataset(dsloc_control)
        dsd = netCDF4.Dataset(dsloc_2xCO2)
        dsq = netCDF4.Dataset(dsloc_4xCO2)
        psc = netCDF4.Dataset(precip_control)
        psd = netCDF4.Dataset(precip_2xCO2)
        psq = netCDF4.Dataset(precip_4xCO2)
        C1 = netCDF4.Dataset(outfile_1C)
        D1 = netCDF4.Dataset(outfile_1D)
        Q1 = netCDF4.Dataset(outfile_1Q)
        PC1 = netCDF4.Dataset(outfile_1PC)
        PD1 = netCDF4.Dataset(outfile_1PD)
        PQ1 = netCDF4.Dataset(outfile_1PQ)
        lat = dsc.variables['lat'][:]
        
        # Function defining calculation for poleward heat transport
        def inferred_heat_transport( energy_in, lat_deg ):
            '''Returns the inferred heat transport (in PW) by integrating the net energy imbalance from pole to pole.'''
            from scipy import integrate
            from climlab import constants as const
            lat_rad = np.deg2rad( lat_deg )
            energy_in = energy_in.squeeze()
            return ( 1E-15 * 2 * np.math.pi * const.a**2 * integrate.cumtrapz( np.cos(lat_rad)*energy_in, x=lat_rad, initial=0. ) )
        
        # Function for applying inferred_heat_tranport() to CESM model output
        def CESM_heat_transport(ncdata,ncdata2):
            import numpy as np
            import matplotlib.pyplot as plt
            import netCDF4 as nc
            import climlab
            from climlab import constants as const
            lat = ncdata.variables['lat'][:]
            # TOA radiation
            OLR = ncdata.variables['FLNT'][:]
            ASR = ncdata.variables['FSNT'][:]
            Rtoa = ASR - OLR  # net downwelling radiation
            #  surface fluxes  (all positive UP)
            LHF = ncdata.variables['LHFLX'][:]  # latent heat flux (evaporation)
            SHF = ncdata.variables['SHFLX'][:] # sensible heat flux
            LWsfc = ncdata.variables['FLNS'][:]  # net longwave radiation at surface
            SWsfc = -ncdata.variables['FSNS'][:] # net shortwave radiation at surface
            #  energy flux due to snowfall
            SnowFlux =  np.mean(ncdata2.variables['PRECSC'][:]+ncdata2.variables['PRECSL'][:], axis=2)*const.rho_w*const.Lhfus
            #  hydrological cycle
            Evap = np.mean(ncdata.variables['QFLX'][:], axis=2)  # kg/m2/s or mm/s
            Precip = np.mean(ncdata2.variables['PRECC'][:]+ncdata2.variables['PRECL'][:], axis=2)*const.rho_w  # kg/m2/s or mm/s
            EminusP = Evap.squeeze() - Precip.squeeze()  # kg/m2/s or mm/s
            SurfaceRadiation = LWsfc + SWsfc  # net upward radiation from surface
            SurfaceHeatFlux = SurfaceRadiation.squeeze() + LHF.squeeze() + SHF.squeeze() + SnowFlux.squeeze()  # net upward surface heat flux
            Fatmin = Rtoa.squeeze() + SurfaceHeatFlux.squeeze()  # net heat flux in to atmosphere
            
            # heat transport terms
            HT = {}
            HT['total'] = inferred_heat_transport(Rtoa, lat)
            HT['atm'] = inferred_heat_transport(Fatmin, lat)
            HT['ocean'] = inferred_heat_transport(-SurfaceHeatFlux, lat)
            HT['latent'] = inferred_heat_transport(EminusP*const.Lhvap, lat) # atm. latent heat transport from moisture imbal.
            HT['dse'] = HT['atm'] - HT['latent']  # dry static energy transport as residual

            #  annual averages
           # HTann = {}
            #for name, value in HTmonthly.iteritems():
            #    HTann[name] = np.mean(value, axis=0)
        
            return HT        
        
        HTC = CESM_heat_transport(dsc,psc)
        HTD = CESM_heat_transport(dsd,psd)
        HTQ = CESM_heat_transport(dsq,psq)
        HTC1 = CESM_heat_transport(C1,PC1)
        HTD1 = CESM_heat_transport(D1,PD1)
        HTQ1 = CESM_heat_transport(Q1,PQ1)
        
        
        dsc.close()
        dsd.close()
        dsq.close()
        psc.close()
        psd.close()
        psq.close()
        C1.close()
        D1.close()
        Q1.close()
        PC1.close()
        PD1.close()
        PQ1.close()
        
        # Normalize responses with respect to their radiative forcings. If you wish to see raw doubling responses, set "normalize" to False.
        if normalize is True:
            
            HT_2x = {}
            HT_2x['total'] = rf[2]*(HTD['total'] - HTC['total'])/rf[i]
            HT_2x['atm'] = rf[2]*(HTD['atm'] - HTC['atm'])/rf[i]
            HT_2x['ocean'] = rf[2]*(HTD['ocean'] - HTC['ocean'])/rf[i]
            HT_2x['latent'] = rf[2]*(HTD['latent'] - HTC['latent'])/rf[i]
            HT_2x['dse'] = rf[2]*(HTD['dse'] - HTC['dse'])/rf[i]
            
            HT_4x = {}
            HT_4x['total'] = rfq[2]*(HTQ['total'] - HTC['total'])/rfq[i]
            HT_4x['atm'] = rfq[2]*(HTQ['atm'] - HTC['atm'])/rfq[i]
            HT_4x['ocean'] = rfq[2]*(HTQ['ocean'] - HTC['ocean'])/rfq[i]
            HT_4x['latent'] = rfq[2]*(HTQ['latent'] - HTC['latent'])/rfq[i]
            HT_4x['dse'] = rfq[2]*(HTQ['dse'] - HTC['dse'])/rfq[i]
            
            HT_diff = {}
            HT_diff['total'] = HT_2x['total'] - (HTD1['total'] - HTC1['total'])
            HT_diff['atm'] = HT_2x['atm'] - (HTD1['atm'] - HTC1['atm'])
            HT_diff['ocean'] = HT_2x['ocean'] - (HTD1['ocean'] - HTC1['ocean'])
            HT_diff['latent'] = HT_2x['latent'] - (HTD1['latent'] - HTC1['latent'])
            HT_diff['dse'] = HT_2x['dse'] - (HTD1['dse'] - HTC1['dse'])
            
            HTQ_diff = {}
            HTQ_diff['total'] = HT_4x['total'] - (HTQ1['total'] - HTC1['total'])
            HTQ_diff['atm'] = HT_4x['atm'] - (HTQ1['atm'] - HTC1['atm'])
            HTQ_diff['ocean'] = HT_4x['ocean'] - (HTQ1['ocean'] - HTC1['ocean'])
            HTQ_diff['latent'] = HT_4x['latent'] - (HTQ1['latent'] - HTC1['latent'])
            HTQ_diff['dse'] = HT_4x['dse'] - (HTQ1['dse'] - HTC1['dse'])
            
        else:
            
            HT_2x = {}
            HT_2x['total'] = HTD['total'] - HTC['total']
            HT_2x['atm'] = HTD['atm'] - HTC['atm']
            HT_2x['ocean'] = HTD['ocean'] - HTC['ocean']
            HT_2x['latent'] = HTD['latent'] - HTC['latent']
            HT_2x['dse'] = HTD['dse'] - HTC['dse']
            
            HT_4x = {}
            HT_4x['total'] = HTQ['total'] - HTC['total']
            HT_4x['atm'] = HTQ['atm'] - HTC['atm']
            HT_4x['ocean'] = HTQ['ocean'] - HTC['ocean']
            HT_4x['latent'] = HTQ['latent'] - HTC['latent']
            HT_4x['dse'] = HTQ['dse'] - HTC['dse']
            
            HT_diff = {}
            HT_diff['total'] = HT_2x['total'] - (HTD1['total'] - HTC1['total'])
            HT_diff['atm'] = HT_2x['atm'] - (HTD1['atm'] - HTC1['atm'])
            HT_diff['ocean'] = HT_2x['ocean'] - (HTD1['ocean'] - HTC1['ocean'])
            HT_diff['latent'] = HT_2x['latent'] - (HTD1['latent'] - HTC1['latent'])
            HT_diff['dse'] = HT_2x['dse'] - (HTD1['dse'] - HTC1['dse'])
            
            HTQ_diff = {}
            HTQ_diff['total'] = HT_4x['total'] - (HTQ1['total'] - HTC1['total'])
            HTQ_diff['atm'] = HT_4x['atm'] - (HTQ1['atm'] - HTC1['atm'])
            HTQ_diff['ocean'] = HT_4x['ocean'] - (HTQ1['ocean'] - HTC1['ocean'])
            HTQ_diff['latent'] = HT_4x['latent'] - (HTQ1['latent'] - HTC1['latent'])
            HTQ_diff['dse'] = HT_4x['dse'] - (HTQ1['dse'] - HTC1['dse'])
    
    Casenames = ['90% $S_0$','95% $S_0$','$S_0$','105% $S_0$']
    
    ticks = [-90, -60, -30, 0, 30, 60, 90]
    ax = fig.add_subplot(4,5,5*i+1)
    total = ax.plot(lat, HTC['total'], 'k-', label='total', linewidth=2)
    atm = ax.plot(lat, HTC['atm'], 'r-', label='atm', linewidth=2)
    dry = ax.plot(lat, HTC['dse'], 'r--', label='dry')
    latent = ax.plot(lat, HTC['latent'], 'r:', label='latent')
    ocean = ax.plot(lat, HTC['ocean'], 'b-', label='ocean', linewidth=2)
    
    yticks = np.linspace(-6,6,7)
    ydiffticks = np.linspace(-1.5,1.5,7)
    
    ax.set_xlim(-90,90)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks,fontsize=14)
    ax.text(0.05, 0.9, Casenames[i], va='center', rotation='horizontal',transform=ax.transAxes,fontsize=16,color='black')
    if i == 0:
        ax.set_title('Control Climates',fontsize=16)
    if i == 3:
        ax.set_xlabel('Latitude (Deg N)',fontsize=16)
    ax.set_ylabel('Heat Tranport (PW)',fontsize=16)
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticks,fontsize=14)
    #ax.legend(loc='upper left')
    ax.grid()
    
    ax2 = fig.add_subplot(4, 5, 5*i+2)
    ax2.plot(lat, HT_2x['total'], 'k-', label='total', linewidth=2)
    ax2.plot(lat, HT_2x['atm'], 'r-', label='atm', linewidth=2)
    ax2.plot(lat, HT_2x['dse'], 'r--', label='dry')
    ax2.plot(lat, HT_2x['latent'], 'r:', label='latent')
    ax2.plot(lat, HT_2x['ocean'], 'b-', label='ocean', linewidth=2)

    ax2.set_xlim(-90,90)
    ax2.set_yticks(ydiffticks)
    ax2.set_yticklabels(ydiffticks,fontsize=14)
    if i == 0:
        #ax2.set_title('$F_2$$_x$ Response',fontsize=16)
        ax2.set_title('$CO_2$ Doubling Response',fontsize=16)
    if i == 3:
        ax2.set_xlabel('Latitude (Deg N)',fontsize=16)
    ax2.set_xticks(ticks)
    ax2.set_xticklabels(ticks,fontsize=14)
    #ax2.legend(loc='upper left')
    ax2.grid()
    
    ax3 = fig.add_subplot(4, 5, 5*i+3)
    ax3.plot(lat, HT_diff['total'], 'k-', label='total', linewidth=2)
    ax3.plot(lat, HT_diff['atm'], 'r-', label='atm', linewidth=2)
    ax3.plot(lat, HT_diff['dse'], 'r--', label='dry')
    ax3.plot(lat, HT_diff['latent'], 'r:', label='latent')
    ax3.plot(lat, HT_diff['ocean'], 'b-', label='ocean', linewidth=2)

    ax3.set_xlim(-90,90)
    ax3.set_yticks(ydiffticks)
    ax3.set_yticklabels(ydiffticks,fontsize=14)
    if i == 0:
        ax3.set_title('Difference W.R.T. $S_0$',fontsize=16)
    if i == 3:
        ax3.set_xlabel('Latitude (Deg N)',fontsize=16)
    ax3.set_xticks(ticks)
    ax3.set_xticklabels(ticks,fontsize=14)
    #ax3.legend(loc='upper left')
    ax3.grid()
    
    ax4 = fig.add_subplot(4, 5, 5*i+4)
    ax4.plot(lat, HT_4x['total'], 'k-', label='total', linewidth=2)
    ax4.plot(lat, HT_4x['atm'], 'r-', label='atm', linewidth=2)
    ax4.plot(lat, HT_4x['dse'], 'r--', label='dry')
    ax4.plot(lat, HT_4x['latent'], 'r:', label='latent')
    ax4.plot(lat, HT_4x['ocean'], 'b-', label='ocean', linewidth=2)

    ax4.set_xlim(-90,90)
    ax4.set_yticks(ydiffticks)
    ax4.set_yticklabels(ydiffticks,fontsize=14)
    if i == 0:
        #ax4.set_title('$F_4$$_x$ Response',fontsize=16)
        ax4.set_title('$CO_2$ Quadrupling Response',fontsize=16)
    if i == 3:
        ax4.set_xlabel('Latitude (Deg N)',fontsize=16)
    ax4.set_xticks(ticks)
    ax4.set_xticklabels(ticks,fontsize=14)
    #ax4.legend(loc='upper left')
    ax4.grid()
    
    ax5 = fig.add_subplot(4, 5, 5*i+5)
    ax5.plot(lat, HTQ_diff['total'], 'k-', label='total', linewidth=2)
    ax5.plot(lat, HTQ_diff['atm'], 'r-', label='atm', linewidth=2)
    ax5.plot(lat, HTQ_diff['dse'], 'r--', label='dry')
    ax5.plot(lat, HTQ_diff['latent'], 'r:', label='latent')
    ax5.plot(lat, HTQ_diff['ocean'], 'b-', label='ocean', linewidth=2)

    ax5.set_xlim(-90,90)
    ax5.set_yticks(ydiffticks)
    ax5.set_yticklabels(ydiffticks,fontsize=14)
    if i == 0:
        ax5.set_title('Difference W.R.T. $S_0$',fontsize=16)
    if i == 3:
        ax5.set_xlabel('Latitude (Deg N)',fontsize=16)
    ax5.set_xticks(ticks)
    ax5.set_xticklabels(ticks,fontsize=14)
    #ax5.legend(loc='upper left')
    ax5.grid()

    i +=1
    
fig.suptitle('Poleward Heat Transport Response to $CO_2$ Doubling',fontsize=24)
#fig.suptitle('Poleward Heat Transport Response to Equivalent Radiative Forcing',fontsize=24)

handles, labels = ax.get_legend_handles_labels()
fig.legend(handles,labels,loc = (0.35, 0.92), ncol=5 )

fig.savefig(figure_path+'Poleward_Heat_Transport_2x_Response.pdf',bbox_inches='tight')