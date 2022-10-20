#!/usr/bin/env python3

#
# Script for plotting figure 3.21. for altercation, make necessary comments and uncomments
# of particular axes settings in the plotting section.
#

import os
import sys
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.colors as colors
figure_path = '/home/brandonsmith/climate_state_thesis/figures/'
casenames = ['0.9','0.95','1.0','1.05']

field = 'SWCF'
field2 = 'LWCF'
for CASENAME in casenames:
    run = 'Control'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_control = '/CRF_RM_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_control):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo timmean -seltimestep,-360/-1 -select,name='+field+','+field2+',AREA '+inpath+infile +' '+ outpath+outfile_control
            print(syscall)
    run = '2xCO2'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'_.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_2xCO2 = '/CRF_RM_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_2xCO2):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo timmean -seltimestep,-360/-1 -select,name='+field+','+field2+',AREA '+inpath+infile +' '+ outpath+outfile_2xCO2
            print(syscall)
    run = '4xCO2'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'_.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_4xCO2 = '/CRF_RM_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_4xCO2):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo timmean -seltimestep,-360/-1 -select,name='+field+','+field2+',AREA '+inpath+infile +' '+ outpath+outfile_4xCO2
            print(syscall)
# ONCE LINES ARE PRINTED, COPY OVER TO A SHELL SCRIPT ON A SERVER WHERE CLIMATE DATA OPERATORS IS INSTALLED AND RUN TO CREATE THE FILES.
# ALTERNATIVELY, REPLACE print() STATEMENTS WITH os.system() FUNCTION CALL TO CALL DIRECTLY TO COMMAND LINE.

# Load variables and perform calculations from outfiles
i = 0
outfile_1C = '/home/brandonsmith/modeloutput/Control/1.0/CRF_RM_Control_1.0.nc'
outfile_1D = '/home/brandonsmith/modeloutput/2xCO2/1.0/CRF_RM_2xCO2_1.0.nc'
outfile_1Q = '/home/brandonsmith/modeloutput/4xCO2/1.0/CRF_RM_4xCO2_1.0.nc'
cn = []
rf = [5.53394852, 4.72283731, 3.88779531, 2.78457941] #defined radiative forcing from doubling CO2 (Byrne and Goldblatt (2014))
rfq = [11.4426504, 9.82042797, 8.15034396, 5.94391217] #defined radiative forcing from quadrupling CO2 (Byrne and Goldblatt (2014))
Tropical = []
Subtropical = []
Mid_latitude = []
Polar = []
for CASENAME in casenames: #Loop over solar cases
    outfile_2xCO2 = '/CRF_RM_2xCO2_'+CASENAME+'.nc'
    outfile_4xCO2 = '/CRF_RM_4xCO2_'+CASENAME+'.nc'
    outfile_control = '/CRF_RM_Control_'+CASENAME+'.nc'
    outfile2_control = '/CRF_RM_Control_'+CASENAME+'.nc'
    outfile2_2xCO2 = '/CRF_RM_2xCO2_'+CASENAME+'.nc'
    outfile2_4xCO2 = '/CRF_RM_4xCO2_'+CASENAME+'.nc'
    outpath_control = '/home/brandonsmith/modeloutput/Control/'+CASENAME
    outpath_2xCO2 = '/home/brandonsmith/modeloutput/2xCO2/'+CASENAME
    outpath_4xCO2 = '/home/brandonsmith/modeloutput/4xCO2/'+CASENAME
    dsloc_control = outpath_control+outfile_control
    dsloc_2xCO2 = outpath_2xCO2+outfile_2xCO2
    dsloc_4xCO2 = outpath_4xCO2+outfile_4xCO2
    dsloc_control_T = outpath_control+outfile2_control
    dsloc_2xCO2_T = outpath_2xCO2+outfile2_2xCO2
    dsloc_4xCO2_T = outpath_4xCO2+outfile2_4xCO2
    if os.path.isfile(dsloc_control) and os.path.isfile(dsloc_2xCO2):
        # Load in Variables
        dsc = netCDF4.Dataset(dsloc_control)
        SWCF = dsc.variables['SWCF'][:]
        LWCF = dsc.variables['LWCF'][:]
        AREA = dsc.variables['AREA'][:]
        lat = dsc.variables['lat'][:]
        lon = dsc.variables['lon'][:]
        dsc.close()

        dsd = netCDF4.Dataset(dsloc_2xCO2)
        SWCFd = dsd.variables['SWCF'][:]
        LWCFd = dsd.variables['LWCF'][:]
        dsd.close()

        dsq = netCDF4.Dataset(dsloc_4xCO2)
        SWCFq = dsq.variables['SWCF'][:]
        LWCFq = dsq.variables['LWCF'][:]
        dsq.close()

        C1 = netCDF4.Dataset(outfile_1C)
        SWCF1 = C1.variables['SWCF'][:]
        LWCF1 = C1.variables['LWCF'][:]
        C1.close()

        D1 = netCDF4.Dataset(outfile_1D)
        SWCF1d = D1.variables['SWCF'][:]
        LWCF1d = D1.variables['LWCF'][:]
        D1.close()

        Q1 = netCDF4.Dataset(outfile_1Q)
        SWCF1q = Q1.variables['SWCF'][:]
        LWCF1q = Q1.variables['LWCF'][:]
        Q1.close()
        
        # adjustment for lower solar constant must be made by dividing by fraction of present day solar constant
        SWCF = SWCF.squeeze()/float(CASENAME)
        SWCFd = SWCFd.squeeze()/float(CASENAME)
        SWCFq = SWCFq.squeeze()/float(CASENAME)
        SWCF1 = SWCF1.squeeze()/float(CASENAME)
        SWCF1d = SWCF1d.squeeze()/float(CASENAME)
        SWCF1q = SWCF1q.squeeze()/float(CASENAME)
        LWCF = LWCF.squeeze()
        LWCFd = LWCFd.squeeze()
        LWCFq = LWCFq.squeeze()
        LWCF1 = LWCF1.squeeze()
        LWCF1d = LWCF1d.squeeze()
        LWCF1q = LWCF1q.squeeze()
        netCRF = SWCF + LWCF
        netCRFd = SWCFd + LWCFd
        netCRFq = SWCFq + LWCFq
        netCRF1 = SWCF1 + LWCF1
        netCRF1d = SWCF1d + LWCF1d
        netCRF1q = SWCF1q + LWCF1q
        netCRF = netCRF.squeeze()
        netCRFd = netCRFd.squeeze()
        netCRFq = netCRFq.squeeze()
        netCRF1 = netCRF1.squeeze()
        netCRF1d = netCRF1d.squeeze()
        netCRF1q = netCRF1q.squeeze()

        Cvar = np.stack((SWCF,LWCF,netCRF))
        Dvar = np.stack((SWCFd,LWCFd,netCRFd))
        Qvar = np.stack((SWCFq,LWCFq,netCRFq))
        Cvar1 = np.stack((SWCF1,LWCF1,netCRF1))
        Dvar1 = np.stack((SWCF1d,LWCF1d,netCRF1d))
        Qvar1 = np.stack((SWCF1q,LWCF1q,netCRF1q))

        Tropical_forcings = []
        Subtropical_forcings = []
        Midlatitude_forcings = []
        Polar_forcings = []
        R_earth = 6371000

        for k in range(0,len(Cvar)): #Loop over shortwave, then longwave, then net
            doubling_response_n = rf[2]*(Dvar[k] - Cvar[k])/rf[i]
            diffn = doubling_response_n - (Dvar1[k] - Cvar1[k])
            Quad_response_n = rfq[2]*(Qvar[k]-Cvar[k])/rfq[i]
            Qdiffn = Quad_response_n - (Qvar1[k]-Cvar1[k])
            doubling_response = Dvar[k] - Cvar[k]
            diff = doubling_response - (Dvar1[k] - Cvar1[k])
            Quad_response = (Qvar[k]-Cvar[k])
            Qdiff = Quad_response - (Qvar1[k]-Cvar1[k])

            AREA = AREA.squeeze()
            doubling_response = doubling_response.squeeze()
            Quad_response = Quad_response.squeeze()
            Cvar[k] = Cvar[k].squeeze()
            Dvar[k] = Dvar[k].squeeze()
            Qvar[k] = Qvar[k].squeeze()
            CRF = np.stack((Cvar[k],Dvar[k],Qvar[k],doubling_response,Quad_response,doubling_response_n,Quad_response_n))

            tropical = []
            subtropical = []
            mid_Latitude = []
            polar = []

            for j in range(0,len(CRF)): #Loop over Control Values, responses, and normalized responses.
                crf = CRF[j].squeeze()
                Tropical_SWCF = crf[40:55,:]
                Subtropical_SWCF_N = crf[32:39,:]
                Mid_Latitude_SWCF_S = crf[22:31,:]
                Polar_SWCF_S = crf[:21,:]
                Subtropical_SWCF_S = crf[56:63,:]
                Mid_Latitude_SWCF_N = crf[64:73,:]
                Polar_SWCF_N = crf[74:,:]

                Mid_latitude_SWCF = np.vstack((Mid_Latitude_SWCF_S,Mid_Latitude_SWCF_N))
                Polar_SWCF = np.vstack((Polar_SWCF_S,Polar_SWCF_N))
                Subtropical_SWCF = np.vstack((Subtropical_SWCF_S,Subtropical_SWCF_N))
                Tropical_SWCF = np.sum(Tropical_SWCF*AREA[38:53,:])/np.sum(AREA[38:53,:])
                Subtropical_SWCF = np.sum(Subtropical_SWCF*(np.vstack((AREA[32:39,:],AREA[56:63,:]))))/np.sum(AREA[32:39,:]+AREA[56:63,:])
                Mid_latitude_SWCF = np.sum(Mid_latitude_SWCF*(np.vstack((AREA[22:31,:],AREA[64:73,:]))))/np.sum(AREA[22:31,:]+AREA[64:73,:])
                Polar_SWCF = np.sum(Polar_SWCF*(np.vstack((AREA[:21,:],AREA[74:,:]))))/np.sum(np.vstack((AREA[:21,:],AREA[74:,:])))


                Tropical_SWCF = Tropical_SWCF.squeeze()
                Subtropical_SWCF = Subtropical_SWCF.squeeze()
                Mid_latitude_SWCF = Mid_latitude_SWCF.squeeze()
                Polar_SWCF = Polar_SWCF.squeeze()


#        add_cyclic_pointA_tropical = 2*np.pi*R_earth*(np.sin(lat[65])-np.sin(lat[31]))
#        A_midlat = 2*np.pi*R_earth*((np.sin(lat[31])-np.sin(lat[17])) + (np.sin(lat[79])-np.sin(lat[65])))
#        A_polar = 2*np.pi*R_earth*((np.sin(lat[17])-np.sin(lat[0])) + (np.sin(lat[95])-np.sin(lat[79])))

                A_tropical = np.sum(AREA[40:55,:])
                A_subtropical = np.sum(AREA[32:39]+AREA[56:63])
                A_midlat = np.sum(AREA[22:31,:]+AREA[64:73,:])
                A_polar = np.sum(np.vstack((AREA[:21,:],AREA[74:,:])))

                A_Earth = 4*np.pi*(R_earth**2)

                Tropical_contribution = np.around((Tropical_SWCF * (A_tropical/A_Earth)),2)
                Subtropical_contribution = np.around((Subtropical_SWCF * (A_subtropical/A_Earth)),2)
                Midlat_contribution = np.around((Mid_latitude_SWCF * (A_midlat/A_Earth)),2)
                Polar_contribution = np.around((Polar_SWCF * (A_polar/A_Earth)),2)

                tropical = np.append(tropical,Tropical_contribution)
                subtropical = np.append(subtropical,Subtropical_contribution)
                mid_Latitude = np.append(mid_Latitude,Midlat_contribution)
                polar = np.append(polar,Polar_contribution)

                field_mean = np.around((np.sum(crf*AREA)/np.sum(AREA))-Polar_contribution,2)

                T_percent = np.around(((Tropical_contribution / field_mean)*100),0)
                S_percent = np.around(((Subtropical_contribution / field_mean)*100),0)
                M_percent = np.around(((Midlat_contribution / field_mean)*100),0)
                P_percent = np.around(((Polar_contribution / field_mean)*100),0)

                tropical = np.array(tropical)
                subtropical = np.array(subtropical)
                mid_Latitude = np.array(mid_Latitude)
                polar = np.array(polar)

            if k == 0:
                Tropical_forcings = np.append(Tropical_forcings,tropical)
                Subtropical_forcings = np.append(Subtropical_forcings,subtropical)
                Midlatitude_forcings = np.append(Midlatitude_forcings,mid_Latitude)
                Polar_forcings = np.append(Polar_forcings,polar)
            else:
                Tropical_forcings = np.vstack((Tropical_forcings,tropical))
                Subtropical_forcings = np.vstack((Subtropical_forcings,subtropical))
                Midlatitude_forcings = np.vstack((Midlatitude_forcings,mid_Latitude))
                Polar_forcings = np.vstack((Polar_forcings,polar))

            print(CASENAME+' Tropical Contribution: '+str(T_percent))
            print(CASENAME+' Subtropical Contribution: '+str(S_percent))
            print(CASENAME+' Mid-Latitude Contribution: '+str(M_percent))
            print(CASENAME+' High-Latitude Contribution: '+str(P_percent))

        if i == 0:
            Tropical = Tropical_forcings.swapaxes(0,1)
            Subtropical = Subtropical_forcings.swapaxes(0,1)
            Mid_latitude = Midlatitude_forcings.swapaxes(0,1)
            Polar = Polar_forcings.swapaxes(0,1)
        elif i == 1:
            Tropical = np.array([Tropical,Tropical_forcings.swapaxes(0,1)]).swapaxes(0,2)
            Subtropical = np.array([Subtropical,Subtropical_forcings.swapaxes(0,1)]).swapaxes(0,2)
            Mid_latitude = np.array([Mid_latitude,Midlatitude_forcings.swapaxes(0,1)]).swapaxes(0,2)
            Polar = np.array([Polar,Polar_forcings.swapaxes(0,1)]).swapaxes(0,2)
        else:
            Tropical = np.dstack((Tropical,Tropical_forcings))
            Subtropical = np.dstack((Subtropical,Subtropical_forcings))
            Mid_latitude = np.dstack((Mid_latitude,Midlatitude_forcings))
            Polar = np.dstack((Polar,Polar_forcings))
    i+=1

fig = plt.figure(figsize=(20,10))
Casenames = ['90%','95%','100%','105%']
i = 0
while i < 3:
    ax = fig.add_subplot(3,3,1+i)
    #plt.plot(cn,ts,'o',color='b')
    ax.plot(casenames,Tropical[i,0,:].squeeze(),'o',color='r',linestyle='solid',label='Tropics (0$^o$-14.5$^o$)')
    ax.plot(casenames,Subtropical[i,0,:].squeeze(),'o',color='orange',linestyle='solid',label='Subtropics (14.5$^o$-30$^o$)')
    ax.plot(casenames,Mid_latitude[i,0,:].squeeze(),'o',color='g',linestyle='solid',label='Mid-Latitudes (30$^o$-48.6$^o$)')
    ax.plot(casenames,Polar[i,0,:].squeeze(),'o',color='b',linestyle='solid',label='High-Latitudes (48.6$^o$-90$^o$)')
    ax.plot(casenames,Tropical[i,1,:].squeeze(),'x',color='r',linestyle='dashed')
    ax.plot(casenames,Subtropical[i,1,:].squeeze(),'x',color='orange',linestyle='dashed')
    ax.plot(casenames,Mid_latitude[i,1,:].squeeze(),'x',color='g',linestyle='dashed')
    ax.plot(casenames,Polar[i,1,:].squeeze(),'x',color='b',linestyle='dashed')
    ax.plot(casenames,Tropical[i,2,:].squeeze(),'+',color='r',linestyle='dotted')
    ax.plot(casenames,Subtropical[i,2,:].squeeze(),'+',color='orange',linestyle='dotted')
    ax.plot(casenames,Mid_latitude[i,2,:].squeeze(),'+',color='g',linestyle='dotted')
    ax.plot(casenames,Polar[i,2,:].squeeze(),'+',color='b',linestyle='dotted')
    #ax.set_xticks(cn)
    #plt.yticks(np.linspace(0,5,6))
    #plt.xticks(cn)
    #plt.ylabel('Forcing ($W/m^2$)')

    if i == 0:
        yticks = np.linspace(-15,-7,5)
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks.astype(int).astype(str),fontsize=14)
        ax.set_ylabel('Cloud Forcing $(W/m^2)$',fontsize=16)
        ax.set_title('Adjusted Shortwave Cloud Forcing',fontsize=16)
    if i == 1:
        yticks = np.linspace(3,9,7)
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks.astype(int).astype(str),fontsize=14)
        ax.set_title('Longwave Cloud Forcing',fontsize=16)
    if i == 2:
        yticks = np.linspace(-9,-3,7)
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks.astype(int).astype(str),fontsize=14)
        ax.set_title('Net Cloud Forcing',fontsize=16)

    #ax.legend(['Control','$2xCO_2$','$4xCO_2$'],ncol=3,prop={'size':6},loc='upper center')
    #ax.set_xticklabels(Casenames,fontsize=16)
    ax.get_xaxis().set_ticks([])
    ax.grid(True)

    ax2 = fig.add_subplot(3,3,4+i)
    ax2.plot(casenames,Tropical[i,3,:].squeeze(),'x',color='r',linestyle='dashed')
    ax2.plot(casenames,Subtropical[i,3,:].squeeze(),'x',color='orange',linestyle='dashed')
    ax2.plot(casenames,Mid_latitude[i,3,:].squeeze(),'x',color='g',linestyle='dashed')
    ax2.plot(casenames,Polar[i,3,:].squeeze(),'x',color='b',linestyle='dashed')
    ax2.plot(casenames,Tropical[i,4,:].squeeze(),'+',color='r',linestyle='dotted')
    ax2.plot(casenames,Subtropical[i,4,:].squeeze(),'+',color='orange',linestyle='dotted')
    ax2.plot(casenames,Mid_latitude[i,4,:].squeeze(),'+',color='g',linestyle='dotted')
    ax2.plot(casenames,Polar[i,4,:].squeeze(),'+',color='b',linestyle='dotted')
    #plt.plot(cn,lwcf_4x,'o',color='cyan')
    #plt.ylim([-5,0])
    #plt.yticks(np.linspace(-5,0,6))
    #plt.yticks(np.linspace(0.25,0.5,6))
    #ax2.set_xticks(cn)
    #plt.ylabel('Forcing ($W/m^2$)')
    if i == 0:
        ax2.set_ylabel('Response $(W/m^2)$',fontsize=16)
        yticks = np.linspace(-1,3,5)
        ax2.set_yticks(yticks)
        ax2.set_yticklabels(yticks.astype(int).astype(str),fontsize=14)
    if i == 1:
        yticks = np.linspace(-2,1,4)
        ax2.set_yticks(yticks)
        ax2.set_yticklabels(yticks.astype(int).astype(str),fontsize=14)
    if i == 2:
        yticks = np.linspace(-1,2,4)
        ax2.set_yticks(yticks)
        ax2.set_yticklabels(yticks.astype(int).astype(str),fontsize=14)

    #ax2.legend(['$2xCO_2$','$4xCO_2$'],ncol=2,prop={'size':6},loc='upper center')
    #ax2.set_xticklabels(Casenames,fontsize=16)
    ax2.get_xaxis().set_ticks([])
    ax2.grid(True)
    ax3 = fig.add_subplot(3,3,7+i)
    ax3.plot(casenames,Tropical[i,5,:].squeeze(),'x',color='r',linestyle='dashed')
    ax3.plot(casenames,Subtropical[i,5,:].squeeze(),'x',color='orange',linestyle='dashed')
    ax3.plot(casenames,Mid_latitude[i,5,:].squeeze(),'x',color='g',linestyle='dashed')
    ax3.plot(casenames,Polar[i,5,:].squeeze(),'x',color='b',linestyle='dashed')
    ax3.plot(casenames,Tropical[i,6,:].squeeze(),'+',color='r',linestyle='dotted')
    ax3.plot(casenames,Subtropical[i,6,:].squeeze(),'+',color='orange',linestyle='dotted')
    ax3.plot(casenames,Mid_latitude[i,6,:].squeeze(),'+',color='g',linestyle='dotted')
    ax3.plot(casenames,Polar[i,6,:].squeeze(),'+',color='b',linestyle='dotted')
    #plt.plot(cn,lwcf_4x,'o',color='cyan')
    #plt.ylim([-0.2,3])
    #plt.yticks(np.linspace(0,3,4))
    #plt.yticks(np.linspace(0.25,0.5,6))
    #ax3.set_xticks(cn)
    #plt.ylabel('Forcing ($W/m^2$)')
    if i == 0:
        ax3.set_ylabel('Normalized Response $(W/m^2)$',fontsize=16)
        yticks = np.linspace(-2,3,6)
        ax3.set_yticks(yticks)
        ax3.set_yticklabels(yticks.astype(int).astype(str),fontsize=14)
    if i == 1:
        yticks = np.linspace(-2,1,4)
        ax3.set_yticks(yticks)
        ax3.set_yticklabels(yticks.astype(int).astype(str),fontsize=14)
    if i == 2:
        yticks = np.linspace(-1,2,4)
        ax3.set_yticks(yticks)
        ax3.set_yticklabels(yticks.astype(int).astype(str),fontsize=14)
        #ax3.set_ylim([-0.05,0.2])
    ax3.set_xlabel('Solar Constant Multiplier',fontsize=16)
    #ax3.legend(['$2xCO_2$','$4xCO_2$'],ncol=2,prop={'size':6},loc='upper center')
    ax3.set_xticklabels(Casenames,fontsize=16)
    ax3.grid(True)
    i+=1
###################################################################################################################################################
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles,labels,loc = (0.25, 0.92), ncol=4, fontsize=12)
fig.suptitle('Global Mean Cloud Radiative Forcing Contributions by Latitude Band',y=1.02,fontsize=24)
fig.tight_layout()
plt.show()
fig.savefig(figure_path+'Global_mean_CRF_contributions.pdf',bbox_inches='tight')


