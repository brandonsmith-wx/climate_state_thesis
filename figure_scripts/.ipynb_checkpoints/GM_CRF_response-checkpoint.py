#!/usr/bin/env python3

#
# Script for plotting figure 3.22. for altercation, make necessary comments and uncomments
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
FSNT = 'FSNT'
FLNT = 'FLNT'
FSNTC = 'FSNTC'
FLNTC = 'FLNTC'
for CASENAME in casenames:
    run = 'Control'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_control = '/CRFGM_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_control):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo timmean -fldmean -seltimestep,-360/-1 -select,name='+field+','+field2+',FSNT,FSNTC,FLNT,FLNTC '+inpath+infile +' '+ outpath+outfile_control
            print(syscall)
    run = '2xCO2'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'_.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_2xCO2 = '/CRFGM_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_2xCO2):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo timmean -fldmean -seltimestep,-360/-1 -select,name='+field+','+field2+',FSNT,FSNTC,FLNT,FLNTC '+inpath+infile +' '+ outpath+outfile_2xCO2
            print(syscall)
    run = '4xCO2'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'_.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_4xCO2 = '/CRFGM_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_4xCO2):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo timmean -fldmean -seltimestep,-360/-1 -select,name='+field+','+field2+',FSNT,FSNTC,FLNT,FLNTC '+inpath+infile +' '+ outpath+outfile_4xCO2
            print(syscall)
# ONCE LINES ARE PRINTED, COPY OVER TO A SHELL SCRIPT ON A SERVER WHERE CLIMATE DATA OPERATORS IS INSTALLED AND RUN TO CREATE THE FILES.
# ALTERNATIVELY, REPLACE print() STATEMENTS WITH os.system() FUNCTION CALL TO CALL DIRECTLY TO COMMAND LINE.


# Load variables and perform calculations from outfiles
i = 0
outfile_1C = '/home/brandonsmith/modeloutput/Control/1.0/CRFGM_Control_1.0.nc'
outfile_1D = '/home/brandonsmith/modeloutput/2xCO2/1.0/CRFGM_2xCO2_1.0.nc'
colors = ['b','g','k','r']
swcf = []
swcf_d = []
swcf_q = []
lwcf = []
lwcf_d = []
lwcf_q = []
swcf_2x = []
swcf_4x = []
lwcf_2x = []
lwcf_4x = []
swcf_2x_norm = []
swcf_4x_norm = []
lwcf_2x_norm = []
lwcf_4x_norm = []
net = []
net_d = []
net_q = []
net_2x = []
net_4x = []
net_2x_norm = []
net_4x_norm = []
scre = []
lcre = []
delta_s = []
delta_l = []
cn = []
rf = [5.53394852, 4.72283731, 3.88779531, 2.78457941] #defined radiative forcing from doubling CO2 (Byrne and Goldblatt (2014))
rfq = [11.4426504, 9.82042797, 8.15034396, 5.94391217] #defined radiative forcing from quadrupling CO2 (Byrne and Goldblatt (2014))
for CASENAME in casenames:
    outfile_2xCO2 = '/CRFGM_2xCO2_'+CASENAME+'.nc'
    outfile_4xCO2 = '/CRFGM_4xCO2_'+CASENAME+'.nc'
    outfile_control = '/CRFGM_Control_'+CASENAME+'.nc'
    outfile2_control = '/CRFGM_Control_'+CASENAME+'.nc'
    outfile2_2xCO2 = '/CRFGM_2xCO2_'+CASENAME+'.nc'
    outfile2_4xCO2 = '/CRFGM_4xCO2_'+CASENAME+'.nc'
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

        dsc = netCDF4.Dataset(dsloc_control)
        Cvar= dsc.variables[field][:]
        Cvar2 = dsc.variables[field2][:]
        SWAS = dsc.variables[FSNT][:]
        LWAS = dsc.variables[FLNT][:]
        SWCS = dsc.variables[FSNTC][:]
        LWCS = dsc.variables[FLNTC][:]
        #lev = dsc.variables['lev'][:]
        dsc.close()

        dsd = netCDF4.Dataset(dsloc_2xCO2)
        Dvar = dsd.variables[field][:]
        Dvar2 = dsd.variables[field2][:]
        SWAS_2x = dsd.variables[FSNT][:]
        LWAS_2x = dsd.variables[FLNT][:]
        SWCS_2x = dsd.variables[FSNTC][:]
        LWCS_2x = dsd.variables[FLNTC][:]
        # no need to reload lat/lon
        dsd.close()

        C1 = netCDF4.Dataset(outfile_1C)
        Cvar1 = C1.variables[field][:]
        Cvar1_2 = C1.variables[field2][:]
        C1.close()
        D1 = netCDF4.Dataset(outfile_1D)
        Dvar1 = D1.variables[field][:]
        Dvar1_2 = D1.variables[field2][:]
        D1.close()

        dsq = netCDF4.Dataset(dsloc_4xCO2)
        Qvar = dsq.variables[field][:]
        Qvar2 = dsq.variables[field2][:]
        SWAS_4x = dsq.variables[FSNT][:]
        LWAS_4x = dsq.variables[FLNT][:]
        SWCS_4x = dsq.variables[FSNTC][:]
        LWCS_4x = dsq.variables[FLNTC][:]
        # no need to reload lat/lon
        dsq.close()

        Cvar = Cvar/float(CASENAME)
        Cvar = Cvar.squeeze()
        Cvar2 = Cvar2.squeeze()
        Dvar = Dvar/float(CASENAME)
        Dvar = Dvar.squeeze()
        Dvar2 = Dvar2.squeeze()
        Qvar = Qvar/float(CASENAME)
        Qvar = Qvar.squeeze()
        Qvar2 = Qvar2.squeeze()
#            SWAS = SWAS.squeeze()
#            LWAS = LWAS.squeeze()
#            SWCS = SWCS.squeeze()
#            LWCS = LWCS.squeeze()
#            SWAS_2x = SWAS_2x.squeeze()
#            LWAS_2x = LWAS_2x.squeeze()
#            SWCS_2x = SWCS_2x.squeeze()
#            LWCS_2x = LWCS_2x.squeeze()
#            SWAS_4x = SWAS_4x.squeeze()
#            LWAS_4x = LWAS_4x.squeeze()
#            SWCS_4x = SWCS_4x.squeeze()
#            LWCS_4x = LWCS_4x.squeeze()

#            SCRE = SWAS - SWCS
#            LCRE = LWCS - LWAS
#            SCRE_2x = SWAS_2x - SWCS_2x
#            LCRE_2x = LWAS_2x - LWCS_2x
#            SCRE_4x = SWAS_4x - SWCS_4x
#            LCRE_4x = LWAS_4x - LWCS_4x


        diffsn = (Dvar - Cvar)/rf[i]
        Qdiffsn = (Qvar - Cvar)/rfq[i]
        diffln = (Dvar2 - Cvar2)/rf[i]
        Qdiffln = (Qvar2 - Cvar2)/rfq[i]
        diffsn= diffsn.squeeze().squeeze()
        Qdiffsn = Qdiffsn.squeeze().squeeze()
        netCRF_2xn = diffsn + diffln
        netCRF_4xn = Qdiffsn + Qdiffln
        swcf_2x_norm.append(diffsn)
        swcf_4x_norm.append(Qdiffsn)
        lwcf_2x_norm.append(diffln)
        lwcf_4x_norm.append(Qdiffln)
        net_2x_norm.append(netCRF_2xn)
        net_4x_norm.append(netCRF_4xn)
        diffs = (Dvar - Cvar)
        Qdiffs = (Qvar - Cvar)
        diffl = (Dvar2 - Cvar2)
        Qdiffl = (Qvar2 - Cvar2)
        diffs = diffs.squeeze().squeeze()
        Qdiffs = Qdiffs.squeeze().squeeze()
        netCRF_2x = diffs + diffl
        netCRF_4x = Qdiffs + Qdiffl
        swcf_2x.append(diffs)
        swcf_4x.append(Qdiffs)
        lwcf_2x.append(diffl)
        lwcf_4x.append(Qdiffl)
        net_2x.append(netCRF_2x)
        net_4x.append(netCRF_4x)

        netCRF = Cvar + Cvar2
        netCRF_d = Dvar + Dvar2
        netCRF_q = Qvar + Qvar2
        #diff = diff - (Dvar1 - Cvar1)
        swcf.append(Cvar)
        swcf_d.append(Dvar)
        swcf_q.append(Qvar)
        lwcf.append(Cvar2)
        lwcf_d.append(Dvar2)
        lwcf_q.append(Qvar2)
        net.append(netCRF)
        net_d.append(netCRF_d)
        net_q.append(netCRF_q)
#            scre.append(SCRE)
#            lcre.append(LCRE)
#            delta_s.append(diffS)
#            delta_l.append(diffL)
        cn.append(float(casenames[i]))

        i+=1
    else:
        print('No such file or directory')
crf = []
crf.append(swcf)
crf.append(lwcf)
crf.append(net)
dcrf = []
dcrf.append(swcf_d)
dcrf.append(lwcf_d)
dcrf.append(net_d)
qcrf = []
qcrf.append(swcf_q)
qcrf.append(lwcf_q)
qcrf.append(net_q)
Doubling_response = []
Doubling_response.append(swcf_2x)
Doubling_response.append(lwcf_2x)
Doubling_response.append(net_2x)
Quad_response = []
Quad_response.append(swcf_4x)
Quad_response.append(lwcf_4x)
Quad_response.append(net_4x)
normalized_Doubling_response = []
normalized_Doubling_response.append(swcf_2x_norm)
normalized_Doubling_response.append(lwcf_2x_norm)
normalized_Doubling_response.append(net_2x_norm)
normalized_Quad_response = []
normalized_Quad_response.append(swcf_4x_norm)
normalized_Quad_response.append(lwcf_4x_norm)
normalized_Quad_response.append(net_4x_norm)


fig = plt.figure(figsize=(10,5))
i = 0
Casenames = ['90%','95%','100%','105%']
while i < 3:
    ax = fig.add_subplot(3,3,1+i)
    #plt.plot(cn,ts,'o',color='b')
    ax.plot(cn,crf[i],'.',color='k',label='Control',alpha=0.7)
    ax.plot(cn,dcrf[i],'.',color='b',label='$2xCO_2$',alpha=0.7)
    ax.plot(cn,qcrf[i],'.',color='r',label='$4xCO_2$',alpha=0.7)
    ax.set_xticks(cn)
    #plt.yticks(np.linspace(0,5,6))
    #plt.xticks(cn)
    #plt.ylabel('Forcing ($W/m^2$)')

    if i == 0:
        ax.set_yticks(np.linspace(-55,-40,4))
        ax.set_ylabel('Cloud Forcing $(W/m^2)$')
        ax.set_title('Adjusted Shortwave Cloud Forcing')
    if i == 1:
        ax.set_yticks(np.linspace(18,25,8))
        ax.set_title('Longwave Cloud Forcing')
    if i == 2:
        ax.set_yticks(np.linspace(-25,-20,6))
        ax.set_title('Net Cloud Forcing')

    #ax.legend(['Control','$2xCO_2$','$4xCO_2$'],ncol=3,prop={'size':6},loc='upper center')
    ax.set_xticklabels(Casenames)
    ax.grid(True)

    ax2 = fig.add_subplot(3,3,4+i)
    ax2.plot(cn,Doubling_response[i],'.',color='b',alpha=0.7)
    ax2.plot(cn,Quad_response[i],'.',color='r',alpha=0.7)
    #plt.plot(cn,lwcf_4x,'o',color='cyan')
    #plt.ylim([-5,0])
    #plt.yticks(np.linspace(-5,0,6))
    #plt.yticks(np.linspace(0.25,0.5,6))
    ax2.set_xticks(cn)
    #plt.ylabel('Forcing ($W/m^2$)')
    if i == 0:
        ax2.set_ylabel('Difference $(W/m^2)$')
        ax2.set_yticks(np.linspace(0,5,6))
    if i == 1:
        ax2.set_yticks(np.linspace(-4,0,5))
    if i == 2:

        ax2.set_yticks(np.linspace(0,3,4))

    #ax2.legend(['$2xCO_2$','$4xCO_2$'],ncol=2,prop={'size':6},loc='upper center')
    ax2.set_xticklabels(Casenames)
    ax2.grid(True)

    ax3 = fig.add_subplot(3,3,7+i)
    ax3.plot(cn,normalized_Doubling_response[i],'.',color='b',alpha=0.7)
    ax3.plot(cn,normalized_Quad_response[i],'.',color='r',alpha=0.7)
    #plt.plot(cn,lwcf_4x,'o',color='cyan')
    #plt.ylim([-0.2,3])
    #plt.yticks(np.linspace(0,3,4))
    #plt.yticks(np.linspace(0.25,0.5,6))
    ax3.set_xticks(cn)
    #plt.ylabel('Forcing ($W/m^2$)')
    if i == 0:

        ax3.set_yticks(np.linspace(0,0.5,6))
        ax3.set_ylabel('Difference/RF')
    if i == 1:
        ax3.set_yticks(np.linspace(-0.4,0,5))

    if i == 2:
        ax3.set_yticks(np.linspace(-0.05,0.2,6))
        #ax3.set_ylim([-0.05,0.2])
    ax3.set_xlabel('Solar Constant Multiplier')
    #ax3.legend(['$2xCO_2$','$4xCO_2$'],ncol=2,prop={'size':6},loc='upper center')
    ax3.set_xticklabels(Casenames)
    ax3.grid(True)
    i+=1
###################################################################################################################################################
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles,labels,loc = (0.35, 0.86), ncol=5 )
fig.suptitle('Global Mean Cloud Radiative Forcing Comparison',y=1.02,fontsize=20)
fig.tight_layout()
plt.show()
fig.savefig(figure_path+'Global_mean_CRF.pdf',bbox_inches='tight')


