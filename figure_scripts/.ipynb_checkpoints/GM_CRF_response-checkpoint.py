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
from scipy.stats import sem
figure_path = '/home/brandonsmith/thesis/figures/'
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
            syscall = '/usr/bin/cdo fldmean -seltimestep,-360/-1 -select,name='+field+','+field2+',FSNT,FSNTC,FLNT,FLNTC '+inpath+infile +' '+ outpath+outfile_control
            os.system(syscall)

    run = '2xCO2'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_2xCO2 = '/CRFGM_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_2xCO2):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo fldmean -seltimestep,-360/-1 -select,name='+field+','+field2+',FSNT,FSNTC,FLNT,FLNTC '+inpath+infile +' '+ outpath+outfile_2xCO2
            os.system(syscall)

    run = '4xCO2'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_4xCO2 = '/CRFGM_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_4xCO2):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo fldmean -seltimestep,-360/-1 -select,name='+field+','+field2+',FSNT,FSNTC,FLNT,FLNTC '+inpath+infile +' '+ outpath+outfile_4xCO2
            os.system(syscall)

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
swcf_sem = []
swcf_dsem = []
swcf_qsem = []
lwcf_sem = []
lwcf_dsem = []
lwcf_qsem = []
net_sem = []
net_dsem = []
net_qsem = []

swcf_doubling_sem = []
swcf_quad_sem = []
norm_swcf_doubling_sem = []
norm_swcf_quad_sem = []

lwcf_doubling_sem = []
lwcf_quad_sem = []
norm_lwcf_doubling_sem = []
norm_lwcf_quad_sem = []

net_doubling_sem = []
net_quad_sem = []
norm_net_doubling_sem = []
norm_net_quad_sem = []
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
        #SWAS = dsc.variables[FSNT][:]
        #LWAS = dsc.variables[FLNT][:]
        #SWCS = dsc.variables[FSNTC][:]
        #LWCS = dsc.variables[FLNTC][:]
        #lev = dsc.variables['lev'][:]
        dsc.close()

        dsd = netCDF4.Dataset(dsloc_2xCO2)
        Dvar = dsd.variables[field][:]
        Dvar2 = dsd.variables[field2][:]
        #SWAS_2x = dsd.variables[FSNT][:]
        #LWAS_2x = dsd.variables[FLNT][:]
        #SWCS_2x = dsd.variables[FSNTC][:]
        #LWCS_2x = dsd.variables[FLNTC][:]
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
        #SWAS_4x = dsq.variables[FSNT][:]
        #LWAS_4x = dsq.variables[FLNT][:]
        #SWCS_4x = dsq.variables[FSNTC][:]
        #LWCS_4x = dsq.variables[FLNTC][:]
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
        
        netCRF = Cvar + Cvar2
        netCRF_d = Dvar + Dvar2
        netCRF_q = Qvar + Qvar2
        
        netCRF_sem = sem(netCRF)
        netCRF_d_sem = sem(netCRF_d)
        netCRF_q_sem = sem(netCRF_q)
        
        Cvar_sem = sem(Cvar)
        Dvar_sem = sem(Dvar)
        Qvar_sem = sem(Qvar)
        Cvar2_sem = sem(Cvar2)
        Dvar2_sem = sem(Dvar2)
        Qvar2_sem = sem(Qvar2)
        
        diffsn = (Dvar - Cvar)/rf[i]
        Qdiffsn = (Qvar - Cvar)/rfq[i]
        diffln = (Dvar2 - Cvar2)/rf[i]
        Qdiffln = (Qvar2 - Cvar2)/rfq[i]
        diffsn= diffsn.squeeze().squeeze()
        Qdiffsn = Qdiffsn.squeeze().squeeze()
        netCRF_2xn = diffsn + diffln
        netCRF_4xn = Qdiffsn + Qdiffln
        
        diffsn_sem = sem(diffsn)
        Qdiffsn_sem = sem(Qdiffsn)
        diffln_sem = sem(diffln)
        Qdiffln_sem = sem(Qdiffln)
        netCRF_2xn_sem = sem(netCRF_2xn)
        netCRF_4xn_sem = sem(netCRF_4xn)
        
        diffs = (Dvar - Cvar)
        Qdiffs = (Qvar - Cvar)
        diffl = (Dvar2 - Cvar2)
        Qdiffl = (Qvar2 - Cvar2)
        diffs = diffs.squeeze().squeeze()
        Qdiffs = Qdiffs.squeeze().squeeze()
        netCRF_2x = diffs + diffl
        netCRF_4x = Qdiffs + Qdiffl
        
        diffs_sem = sem(diffs)
        Qdiffs_sem = sem(Qdiffs)
        diffl_sem = sem(diffl)
        Qdiffl_sem = sem(Qdiffl)
        netCRF_2x_sem = sem(netCRF_2x)
        netCRF_4x_sem = sem(netCRF_4x)
        
        Cvar = np.mean(Cvar)
        Dvar = np.mean(Dvar)
        Qvar = np.mean(Qvar)
        Cvar2 = np.mean(Cvar2)
        Dvar2 = np.mean(Dvar2)
        Qvar2 = np.mean(Qvar2)
        netCRF = np.mean(netCRF)
        netCRF_d = np.mean(netCRF_d)
        netCRF_q = np.mean(netCRF_q)
        
        diffs = np.mean(diffs)
        Qdiffs = np.mean(Qdiffs)
        diffl = np.mean(diffl)
        Qdiffl = np.mean(Qdiffl)
        netCRF_2x = np.mean(netCRF_2x)
        netCRF_4x = np.mean(netCRF_4x)
        
        diffsn = np.mean(diffsn)
        Qdiffsn = np.mean(Qdiffsn)
        diffln = np.mean(diffln)
        Qdiffln = np.mean(Qdiffln)
        netCRF_2xn = np.mean(netCRF_2xn)
        netCRF_4xn = np.mean(netCRF_4xn)
        
        swcf_2x.append(diffs)
        swcf_4x.append(Qdiffs)
        lwcf_2x.append(diffl)
        lwcf_4x.append(Qdiffl)
        net_2x.append(netCRF_2x)
        net_4x.append(netCRF_4x)
        
        swcf_2x_norm.append(diffsn)
        swcf_4x_norm.append(Qdiffsn)
        lwcf_2x_norm.append(diffln)
        lwcf_4x_norm.append(Qdiffln)
        net_2x_norm.append(netCRF_2xn)
        net_4x_norm.append(netCRF_4xn)

        swcf.append(Cvar)
        swcf_d.append(Dvar)
        swcf_q.append(Qvar)
        lwcf.append(Cvar2)
        lwcf_d.append(Dvar2)
        lwcf_q.append(Qvar2)
        net.append(netCRF)
        net_d.append(netCRF_d)
        net_q.append(netCRF_q)
        cn.append(float(casenames[i]))
        
        swcf_sem.append(Cvar_sem)
        swcf_dsem.append(Dvar_sem)
        swcf_qsem.append(Qvar_sem)
        lwcf_sem.append(Cvar2_sem)
        lwcf_dsem.append(Dvar2_sem)
        lwcf_qsem.append(Qvar2_sem)
        net_sem.append(netCRF_sem)
        net_dsem.append(netCRF_d_sem)
        net_qsem.append(netCRF_q_sem)
        swcf_doubling_sem.append(diffs_sem)
        swcf_quad_sem.append(Qdiffs_sem)
        lwcf_doubling_sem.append(diffl_sem)
        lwcf_quad_sem.append(Qdiffl_sem)
        net_doubling_sem.append(netCRF_2x_sem)
        net_quad_sem.append(netCRF_4x_sem)
        norm_swcf_doubling_sem.append(diffsn_sem)
        norm_swcf_quad_sem.append(Qdiffsn_sem)
        norm_lwcf_doubling_sem.append(diffln_sem)
        norm_lwcf_quad_sem.append(Qdiffln_sem)
        norm_net_doubling_sem.append(netCRF_2xn_sem)
        norm_net_quad_sem.append(netCRF_4xn_sem)

        i+=1
    else:
        print('No such file or directory')
crf = []
crf.append(swcf)
crf.append(lwcf)
crf.append(net)
crf_sem = []
crf_sem.append(swcf_sem)
crf_sem.append(lwcf_sem)
crf_sem.append(net_sem)
dcrf = []
dcrf.append(swcf_d)
dcrf.append(lwcf_d)
dcrf.append(net_d)
dcrf_sem = []
dcrf_sem.append(swcf_dsem)
dcrf_sem.append(lwcf_dsem)
dcrf_sem.append(net_dsem)
qcrf = []
qcrf.append(swcf_q)
qcrf.append(lwcf_q)
qcrf.append(net_q)
qcrf_sem = []
qcrf_sem.append(swcf_qsem)
qcrf_sem.append(lwcf_qsem)
qcrf_sem.append(net_qsem)
Doubling_response = []
Doubling_response.append(swcf_2x)
Doubling_response.append(lwcf_2x)
Doubling_response.append(net_2x)
Doubling_sem = []
Doubling_sem.append(swcf_doubling_sem)
Doubling_sem.append(lwcf_doubling_sem)
Doubling_sem.append(net_doubling_sem)
Quad_response = []
Quad_response.append(swcf_4x)
Quad_response.append(lwcf_4x)
Quad_response.append(net_4x)
Quad_sem = []
Quad_sem.append(swcf_quad_sem)
Quad_sem.append(lwcf_quad_sem)
Quad_sem.append(net_quad_sem)
normalized_Doubling_response = []
normalized_Doubling_response.append(swcf_2x_norm)
normalized_Doubling_response.append(lwcf_2x_norm)
normalized_Doubling_response.append(net_2x_norm)
norm_Doubling_sem = []
norm_Doubling_sem.append(norm_swcf_doubling_sem)
norm_Doubling_sem.append(norm_lwcf_doubling_sem)
norm_Doubling_sem.append(norm_net_doubling_sem)
normalized_Quad_response = []
normalized_Quad_response.append(swcf_4x_norm)
normalized_Quad_response.append(lwcf_4x_norm)
normalized_Quad_response.append(net_4x_norm)
norm_Quad_sem = []
norm_Quad_sem.append(norm_swcf_quad_sem)
norm_Quad_sem.append(norm_lwcf_quad_sem)
norm_Quad_sem.append(norm_net_quad_sem)



fig = plt.figure(figsize=(10,5))
i = 0
Casenames = ['90%','95%','100%','105%']
while i < 3:
    print(crf_sem[i])
    ax = fig.add_subplot(3,3,1+i)
    cyerr = crf_sem[i]
    dyerr = dcrf_sem[i]
    qyerr = qcrf_sem[i]
    yerr_2x = Doubling_sem[i]
    yerr_4x = Quad_sem[i]
    yerr_2xn = norm_Doubling_sem[i]
    yerr_4xn = norm_Quad_sem[i]
    #plt.plot(cn,ts,'o',color='b')
    #ax.plot(cn,crf[i],'.',color='k',label='Control',alpha=0.7)
    ax.errorbar(cn,crf[i],yerr=cyerr,fmt='.',capsize=3,color='k',label='Control',alpha=0.7)
    ax.errorbar(cn,dcrf[i],yerr=dyerr,fmt='.',capsize=3,color='b',label='$2xCO_2$',alpha=0.7)
    ax.errorbar(cn,qcrf[i],yerr=qyerr,fmt='.',capsize=3,color='r',label='$4xCO_2$',alpha=0.7)
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
    ax2.errorbar(cn,Doubling_response[i],yerr=yerr_2x,fmt='.',capsize=3,color='b',alpha=0.7)
    ax2.errorbar(cn,Quad_response[i],yerr=yerr_4x,fmt='.',capsize=3,color='r',alpha=0.7)
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
    ax3.errorbar(cn,normalized_Doubling_response[i],yerr=yerr_2xn,fmt='.',capsize=3,color='b',alpha=0.7)
    ax3.errorbar(cn,normalized_Quad_response[i],yerr=yerr_4xn,fmt='.',capsize=3,color='r',alpha=0.7)
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


