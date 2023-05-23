#!/usr/bin/env python3

#
# Script for plotting figure 3.19. for altercation, make necessary comments and uncomments
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

field = 'TGCLDLWP'
field2 = 'TGCLDIWP'
for CASENAME in casenames:
    run = 'Control'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME #change to location of parent data set
    infile = '/Merged_'+run+'_'+CASENAME+'.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_control = '/cloud_phase_GM_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_control):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo fldmean -seltimestep,-360/-1 -select,name='+field+','+field2+',TGCLDCWP '+inpath+infile +' '+ outpath+outfile_control
            os.system(syscall)
    run = '2xCO2'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_2xCO2 = '/cloud_phase_GM_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_2xCO2):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo fldmean -seltimestep,-360/-1 -select,name='+field+','+field2+',AREL,AREI,TGCLDCWP '+inpath+infile +' '+ outpath+outfile_2xCO2
            os.system(syscall)
    run = '4xCO2'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_4xCO2 = '/cloud_phase_GM_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_4xCO2):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo fldmean -seltimestep,-360/-1 -select,name='+field+','+field2+',AREL,AREI,TGCLDCWP '+inpath+infile +' '+ outpath+outfile_4xCO2
            os.system(syscall)
# ONCE LINES ARE PRINTED, COPY OVER TO A SHELL SCRIPT ON A SERVER WHERE CLIMATE DATA OPERATORS IS INSTALLED AND RUN TO CREATE THE FILES.
# ALTERNATIVELY, REPLACE print() STATEMENTS WITH os.system() FUNCTION CALL TO CALL DIRECTLY TO COMMAND LINE.

# Load variables and perform calculations from outfiles
i = 0
outfile_1C = '/home/brandonsmith/modeloutput/Control/1.0/cloud_phase_GM_Control_1.0.nc'
outfile_1D = '/home/brandonsmith/modeloutput/2xCO2/1.0/cloud_phase_GM_2xCO2_1.0.nc'
colors = ['b','g','k','r']
total = []
total_d = []
total_q = []
total_2xCO2 = []
total_4xCO2 = []
total_2xCO2_norm = []
total_4xCO2_norm = []
total_sem = []
total_dsem = []
total_qsem = []
total_2xCO2_sem = []
total_4xCO2_sem = []
total_2xCO2_norm_sem = []
total_4xCO2_norm_sem = []
difference = []
Liquid = []
Liquid_d = []
Liquid_q = []
Liquid_2xCO2 = []
Liquid_4xCO2 = []
Liquid_2xCO2_norm = []
Liquid_4xCO2_norm = []
Liquid_sem = []
Liquid_dsem = []
Liquid_qsem = []
Liquid_2xCO2_sem = []
Liquid_4xCO2_sem = []
Liquid_2xCO2_norm_sem = []
Liquid_4xCO2_norm_sem = []
Ice = []
Ice_d = []
Ice_q = []
Ice_2xCO2 = []
Ice_4xCO2 = []
Ice_2xCO2_norm = []
Ice_4xCO2_norm = []
Ice_sem = []
Ice_dsem = []
Ice_qsem = []
Ice_2xCO2_sem = []
Ice_4xCO2_sem = []
Ice_2xCO2_norm_sem = []
Ice_4xCO2_norm_sem = []
cn = []
rf = [5.53394852, 4.72283731, 3.88779531, 2.78457941] #defined radiative forcing from doubling CO2 (Byrne and Goldblatt (2014))
rfq = [11.4426504, 9.82042797, 8.15034396, 5.94391217] #defined radiative forcing from quadrupling CO2 (Byrne and Goldblatt (2014))
for CASENAME in casenames:
    outfile_2xCO2 = '/cloud_phase_GM_2xCO2_'+CASENAME+'.nc'
    outfile_4xCO2 = '/cloud_phase_GM_4xCO2_'+CASENAME+'.nc'
    outfile_control = '/cloud_phase_GM_Control_'+CASENAME+'.nc'
    outfile2_control = '/cloud_phase_GM_Control_'+CASENAME+'.nc'
    outfile2_2xCO2 = '/cloud_phase_GM_2xCO2_'+CASENAME+'.nc'
    outfile2_4xCO2 = '/cloud_phase_GM_4xCO2_'+CASENAME+'.nc'
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
        cloud = dsc.variables['TGCLDCWP'][:]
        dsc.close()

        dsd = netCDF4.Dataset(dsloc_2xCO2)
        Dvar = dsd.variables[field][:]
        Dvar2 = dsd.variables[field2][:]
        Dcloud = dsd.variables['TGCLDCWP'][:]
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
        Qcloud = dsq.variables['TGCLDCWP'][:]
        # no need to reload lat/lon
        dsq.close()

        Cvar = Cvar*1000
        Cvar = Cvar.squeeze()
        Cvar2 = Cvar2*1000
        Cvar2 = Cvar2.squeeze()
        Dvar = Dvar*1000
        Dvar = Dvar.squeeze()
        Dvar2 = Dvar2*1000
        Dvar2 = Dvar2.squeeze()
        Qvar = Qvar*1000
        Qvar = Qvar.squeeze()
        Qvar2 = Qvar2*1000
        Qvar2 = Qvar2.squeeze()

        cloud = cloud.squeeze()*1000
        Dcloud = Dcloud.squeeze()*1000
        Qcloud = Qcloud.squeeze()*1000

        Liquid_diff = Dvar - Cvar
        Liquid_Qdiff = Qvar - Cvar
        Ice_diff = Dvar2 - Cvar2
        Ice_Qdiff = Qvar2 - Cvar2
        total_diff = Dcloud - cloud
        total_Qdiff = Qcloud - cloud

        Liquid_diff_norm = (Dvar - Cvar)/rf[i]
        Liquid_Qdiff_norm = (Qvar - Cvar)/rfq[i]
        Ice_diff_norm = (Dvar2 - Cvar2)/rf[i]
        Ice_Qdiff_norm = (Qvar2 - Cvar2)/rfq[i]
        total_diff_norm = (Dcloud - cloud)/rf[i]
        total_Qdiff_norm = (Qcloud - cloud)/rfq[i]
        
        Cvar_sem = sem(Cvar)
        Dvar_sem = sem(Dvar)
        Qvar_sem = sem(Qvar)
        Cvar2_sem = sem(Cvar2)
        Dvar2_sem = sem(Dvar2)
        Qvar2_sem = sem(Qvar2)
        cloud_sem = sem(cloud)
        Dcloud_sem = sem(Dcloud)
        Qcloud_sem = sem(Qcloud)
        
        Liquid_diff_sem = sem(Liquid_diff)
        Liquid_Qdiff_sem = sem(Liquid_Qdiff)
        Ice_diff_sem = sem(Ice_diff)
        Ice_Qdiff_sem = sem(Ice_Qdiff)
        total_diff_sem = sem(total_diff)
        total_Qdiff_sem = sem(total_Qdiff)
        
        Liquid_diff_norm_sem = sem(Liquid_diff_norm)
        Liquid_Qdiff_norm_sem = sem(Liquid_Qdiff_norm)
        Ice_diff_norm_sem = sem(Ice_diff_norm)
        Ice_Qdiff_norm_sem = sem(Ice_Qdiff_norm)
        total_diff_norm_sem = sem(total_diff_norm)
        total_Qdiff_norm_sem = sem(total_Qdiff_norm)
        
        Cvar = np.mean(Cvar)
        Dvar = np.mean(Dvar)
        Qvar = np.mean(Qvar)
        Cvar2 = np.mean(Cvar2)
        Dvar2 = np.mean(Dvar2)
        Qvar2 = np.mean(Qvar2)
        cloud = np.mean(cloud)
        Dcloud = np.mean(Dcloud)
        Qcloud = np.mean(Qcloud)
        Liquid_diff = np.mean(Liquid_diff)
        Liquid_Qdiff = np.mean(Liquid_Qdiff)
        Ice_diff = np.mean(Ice_diff)
        Ice_Qdiff = np.mean(Ice_Qdiff)
        total_diff = np.mean(total_diff)
        total_Qdiff = np.mean(total_Qdiff)
        Liquid_diff_norm = np.mean(Liquid_diff_norm)
        Liquid_Qdiff_norm = np.mean(Liquid_Qdiff_norm)
        Ice_diff_norm = np.mean(Ice_diff_norm)
        Ice_Qdiff_norm = np.mean(Ice_Qdiff_norm)
        total_diff_norm = np.mean(total_diff_norm)
        total_Qdiff_norm = np.mean(total_Qdiff_norm)
        
        total.append(cloud)
        total_d.append(Dcloud)
        total_q.append(Qcloud)
        total_2xCO2.append(total_diff)
        total_4xCO2.append(total_Qdiff)
        total_2xCO2_norm.append(total_diff_norm)
        total_4xCO2_norm.append(total_Qdiff_norm)
        total_sem.append(cloud_sem)
        total_dsem.append(Dcloud_sem)
        total_qsem.append(Qcloud_sem)
        Liquid_sem.append(Cvar_sem)
        Liquid_dsem.append(Dvar_sem)
        Liquid_qsem.append(Qvar_sem)
        Ice_sem.append(Cvar2_sem)
        Ice_dsem.append(Dvar2_sem)
        Ice_qsem.append(Qvar2_sem)
        total_2xCO2_sem.append(total_diff_sem)
        total_4xCO2_sem.append(total_Qdiff_sem)
        total_2xCO2_norm_sem.append(total_diff_norm_sem)
        total_4xCO2_norm_sem.append(total_Qdiff_norm_sem)
        Liquid_2xCO2_sem.append(Liquid_diff_sem)
        Liquid_4xCO2_sem.append(Liquid_Qdiff_sem)
        Liquid_2xCO2_norm_sem.append(Liquid_diff_norm_sem)
        Liquid_4xCO2_norm_sem.append(Liquid_Qdiff_norm_sem)
        Liquid.append(Cvar)
        Liquid_d.append(Dvar)
        Liquid_q.append(Qvar)
        Liquid_2xCO2.append(Liquid_diff)
        Liquid_4xCO2.append(Liquid_Qdiff)
        Liquid_2xCO2_norm.append(Liquid_diff_norm)
        Liquid_4xCO2_norm.append(Liquid_Qdiff_norm)
        Ice.append(Cvar2)
        Ice_d.append(Dvar2)
        Ice_q.append(Qvar2)
        Ice_2xCO2.append(Ice_diff)
        Ice_4xCO2.append(Ice_Qdiff)
        Ice_2xCO2_norm.append(Ice_diff_norm)
        Ice_4xCO2_norm.append(Ice_Qdiff_norm)
        Ice_2xCO2_sem.append(Ice_diff_sem)
        Ice_4xCO2_sem.append(Ice_Qdiff_sem)
        Ice_2xCO2_norm_sem.append(Ice_diff_norm_sem)
        Ice_4xCO2_norm_sem.append(Ice_Qdiff_norm_sem)

        cn.append(float(casenames[i]))

        i+=1
    else:
        print('No such file or directory')

field = []
field.append(Liquid)
field.append(Ice)
field.append(total)
field_sem = []
field_sem.append(Liquid_sem)
field_sem.append(Ice_sem)
field_sem.append(total_sem)
dfield = []
dfield.append(Liquid_d)
dfield.append(Ice_d)
dfield.append(total_d)
dfield_sem = []
dfield_sem.append(Liquid_dsem)
dfield_sem.append(Ice_dsem)
dfield_sem.append(total_dsem)
qfield = []
qfield.append(Liquid_q)
qfield.append(Ice_q)
qfield.append(total_q)
qfield_sem = []
qfield_sem.append(Liquid_qsem)
qfield_sem.append(Ice_qsem)
qfield_sem.append(total_qsem)
Doubling_response = []
Doubling_response.append(Liquid_2xCO2)
Doubling_response.append(Ice_2xCO2)
Doubling_response.append(total_2xCO2)
Doubling_response_sem = []
Doubling_response_sem.append(Liquid_2xCO2_sem)
Doubling_response_sem.append(Ice_2xCO2_sem)
Doubling_response_sem.append(total_2xCO2_sem)
Quad_response = []
Quad_response.append(Liquid_4xCO2)
Quad_response.append(Ice_4xCO2)
Quad_response.append(total_4xCO2)
Quad_response_sem = []
Quad_response_sem.append(Liquid_4xCO2_sem)
Quad_response_sem.append(Ice_4xCO2_sem)
Quad_response_sem.append(total_4xCO2_sem)
normalized_Doubling_response = []
normalized_Doubling_response.append(Liquid_2xCO2_norm)
normalized_Doubling_response.append(Ice_2xCO2_norm)
normalized_Doubling_response.append(total_2xCO2_norm)
normalized_Doubling_response_sem = []
normalized_Doubling_response_sem.append(Liquid_2xCO2_norm_sem)
normalized_Doubling_response_sem.append(Ice_2xCO2_norm_sem)
normalized_Doubling_response_sem.append(total_2xCO2_norm_sem)
normalized_Quad_response = []
normalized_Quad_response.append(Liquid_4xCO2_norm)
normalized_Quad_response.append(Ice_4xCO2_norm)
normalized_Quad_response.append(total_4xCO2_norm)
normalized_Quad_response_sem = []
normalized_Quad_response_sem.append(Liquid_4xCO2_norm_sem)
normalized_Quad_response_sem.append(Ice_4xCO2_norm_sem)
normalized_Quad_response_sem.append(total_4xCO2_norm_sem)

#plt.plot(cn,ts,'o',color='b')
fig = plt.figure(figsize=(10,5))
i = 0
Casenames = ['90%','95%','100%','105%']
while i < 3:
    ax = fig.add_subplot(3,3,i+1)
    #plt.plot(cn,ts,'o',color='b')
    ax.errorbar(cn,field[i],yerr=field_sem[i],fmt='.',capsize=4,color='k',label='Control',alpha=0.7)
    ax.errorbar(cn,dfield[i],yerr=dfield_sem[i],fmt='.',capsize=4,color='b',label='$2xCO_2$',alpha=0.7)
    ax.errorbar(cn,qfield[i],yerr=qfield_sem[i],fmt='.',capsize=4,color='r',label='$4xCO_2$',alpha=0.7)
    ax.set_xticks(cn)
    #plt.yticks(np.linspace(0,5,6))
    #plt.xticks(cn)
    #plt.ylabel('Forcing ($W/m^2$)')
    if i == 0:
        ax.set_title('Vertically Integrated Cloud Liquid')
        ax.set_yticks(np.linspace(34,46,7))
        ax.set_ylabel('LWC ($g/m^2$)')
    if i == 1:
        ax.set_title('Vertically Integrated Cloud Ice',y=1.1)
        ax.set_yticks(np.linspace(10,20,6))
        ax.set_ylabel('IWC ($g/m^2$)')
    if i == 2:
        ax.set_title('Total Vertically Integrated Cloud Water')
        ax.set_yticks(np.linspace(50,60,6))
        ax.set_ylabel('CWC ($g/m^2$)')
        #ax.set_xlabel('Solar Constant Multiplier')
    #ax.legend(['Control','$2xCO_2$','$4xCO_2$'],ncol=3,prop={'size':6},loc='upper center')
    ax.set_xticklabels(Casenames)
    ax.grid(True)

    ax2 = fig.add_subplot(3,3,i+4)
    ax2.errorbar(cn,Doubling_response[i],yerr=Doubling_response_sem[i],fmt='.',capsize=4,color='b',alpha=0.7)
    ax2.errorbar(cn,Quad_response[i],yerr=Quad_response_sem[i],fmt='.',capsize=4,color='r',alpha=0.7)
    #plt.plot(cn,lwcf_4x,'o',color='cyan')
    #plt.ylim([-5,0])
    #plt.yticks(np.linspace(-5,0,6))
    #plt.yticks(np.linspace(0.25,0.5,6))
    ax2.set_xticks(cn)
    #plt.ylabel('Forcing ($W/m^2$)')
    if i == 0:
        ax2.set_title('Difference')
        ax2.set_yticks(np.linspace(1,5,5))
        ax2.set_ylabel('LWC ($g/m^2$)')
    if i == 1:
        ax2.set_title('Difference')
        ax2.set_yticks(np.linspace(-4,0,5))
        ax2.set_ylabel('IWC ($g/m^2$)')
    if i == 2:
        ax2.set_title('Difference')
        ax2.set_ylabel('CWC ($g/m^2$)')
        ax2.set_yticks(np.linspace(-3,3,7))
        #ax2.set_xlabel('Solar Constant Multiplier')
    #ax2.legend(['$2xCO_2$','$4xCO_2$'],ncol=2,prop={'size':6},loc='upper center')
    ax2.set_xticklabels(Casenames)
    ax2.grid(True)

    ax3 = fig.add_subplot(3,3,i+7)
    ax3.errorbar(cn,normalized_Doubling_response[i],yerr=normalized_Doubling_response_sem[i],fmt='.',capsize=4,color='b',alpha=0.7)
    ax3.errorbar(cn,normalized_Quad_response[i],yerr=normalized_Quad_response_sem[i],fmt='.',capsize=4,color='r',alpha=0.7)
    #plt.plot(cn,lwcf_4x,'o',color='cyan')
    #plt.ylim([-0.2,3])
    #plt.yticks(np.linspace(0,3,4))
    #plt.yticks(np.linspace(0.25,0.5,6))
    ax3.set_xticks(cn)

    #plt.ylabel('Forcing ($W/m^2$)')
    if i == 0:
        ax3.set_title('Difference per unit RF')
        ax3.set_yticks(np.linspace(0,0.8,5))
        ax3.set_ylabel('LWC/RF $g/m^2 / W/m^2$')
        ax3.set_xlabel('Solar Constant Multiplier')
    if i == 1:
        ax3.set_yticks(np.linspace(-0.8,0,5))
        ax3.set_title('Difference per unit RF')
        ax3.set_ylabel('IWC/RF ($g/m^2 / W/m^2$)')
        ax3.set_xlabel('Solar Constant Multiplier')
    if i == 2:
        ax3.set_title('Difference per unit RF')
        ax3.set_yticks(np.linspace(-0.3,0.3,7))
        ax3.set_ylabel('CWC/RF ($g/m^2 / W/m^2$)')
        #ax3.set_ylim([-0.05,0.2])
        ax3.set_xlabel('Solar Constant Multiplier')
    #ax3.legend(['$2xCO_2$','$4xCO_2$'],ncol=2,prop={'size':6},loc='upper center')
    ax3.set_xticklabels(Casenames)
    ax3.grid(True)
    i+=1
###################################################################################################################################################
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles,labels,loc = (0.35, 0.86), ncol=5 )
fig.suptitle('Global Mean Cloud Phase comparison and response to Added Radiative Forcing',y=1.02,fontsize=18)
fig.tight_layout()
plt.show()
fig.savefig(figure_path+'Global_mean_cloud_phase_change.pdf',bbox_inches='tight')


