#!/usr/bin/env python3

#
# Script for plotting figure 3.20. for altercation, make necessary comments and uncomments
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
field2 = 'CDNUMC'
for CASENAME in casenames:
    run = 'Control'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_control = '/cloud_droplets_GM_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_control):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo fldmean -seltimestep,-360/-1 -select,name='+field+','+field2+',AREA '+inpath+infile +' '+ outpath+outfile_control
            os.system(syscall)
    run = '2xCO2'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_2xCO2 = '/cloud_droplets_GM_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_2xCO2):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo fldmean -seltimestep,-360/-1 -select,name='+field+','+field2+',AREL,AREI,AREA '+inpath+infile +' '+ outpath+outfile_2xCO2
            os.system(syscall)
    run = '4xCO2'
    inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    infile = '/Merged_'+run+'_'+CASENAME+'.nc'
    outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
    outfile_4xCO2 = '/cloud_droplets_GM_'+run+'_'+CASENAME+'.nc'
    if not os.path.isfile(outpath+outfile_4xCO2):
        if os.path.isfile(inpath+infile):
            syscall = '/usr/bin/cdo fldmean -seltimestep,-360/-1 -select,name='+field+','+field2+',AREL,AREI,AREA '+inpath+infile +' '+ outpath+outfile_4xCO2
            os.system(syscall)
# ONCE LINES ARE PRINTED, COPY OVER TO A SHELL SCRIPT ON A SERVER WHERE CLIMATE DATA OPERATORS IS INSTALLED AND RUN TO CREATE THE FILES.
# ALTERNATIVELY, REPLACE print() STATEMENTS WITH os.system() FUNCTION CALL TO CALL DIRECTLY TO COMMAND LINE.

# Load variables and perform calculations from outfiles
i = 0
outfile_1C = '/home/brandonsmith/modeloutput/Control/1.0/cloud_droplets_GM_Control_1.0.nc'
outfile_1D = '/home/brandonsmith/modeloutput/2xCO2/1.0/cloud_droplets_GM_2xCO2_1.0.nc'
colors = ['b','g','k','r']
radius = []
radius_d = []
radius_q = []
radius_2xCO2 = []
radius_4xCO2 = []
radius_2xCO2_norm = []
radius_4xCO2_norm = []
difference = []
CLDLWP = []
CLDLWP_d = []
CLDLWP_q = []
CLDLWP_2xCO2 = []
CLDLWP_4xCO2 = []
CLDLWP_2xCO2_norm = []
CLDLWP_4xCO2_norm = []
CDNUMC = []
CDNUMC_d = []
CDNUMC_q = []
CDNUMC_2xCO2 = []
CDNUMC_4xCO2 = []
CDNUMC_2xCO2_norm = []
CDNUMC_4xCO2_norm = []
radius_sem = []
radius_dsem = []
radius_qsem = []
radius_2xCO2_sem = []
radius_4xCO2_sem = []
radius_2xCO2_norm_sem = []
radius_4xCO2_norm_sem = []
CLDLWP_sem = []
CLDLWP_dsem = []
CLDLWP_qsem = []
CLDLWP_2xCO2_sem = []
CLDLWP_4xCO2_sem = []
CLDLWP_2xCO2_norm_sem = []
CLDLWP_4xCO2_norm_sem = []
CDNUMC_sem = []
CDNUMC_dsem = []
CDNUMC_qsem = []
CDNUMC_2xCO2_sem = []
CDNUMC_4xCO2_sem = []
CDNUMC_2xCO2_norm_sem = []
CDNUMC_4xCO2_norm_sem = []
cn = []
rf = [5.53394852, 4.72283731, 3.88779531, 2.78457941] #defined radiative forcing from doubling CO2 (Byrne and Goldblatt (2014))
rfq = [11.4426504, 9.82042797, 8.15034396, 5.94391217] #defined radiative forcing from quadrupling CO2 (Byrne and Goldblatt (2014))
for CASENAME in casenames:
    outfile_2xCO2 = '/cloud_droplets_GM_2xCO2_'+CASENAME+'.nc'
    outfile_4xCO2 = '/cloud_droplets_GM_4xCO2_'+CASENAME+'.nc'
    outfile_control = '/cloud_droplets_GM_Control_'+CASENAME+'.nc'
    outfile2_control = '/cloud_droplets_GM_Control_'+CASENAME+'.nc'
    outfile2_2xCO2 = '/cloud_droplets_GM_2xCO2_'+CASENAME+'.nc'
    outfile2_4xCO2 = '/cloud_droplets_GM_4xCO2_'+CASENAME+'.nc'
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
        #Droplet_radius_C = dsc.variables['AREL'][:]
        Cvar2 = dsc.variables[field2][:]
        area = dsc.variables['AREA'][:]
        area = area.squeeze()
        #lev = dsc.variables['lev'][:]
        dsc.close()

        dsd = netCDF4.Dataset(dsloc_2xCO2)
        Dvar = dsd.variables[field][:]
        Dvar2 = dsd.variables[field2][:]
        #Droplet_radius_D = dsd.variables['AREL'][:]
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
        # no need to reload lat/lon
        dsq.close()
        rho = 997 #density of liquid water kg/m^3

        #calculate droplet radius (microns)
        Droplet_radius_C = 1000000*np.cbrt(3*Cvar/(4*np.pi*rho*Cvar2))
        Droplet_radius_D = 1000000*np.cbrt(3*Dvar/(4*np.pi*rho*Dvar2))
        Droplet_radius_Q = 1000000*np.cbrt(3*Qvar/(4*np.pi*rho*Qvar2))
        Droplet_radius_C = Droplet_radius_C.squeeze()
        Droplet_radius_D = Droplet_radius_D.squeeze()
        Droplet_radius_Q = Droplet_radius_Q.squeeze()

        Cvar = Cvar*1000
        Cvar = Cvar.squeeze()
        Cvar2 = Cvar2.squeeze()
        Dvar = Dvar*1000
        Dvar = Dvar.squeeze()
        Dvar2 = Dvar2.squeeze()
        Qvar = Qvar*1000
        Qvar = Qvar.squeeze()
        Qvar2 = Qvar2.squeeze()

        droplet_diff = Droplet_radius_D - Droplet_radius_C
        droplet_Qdiff = Droplet_radius_Q - Droplet_radius_C
        CLDLWP_diff = Dvar - Cvar
        CLDLWP_Qdiff = Qvar - Cvar
        CDNUMC_diff = Dvar2 - Cvar2
        CDNUMC_Qdiff = Qvar2 - Cvar2
        droplet_diff_norm = (Droplet_radius_D - Droplet_radius_C)/rf[i]
        droplet_Qdiff_norm = (Droplet_radius_Q - Droplet_radius_C)/rfq[i]
        CLDLWP_diff_norm = (Dvar - Cvar)/rf[i]
        CLDLWP_Qdiff_norm = (Qvar - Cvar)/rfq[i]
        CDNUMC_diff_norm = (Dvar2 - Cvar2)/rf[i]
        CDNUMC_Qdiff_norm = (Qvar2 - Cvar2)/rfq[i]
        
        Cvar_sem = sem(Cvar)
        Dvar_sem = sem(Dvar)
        Qvar_sem = sem(Qvar)
        Cvar2_sem = sem(Cvar2)
        Dvar2_sem = sem(Dvar2)
        Qvar2_sem = sem(Qvar2)
        Droplet_radius_Csem = sem(Droplet_radius_C)
        Droplet_radius_Dsem = sem(Droplet_radius_D)
        Droplet_radius_Qsem = sem(Droplet_radius_Q)
        droplet_diff_sem = sem(droplet_diff)
        droplet_Qdiff_sem = sem(droplet_Qdiff)
        CLDLWP_diff_sem = sem(CLDLWP_diff)
        CLDLWP_Qdiff_sem = sem(CLDLWP_Qdiff)
        CDNUMC_diff_sem = sem(CDNUMC_diff)
        CDNUMC_Qdiff_sem = sem(CDNUMC_Qdiff)
        droplet_diff_norm_sem = sem(droplet_diff_norm)
        droplet_Qdiff_norm_sem = sem(droplet_Qdiff_norm)
        CLDLWP_diff_norm_sem = sem(CLDLWP_diff_norm)
        CLDLWP_Qdiff_norm_sem = sem(CLDLWP_Qdiff_norm)
        CDNUMC_diff_norm_sem = sem(CDNUMC_diff_norm)
        CDNUMC_Qdiff_norm_sem = sem(CDNUMC_Qdiff_norm)
        
        Cvar = np.mean(Cvar)
        Dvar = np.mean(Dvar)
        Qvar = np.mean(Qvar)
        Cvar2 = np.mean(Cvar2)
        Dvar2 = np.mean(Dvar2)
        Qvar2 = np.mean(Qvar2)
        Droplet_radius_C = np.mean(Droplet_radius_C)
        Droplet_radius_D = np.mean(Droplet_radius_D)
        Droplet_radius_Q = np.mean(Droplet_radius_Q)
        droplet_diff = np.mean(droplet_diff)
        droplet_Qdiff = np.mean(droplet_Qdiff)
        CLDLWP_diff = np.mean(CLDLWP_diff)
        CLDLWP_Qdiff = np.mean(CLDLWP_Qdiff)
        CDNUMC_diff = np.mean(CDNUMC_diff)
        CDNUMC_Qdiff = np.mean(CDNUMC_Qdiff)
        droplet_diff_norm = np.mean(droplet_diff_norm)
        droplet_Qdiff_norm = np.mean(droplet_Qdiff_norm)
        CLDLWP_diff_norm = np.mean(CLDLWP_diff_norm)
        CLDLWP_Qdiff_norm = np.mean(CLDLWP_Qdiff_norm)
        CDNUMC_diff_norm = np.mean(CDNUMC_diff_norm)
        CDNUMC_Qdiff_norm = np.mean(CDNUMC_Qdiff_norm)

        #diff = diff - (Dvar1 - Cvar1)
        radius.append(Droplet_radius_C)
        radius_d.append(Droplet_radius_D)
        radius_q.append(Droplet_radius_Q)
        radius_2xCO2.append(droplet_diff)
        radius_4xCO2.append(droplet_Qdiff)
        radius_2xCO2_norm.append(droplet_diff_norm)
        radius_4xCO2_norm.append(droplet_Qdiff_norm)
        radius_sem.append(Droplet_radius_Csem)
        radius_dsem.append(Droplet_radius_Dsem)
        radius_qsem.append(Droplet_radius_Qsem)
        radius_2xCO2_sem.append(droplet_diff_sem)
        radius_4xCO2_sem.append(droplet_Qdiff_sem)
        radius_2xCO2_norm_sem.append(droplet_diff_norm_sem)
        radius_4xCO2_norm_sem.append(droplet_Qdiff_norm_sem)
        CLDLWP.append(Cvar)
        CLDLWP_d.append(Dvar)
        CLDLWP_q.append(Qvar)
        CLDLWP_2xCO2.append(CLDLWP_diff)
        CLDLWP_4xCO2.append(CLDLWP_Qdiff)
        CLDLWP_2xCO2_norm.append(CLDLWP_diff_norm)
        CLDLWP_4xCO2_norm.append(CLDLWP_Qdiff_norm)
        CDNUMC.append(Cvar2)
        CDNUMC_d.append(Dvar2)
        CDNUMC_q.append(Qvar2)
        CDNUMC_2xCO2.append(CDNUMC_diff)
        CDNUMC_4xCO2.append(CDNUMC_Qdiff)
        CDNUMC_2xCO2_norm.append(CDNUMC_diff_norm)
        CDNUMC_4xCO2_norm.append(CDNUMC_Qdiff_norm)
        CLDLWP_sem.append(Cvar_sem)
        CLDLWP_dsem.append(Dvar_sem)
        CLDLWP_qsem.append(Qvar_sem)
        CLDLWP_2xCO2_sem.append(CLDLWP_diff_sem)
        CLDLWP_4xCO2_sem.append(CLDLWP_Qdiff_sem)
        CLDLWP_2xCO2_norm_sem.append(CLDLWP_diff_norm_sem)
        CLDLWP_4xCO2_norm_sem.append(CLDLWP_Qdiff_norm_sem)
        CDNUMC_sem.append(Cvar2_sem)
        CDNUMC_dsem.append(Dvar2_sem)
        CDNUMC_qsem.append(Qvar2_sem)
        CDNUMC_2xCO2_sem.append(CDNUMC_diff_sem)
        CDNUMC_4xCO2_sem.append(CDNUMC_Qdiff_sem)
        CDNUMC_2xCO2_norm_sem.append(CDNUMC_diff_norm_sem)
        CDNUMC_4xCO2_norm_sem.append(CDNUMC_Qdiff_norm_sem)

#            scre.append(SCRE)
#            lcre.append(LCRE)
#            delta_s.append(diffS)
#            delta_l.append(diffL)
        cn.append(float(casenames[i]))

        i+=1
    else:
        print('No such file or directory')


field = []
field.append(CLDLWP)
field.append(CDNUMC)
field.append(radius)
field_sem = []
field_sem.append(CLDLWP_sem)
field_sem.append(CDNUMC_sem)
field_sem.append(radius_sem)
dfield = []
dfield.append(CLDLWP_d)
dfield.append(CDNUMC_d)
dfield.append(radius_d)
dfield_sem = []
dfield_sem.append(CLDLWP_dsem)
dfield_sem.append(CDNUMC_dsem)
dfield_sem.append(radius_dsem)
qfield = []
qfield.append(CLDLWP_q)
qfield.append(CDNUMC_q)
qfield.append(radius_q)
qfield_sem = []
qfield_sem.append(CLDLWP_qsem)
qfield_sem.append(CDNUMC_qsem)
qfield_sem.append(radius_qsem)
Doubling_response = []
Doubling_response.append(CLDLWP_2xCO2)
Doubling_response.append(CDNUMC_2xCO2)
Doubling_response.append(radius_2xCO2)
Doubling_response_sem = []
Doubling_response_sem.append(CLDLWP_2xCO2_sem)
Doubling_response_sem.append(CDNUMC_2xCO2_sem)
Doubling_response_sem.append(radius_2xCO2_sem)
Quad_response = []
Quad_response.append(CLDLWP_4xCO2)
Quad_response.append(CDNUMC_4xCO2)
Quad_response.append(radius_4xCO2)
Quad_response_sem = []
Quad_response_sem.append(CLDLWP_4xCO2_sem)
Quad_response_sem.append(CDNUMC_4xCO2_sem)
Quad_response_sem.append(radius_4xCO2_sem)
normalized_Doubling_response = []
normalized_Doubling_response.append(CLDLWP_2xCO2_norm)
normalized_Doubling_response.append(CDNUMC_2xCO2_norm)
normalized_Doubling_response.append(radius_2xCO2_norm)
normalized_Doubling_response_sem = []
normalized_Doubling_response_sem.append(CLDLWP_2xCO2_norm_sem)
normalized_Doubling_response_sem.append(CDNUMC_2xCO2_norm_sem)
normalized_Doubling_response_sem.append(radius_2xCO2_norm_sem)
normalized_Quad_response = []
normalized_Quad_response.append(CLDLWP_4xCO2_norm)
normalized_Quad_response.append(CDNUMC_4xCO2_norm)
normalized_Quad_response.append(radius_4xCO2_norm)
normalized_Quad_response_sem = []
normalized_Quad_response_sem.append(CLDLWP_4xCO2_norm_sem)
normalized_Quad_response_sem.append(CDNUMC_4xCO2_norm_sem)
normalized_Quad_response_sem.append(radius_4xCO2_norm_sem)

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
        ax.set_title('Column Droplet Number',y=1.1)
        ax.set_yticks(np.linspace(6e9,1.1e10,6))
        ax.set_ylabel('$N_c$ ($1/m^2$)')
    if i == 2:
        ax.set_title('Average Cloud Droplet Radius')
        ax.set_yticks(np.linspace(9,12,4))
        ax.set_ylabel('radius ($\mu$m)')
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
        ax2.set_yticks(np.linspace(-3e9,0,4))
        ax2.set_ylabel('$N_c$ ($1/m^2$)')
    if i == 2:
        ax2.set_title('Difference')
        ax2.set_ylabel('radius ($\mu$m)')
        ax2.set_yticks(np.linspace(0.4,1.6,7))
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
        ax3.set_yticks(np.linspace(-2.5e8,0,6))
        ax3.set_title('Difference per unit RF')
        ax3.set_ylabel('$N_c$ ($1/m^2 / W/m^2$)')
        ax3.set_xlabel('Solar Constant Multiplier')
    if i == 2:
        ax3.set_title('Difference per unit RF')
        ax3.set_yticks(np.linspace(0,0.16,5))
        ax3.set_ylabel('radius ($\mu$m / $W/m^2$)')
        #ax3.set_ylim([-0.05,0.2])
        ax3.set_xlabel('Solar Constant Multiplier')
    #ax3.legend(['$2xCO_2$','$4xCO_2$'],ncol=2,prop={'size':6},loc='upper center')
    ax3.set_xticklabels(Casenames)
    ax3.grid(True)
    i+=1
###################################################################################################################################################
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles,labels,loc = (0.35, 0.86), ncol=5 )
fig.suptitle('Global Mean Cloud Droplet Radius Calculation and Comparison',y=1.02,fontsize=20)
fig.tight_layout()
plt.show()
fig.savefig(figure_path+'Global_mean_Droplet_radius.pdf',bbox_inches='tight')


