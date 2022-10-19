#!/usr/bin/env python3

import os
import sys
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.colors as colors
figure_path = '/home/brandonsmith/climate-gcm-bps/plots/'
casenames = ['0.9','0.95','1.0','1.05']

field = 'TGCLDLWP'
field2 = 'CDNUMC'
Zonal_Mean = False
Global_Mean = True
stf = False
if Global_Mean is True:
    for CASENAME in casenames:
        run = 'Control'
        inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        infile = '/Merged_'+run+'_'+CASENAME+'.nc'
        outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        outfile_control = '/cloud_droplets_GM_'+run+'_'+CASENAME+'.nc'
        if not os.path.isfile(outpath+outfile_control):
            if os.path.isfile(inpath+infile):
                syscall = 'cdo timmean -fldmean -seltimestep,-360/-1 -select,name='+field+','+field2+',AREA '+inpath+infile +' '+ outpath+outfile_control
                os.system(syscall)
        run = '2xCO2'
        inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        infile = '/Merged_'+run+'_'+CASENAME+'_.nc'
        outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        outfile_2xCO2 = '/cloud_droplets_GM_'+run+'_'+CASENAME+'.nc'
        if not os.path.isfile(outpath+outfile_2xCO2):
            if os.path.isfile(inpath+infile):
                syscall = 'cdo timmean -fldmean -seltimestep,-360/-1 -select,name='+field+','+field2+',AREL,AREI,AREA '+inpath+infile +' '+ outpath+outfile_2xCO2
                os.system(syscall)
        run = '4xCO2'
        inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        infile = '/Merged_'+run+'_'+CASENAME+'_.nc'
        outpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME
        outfile_4xCO2 = '/cloud_droplets_GM_'+run+'_'+CASENAME+'.nc'
        if not os.path.isfile(outpath+outfile_4xCO2):
            if os.path.isfile(inpath+infile):
                syscall = 'cdo timmean -fldmean -seltimestep,-360/-1 -select,name='+field+','+field2+',AREL,AREI,AREA '+inpath+infile +' '+ outpath+outfile_4xCO2
                os.system(syscall)
    
    # Load variables and perform calculations from outfiles
    i = 0
    outfile_1C = '/home/brandonsmith/modeloutput/Control/1.0/cloud_droplets_GM_Control_1.0.nc'
    outfile_1D = '/home/brandonsmith/modeloutput/2xCO2/1.0/cloud_droplets_GM_2xCO2_1.0.nc'
    colors = ['b','g','k','r']
    control_radius = []
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
    cn = []
    rf = [5.53394852, 4.72283731, 3.88779531, 2.78457941]
    rfq = [11.4426504, 9.82042797, 8.15034396, 5.94391217]
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
            
            #diff = diff - (Dvar1 - Cvar1)
            control_radius.append(Droplet_radius_C)
            radius_d.append(Droplet_radius_D)
            radius_q.append(Droplet_radius_Q)
            radius_2xCO2.append(droplet_diff)
            radius_4xCO2.append(droplet_Qdiff)
            radius_2xCO2_norm.append(droplet_diff_norm)
            radius_4xCO2_norm.append(droplet_Qdiff_norm)
            
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
    field.append(control_radius)
    dfield = []
    dfield.append(CLDLWP_d)
    dfield.append(CDNUMC_d)
    dfield.append(radius_d)
    qfield = []
    qfield.append(CLDLWP_q)
    qfield.append(CDNUMC_q)
    qfield.append(radius_q)
    Doubling_response = []
    Doubling_response.append(CLDLWP_2xCO2)
    Doubling_response.append(CDNUMC_2xCO2)
    Doubling_response.append(radius_2xCO2)
    Quad_response = []
    Quad_response.append(CLDLWP_4xCO2)
    Quad_response.append(CDNUMC_4xCO2)
    Quad_response.append(radius_4xCO2)
    normalized_Doubling_response = []
    normalized_Doubling_response.append(CLDLWP_2xCO2_norm)
    normalized_Doubling_response.append(CDNUMC_2xCO2_norm)
    normalized_Doubling_response.append(radius_2xCO2_norm)
    normalized_Quad_response = []
    normalized_Quad_response.append(CLDLWP_4xCO2_norm)
    normalized_Quad_response.append(CDNUMC_4xCO2_norm)
    normalized_Quad_response.append(radius_4xCO2_norm)
    
    #plt.plot(cn,ts,'o',color='b')
    fig = plt.figure(figsize=(10,5))
    i = 0
    Casenames = ['90%','95%','100%','105%']
    while i < 3:
        ax = fig.add_subplot(3,3,i+1)
        #plt.plot(cn,ts,'o',color='b')
        ax.plot(cn,field[i],'.',color='k',label='Control',alpha=0.7)
        ax.plot(cn,dfield[i],'.',color='b',label='$2xCO_2$',alpha=0.7)
        ax.plot(cn,qfield[i],'.',color='r',label='$4xCO_2$',alpha=0.7)
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
        ax2.plot(cn,Doubling_response[i],'.',color='b',alpha=0.7)
        ax2.plot(cn,Quad_response[i],'.',color='r',alpha=0.7)
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
        ax3.plot(cn,normalized_Doubling_response[i],'.',color='b',alpha=0.7)
        ax3.plot(cn,normalized_Quad_response[i],'.',color='r',alpha=0.7)
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


