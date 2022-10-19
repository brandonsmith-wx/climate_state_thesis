#!/usr/bin/env python3

#
# Script for plotting figure 3.3. for altercation, make necessary comments and uncomments
# of particular axes settings in the plotting section.
#

import os
import sys
import numpy as np
import netCDF4
import matplotlib.pyplot as plt


field = 'TS'
SOLIN = 'SOLIN'
co2vmr = 'co2vmr'
FSNT = 'FSNT'
FLNT = 'FLNT'
filebase = '/home/brandonsmith/modeloutput/'
outfilebase = 'ts_co2_'
casenames = ['0.9','0.95','1.0','1.05']
figure_path = '/home/brandonsmith/climate-gcm-bps/plots/'

run = 'Control'
# calc temp control
for CASENAME in casenames:
    outfile = filebase+run+'/'+ CASENAME+'/'+outfilebase+CASENAME+'.nc'
    # check if the outfile exists
    if not os.path.isfile(outfile):
        if os.path.isdir(filebase+run+'/'+CASENAME):
            infile = filebase + run+'/'+CASENAME+'/Merged_'+run+'_'+CASENAME+'.nc'
            infile_control = '/home/brandonsmith/modeloutput/Control/1.0/Merged_Control_1.0.nc'
            syscall = 'cdo fldmean -timmean -seltimestep,-360/-1 -select,name='+field+','+SOLIN+','+co2vmr+','+FSNT+','+FLNT+' '+infile+' '+outfile
            print(syscall)


run = '2xCO2'
for CASENAME in casenames:
    outfile = filebase+run+'/'+ CASENAME+'/'+outfilebase+CASENAME+'.nc'
    # check if the outfile exists
    if not os.path.isfile(outfile):
        if os.path.isdir(filebase+run+'/'+CASENAME):
            infile = filebase + run+'/'+CASENAME+'/Merged_'+run+'_'+CASENAME+'_.nc'
            # merge files and calc TS global average per month
            syscall = 'cdo fldmean -timmean -seltimestep,-360/-1 -select,name='+field+','+SOLIN+','+co2vmr+','+FSNT+','+FLNT+' '+infile+' '+outfile
            print(syscall)
run = '4xCO2'
for CASENAME in casenames:
    outfile = filebase+run+'/'+ CASENAME+'/'+outfilebase+CASENAME+'.nc'
    # check if the outfile exists
    if not os.path.isfile(outfile):
        if os.path.isdir(filebase+run+'/'+CASENAME):
            infile = filebase + run+'/'+CASENAME+'/Merged_'+run+'_'+CASENAME+'_.nc'
            # merge files and calc TS global average per month
            syscall = 'cdo fldmean -timmean -seltimestep,-360/-1 -select,name='+field+','+SOLIN+','+co2vmr+','+FSNT+','+FLNT+' '+infile+' '+outfile
            print(syscall)
    

#create plot
C0 = 0.000278
a = 0.39
b = 5.32
sens = 0.5 #K/w/m^2


cn = []
control = []
doubled = []
quadrupled = []
co2 = []
co2_2x = []
co2_4x = []
Fbg = []
F_net = []
i=0
#plot the data
for CASENAME in casenames:
    dsloc = filebase+'Control'+'/'+CASENAME+'/'+outfilebase+CASENAME+'.nc'
    dsloc_2x = filebase+'2xCO2'+'/'+CASENAME+'/'+outfilebase+CASENAME+'.nc'
    dsloc_4x = filebase+'4xCO2'+'/'+CASENAME+'/'+outfilebase+CASENAME+'.nc'
    if os.path.isfile(dsloc_4x) and os.path.isfile(dsloc_2x):

        # open the merged file and get out the variables
        ds = netCDF4.Dataset(dsloc)
        ds2 = netCDF4.Dataset(dsloc_2x)
        cs = ds.variables[field][:]
        cs2 = ds2.variables[field][:]
        co2vmr = ds.variables['co2vmr'][:]
        co2vmr_2x = ds2.variables['co2vmr'][:]
        ds4 = netCDF4.Dataset(dsloc_4x)
        cs4 = ds4.variables[field][:]
        co2vmr_4x = ds4.variables['co2vmr'][:]
        ds4.close()
        ds.close() #close the file
        ds2.close() #close the file
        
        cs = cs.flatten()
        cs2 = cs2.flatten()
        cs4 = cs4.flatten()
        control = np.append(control,cs)
        doubled = np.append(doubled,cs2)
        quadrupled = np.append(quadrupled,cs4)
        co2 = np.append(co2,co2vmr)
        co2_2x = np.append(co2_2x,co2vmr_2x)
        co2_4x = np.append(co2_4x,co2vmr_4x)
        cn = np.append(cn,float(casenames[i]))
        i=i+1
    else:
        print('could not validate data set location')
        print(casenames[i])
        i=i+1

control=np.array(control)
doubled=np.array(doubled)
quadrupled = np.array(quadrupled)
co2 = np.array(co2)
co2_2x = np.array(co2_2x)
co2_4x = np.array(co2_4x)


diff = doubled - control
diff_4 = quadrupled - control

diff = np.array(diff)
F0 = a*(np.log(co2/C0))**2+b*(np.log(co2/C0))
F0_2 = a*(np.log(co2_2x/C0))**2+b*(np.log(co2_2x/C0))
F0_4 = a*(np.log(co2_4x/C0))**2+b*(np.log(co2_4x/C0))
dF = diff/sens
Fn = F0_2 - F0 # Added Radiative Forcing per CO2 doubling
Fn_4 = F0_4 - F0 # Added Radiative Forcing per CO2 quadrupling

T2_per_F = diff/Fn
T4_per_F = diff_4/Fn_4
sol_str = ['90%','95%','100%','105%']
#plot the data
cn2 = [0.95,1.0,1.05]

fig = plt.figure()
ax1 = fig.add_subplot(3,1,1)
ax1.plot(cn, control, '.', color='k', label='Control', alpha=0.7)
ax1.plot(cn, doubled, '.', color='b', label='$2xCO_2$', alpha=0.7)
ax1.plot(cn, quadrupled, '.', color='r', label='$4xCO2$', alpha=0.7)
#ax1.hlines(0,0.8,1.1)
#ax1.plot(co2, doubled, marker='o', color='blue', alpha=0.7, label='2xCO2')
ax1.set_title('Global Mean Surface Temperature',fontsize=12)
ax1.set_ylabel('Temperature (K)',fontsize=8)
ax1.set_xticks([0.9,0.95,1.0,1.05])
ax1.set_xticklabels(sol_str)
#ax1.set_ylim(-0.3,0.3)
ax1.set_xlim(0.89,1.06)
ax1.set_yticks(np.linspace(286,302,5))
ax1.set_yticklabels(np.linspace(286,302,5).astype(int).astype(str),fontsize=8)
ax1.grid(True)
#ax1.legend(['$2xCO_2$','$4xCO_2$'],loc='upper right',ncol=2,prop={'size':5})
'''
ax2 = plt.subplot(312)
ax2.plot(cn, CO2bg, marker='.',color='blue',alpha=0.7)
ax2.plot(cn, Next_co2, marker='.',color='green',alpha=0.7)
ax2.set_ylabel(' CO2 volume mixing ratio (ppmv)',fontsize=8)
ax2.set_yscale('log')
ax2.set_xticks([0.9,0.95,1.0,1.05])
ax2.set_xticklabels(casenames)
ax2.set_yticks([10,100,1000,10000])
ax2.set_yticklabels(['10¹','10²','10³','10⁴'])
ax2.set_xlim(0.89,1.06)
ax2.set_ylim(0,40000)
ax2.grid(True)
ax2.legend(['Byrne Goldblatt 2014','next guess'])


ax3 = plt.subplot(313)
ax3.plot(cn, ratio, marker='.', color='black',alpha=0.7)
ax3.hlines(0,0.8,1.1)
ax3.set_ylabel('Change in CO2 (ppmv)',fontsize=8)
ax3.set_xticks([0.9,0.95,1.0,1.05])
ax3.set_xticklabels(casenames)
#ax3.set_yticks([-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30])
ax3.set_xlim(0.89,1.06)
#ax3.set_ylim(-35,35)
#ax3.set_yticks(numpy.linspace(-0.04,0.04,9))
ax3.grid(True)

ax3 = plt.subplot(312)
ax3.plot(cn2, F0_4[1:], marker='.',color='magenta',alpha=0.7)
ax3.plot(cn, F0_2, marker='.', color='red',alpha=0.7)
ax3.plot(cn, F0, marker='.', color='blue',alpha=0.7)
#ax3.plot(cn, F0, marker='.', color='black',alpha=0.7)
#ax3.hlines(0,0.8,1.1)
ax3.set_title('Imposed Radiative Forcing from $CO_2$ Doubling (Byrne & Goldblatt 2014)',fontsize=10)
ax3.set_ylabel('Flux (W/m²)',fontsize=8)
ax3.set_xticks([0.9,0.95,1.0,1.05])
ax3.set_xticklabels(co2_str)
ax3.set_yticks(np.linspace(-10,25,8))
ax3.set_yticklabels(np.linspace(-10,25,8).astype(int).astype(str),fontsize=7)
#ax3.set_title('Net Longwave Flux at TOM',fontsize=10)
ax3.set_xlim(0.89,1.06)
ax3.grid(True)
ax3.legend(['$4xCO_2$','$2xCO_2$','Base'],loc='lower left',prop={'size':5})
'''
ax4 = fig.add_subplot(312)
ax4.plot(cn, diff, '.', color='blue',alpha=0.7,label='$2xCO_2$')
ax4.plot(cn,diff_4, '.',color='red', alpha=0.7,label='$4xCO_2$')
#ax4.hlines(0,0.8,1.1)
ax4.set_ylabel('Temperature (K)',fontsize=8)
ax4.set_title('Surface Temperature Response to $CO_2$ Doubling ',fontsize=12)
ax4.set_xticks([0.9,0.95,1.0,1.05])
ax4.set_xticklabels(sol_str)
ax4.set_yticks(np.linspace(0,14,8))
ax4.set_yticklabels(np.linspace(0,14,8).astype(int).astype(str),fontsize=7)
#ax4.set_title('Difference',fontsize=10)
ax4.set_xlim(0.89,1.06)
ax4.grid(True)
#ax4.legend(['$2xCO_2$','$4xCO_2$'],loc='upper right',ncol=2,prop={'size':5})

ax3 = fig.add_subplot(313)
ax3.plot(cn, T2_per_F, '.', color='blue',alpha=0.5,label='$2xCO_2$')
ax3.plot(cn,T4_per_F, '.',color='red', alpha=0.5,label='$4xCO_2$')
ax3.set_ylabel('Difference (K/W/m²)',fontsize=8)
ax3.set_title('Change in Surface Temperature per Radiative Forcing unit ($K/W/m^2$)',fontsize=12)
ax3.set_xticks([0.9,0.95,1.0,1.05])
ax3.set_xticklabels(sol_str)
ax3.set_yticks(np.linspace(0.8,1.2,5))
ax3.set_yticklabels(np.linspace(0.8,1.2,5).astype(str),fontsize=8)
#ax4.set_title('Difference',fontsize=10)
ax3.set_xlim(0.89,1.06)
ax3.grid(True)
#ax3.legend(['$2xCO_2$','$4xCO_2$'],loc='upper right',ncol=2,prop={'size':5})
handles, labels = ax1.get_legend_handles_labels()
fig.legend(handles,labels,loc = (0.25, 0.83), ncol=3)
plt.suptitle('Climate Sensitivity Response to $CO_2$ Increases',fontsize=16,y=1.1)
plt.xlabel('Fraction of Present Day Solar Constant')
fig.tight_layout(pad=0.5)

#plt.axis([0.89, 1.04, -52,-40])
plt.show()



fig2 = plt.figure()

ax5 = plt.subplot(1,1,1)
ax5.plot(cn, Fn, 'o', color='black',alpha=0.7)
ax5.plot(cn,Fn_4, 'o',color='red', alpha=0.7)
#ax4.hlines(0,0.8,1.1)
ax5.set_ylabel('Radiative Forcing (W/m²)',fontsize=12)
ax5.set_title('Added Radiative Forcing from Additional $CO_2$ ',fontsize=12)
ax5.set_xticks([0.9,0.95,1.0,1.05])
ax5.set_xticklabels(sol_str,fontsize=12)
ax5.set_yticks(np.linspace(0,14,8))
ax5.set_yticklabels(np.linspace(0,14,8).astype(int).astype(str),fontsize=12)
#ax4.set_title('Difference',fontsize=10)
ax5.set_xlim([0.89,1.06])
ax5.set_xlabel('Fraction of Present Solar Constant',fontsize=12)
ax5.grid(True)
ax5.legend(['$2xCO_2$','$4xCO_2$'],loc='upper right',ncol=2,prop={'size':10})
plt.show()

fig.savefig(figure_path+"GMST_CS_Response.pdf",bbox_inches='tight')
fig2.savefig(figure_path+'Radiative_Forcing_Comparison.pdf')