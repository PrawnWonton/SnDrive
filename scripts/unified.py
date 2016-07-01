# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 19:43:42 2015

@author: Prawn Wonton
"""

import numpy as np
import yt
import matplotlib.pyplot as plt
from yt.mods import *
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import AxesGrid, make_axes_locatable
import sys
import math
from path import path
#from pylab import *
import os
import yt.visualization.eps_writer as eps
import time

def get_power(fn,pf,dd,fout, year):
    #pf = yt.load('/media/marston/4EC85DD3C85DBA43/Users/rlandeen/Documents/SnDrive/fixed_cool/take5/resume/sndrive_fixcool5_hdf5_plt_cnt_0001')
    #fn = 1
    fin = fn
    M = int(pf.domain_dimensions[0])
    n = np.array([int(M),int(M),int(M)])
    k0 = 0
    kc = M/2



#5-26-2016 ntoes
#-------------------
#add density weight to power spectra sqrt(rho)*v
##double check box sizes
#add density weighted sound speed graph
#check iso/strat equations (driving scale, scale height, etc)
#check integral scale units
#power spectrum unit question

    reg = pf.h.covering_grid(0,left_edge=pf.domain_left_edge,dims=[M,M,M])
    vx = reg['velx']
    vy = reg['vely']
    vz = reg['velz']
    rho = reg['density']
    s = np.array([M,M,M])
    axes = np.array([0,1,2])
    vx = np.abs( np.fft.fftn(vx,s,axes) ) / float(M**3)
    vy = np.abs( np.fft.fftn(vy,s,axes) ) / float(M**3)
    vz = np.abs( np.fft.fftn(vz,s,axes) ) / float(M**3)

    kx = np.zeros(n)
    ky = np.zeros(n)
    kz = np.zeros(n)

    for j in range(0,n[1]):
        for k in range (0,n[2]):
            kx[:,j,k] = n[0]*np.fft.fftfreq(n[0])
            #time.sleep(0.100)
            #print('kx = {:}'.format(k))
    for i in range(0,n[0]):
        for k in range (0,n[2]):
            ky[i,:,k] = n[1]*np.fft.fftfreq(n[1])
            #time.sleep(0.100)
            #print('ky = {:}'.format(k))
    for i in range(0,n[0]):
        for j in range (0,n[1]):
            kz[i,j,:] = n[2]*np.fft.fftfreq(n[2])
            #time.sleep(0.100)
            #print('kz = {:}'.format(k))

    k = np.sqrt( (kx**2+ky**2+kz**2) )

    k1d=[]
    power=[]
    for i in range(k0,kc+1):
        si = np.where( np.logical_and( k>=float(i),k<float(i+1) ) )
        k1d.append(i)
        power.append( np.sum( ((vx[si]**2) + (vy[si]**2) + (vz[si]**2))*(rho[si]) ) )#stick rho in here
        #time.sleep(0.100)
        print('i = {:}'.format(i))
    k1d=np.array(k1d)
    power=np.array(power)

    np.save("plt_"+str(fn)+"_k1d", k1d)
    np.save("plt_"+str(fn)+"_power", power)
    #fig = plt.figure()

    plt.clf()
    plt.loglog(k1d,power,'b',lw=3)
    plt.xlabel(r'k/2$\pi$', fontsize = 15)
    plt.ylabel('E(k)', fontsize = 15)
    plt.xlim([1,M])
    plt.title('Power Spectrum at t = '+str(year)+' Myrs', fontsize = 20)
    plt.gca().xaxis.grid(True)
    #ln = 'power_%s.pdf'%(str(fn[-4:]))
    print('\n')
    print('-----')
    print('writing plot file: {}'.format(fn))
    print('-----')
    plt.savefig('Power'+str(fn)+'.png')

    slopey = k1d[11]-k1d[81]
    print("Slope from k/2pi = 10 to 80: "+str(slopey))

    #----------------------------
    #---Integral Scale Gubbins---
    #----------------------------
    b = float(k1d[-1])
    a = float(k1d[0])
    n = len(k1d)
    h = float((b-a)/n)

    constant = (3*math.pi)/4

    E = float(0)
    E_k = float(0)
    index = 1
    while (index < n):
        E_k = E_k + (power[index]+power[index-1])/(k1d[index]+k1d[index-1])*(k1d[index]-k1d[index-1])
        E = E + ((power[index]+power[index-1])/2.)*(k1d[index]-k1d[index-1])
        index = index + 1
    E_k = float(h*E_k)
    E = float(h*E)
    intscale = float(constant*(E_k/E))
    print("---integral scale---")
    print(intscale)
    print("--------------------")

###-----------------------------------###
###---output integral scale to file---###
###-----------------------------------###

    f = open("integralscale_v2.txt", 'a')
    f.write("t = "+str(year)+" Myrs\n")
    f.write("integral scale: "+str(intscale)+"\n")
    f.close()

    #----------------------------

    #fig = plt.figure()
    #k_slope = k1d[10:81:1]
    #power_slope = np.diff(power)
    #power_slope = power_slope[10:81:1]
    #plt.plot(k_slope, power_slope, 'b', lw=3)
    #plt.xlabel(r'k/2$\pi$')
    #plt.ylabel('Slope of E(k)')
    #plt.savefig('Slope'+str(fn)+'.png')


#---END datplot def---
#---------------------

def datplot(location, datfile):

    file_in = np.loadtxt(location+datfile, delimiter=None, skiprows=1)
    run_number = int(datfile[4:5])
    time = file_in[::,0]
    mass = file_in[::,1]
    x_momentum = file_in[::,2]
    y_momentum = file_in[::,3]
    z_momentum = file_in[::,4]
    E_total = file_in[::,5]
    E_kinetic = file_in[::,6]
    E_internal = file_in[::,7]
    rms_vel_vavg = file_in[::,8]
    rms_vel_mavg = file_in[::,9]
    rms_sndspd_dw = file_in[::,10]
    rms_mach_mw = file_in[::,11]
    mtemp_vavg = file_in[::,12]
    mtemp_mavg = file_in[::,13]
    Mag_vorticity = file_in[::,14]
    min_dens = file_in[::,15]
    max_dens = file_in[::,16]

    sqrt_rms_v = np.sqrt(file_in[::,8])/1000000
    sqrt_rms_m = np.sqrt(file_in[::,9])/1000000

    ###-------------------------###
    ###---Velocity Dispersion---###
    ###-------------------------###

    vel_disp_KE = (np.sqrt((2*E_kinetic)/mass))*1e-5
    #vel_disp_Etot = (np.sqrt((2*E_total)/mass))*1e-5
    #vel_disp_Eint = (np.sqrt((2*E_internal)/mass))*1e-5

    mid = len(E_kinetic)/2

    avg_vel_disp_KE = np.mean(vel_disp_KE[mid:])
    #avg_vel_disp_Etot = np.mean(vel_disp_Etot[mid:])
    #avg_vel_disp_Eint = np.mean(vel_disp_Eint[mid:])

    print('---')
    print(avg_vel_disp_KE)
    #print(avg_vel_disp_Etot)
    #print(avg_vel_disp_Eint)
    print('---')

    plot_list = [mass, E_total, E_kinetic, E_internal, rms_vel_vavg, rms_vel_mavg, rms_sndspd_dw, rms_mach_mw, mtemp_vavg, mtemp_mavg, Mag_vorticity, min_dens, max_dens, sqrt_rms_v, sqrt_rms_m, vel_disp_KE]
    plot_labels_file = ["mass", "E_total", "E_kinetic", "E_internal", "rms_vel_vavg", "rms_vel_mavg", "rms_sndspd_dw", "rms_mach_mw", "mtemp_vavg", "mtemp_mavg", "Mag_vorticity", "min_dens", "max_dens", "sqrt rms v", "sqrt rms m", "vel disp km per s (KE)"]
    plot_labels_title = ["Mass", "Energy (total)", "Energy (kinetic)", "Energy (internal)", "rms_vel_vavg", "rms_vel_mavg", "Density Weighted Sound Speed", "rms_mach_mw", "mtemp_vavg", "mtemp_mavg", "Mag_vorticity", "min_dens", "max_dens", "sqrt rms v", "sqrt rms m", "Velocity Dispersion"]
    plot_labels_units = ["Mass", "Energy (total)", "Energy (kinetic)", "Energy (internal)", "rms_vel_vavg", "rms_vel_mavg", "Density Weighted Sound Speed", "rms_mach_mw", "mtemp_vavg", "mtemp_mavg", "Mag_vorticity", "min_dens", "max_dens", "sqrt rms v", "sqrt rms m", "Velocity Dispersion (km/s)"]
    log_list = np.log(plot_list)

    print("hrmpf")

    ##---------------###
    ## Plotting Time ###
    ##---------------###

    x = time
    start = 1000
    i=0 #initialize counter
    while (i < len(plot_list) ):

        plt.plot((x)[start:], (plot_list[i][start:]))
        plt.ylabel(plot_labels_units[i], fontsize = 16)
        plt.xlabel("Time (s)", fontsize = 16)
        plt.title("Run {} ".format(run_number)+plot_labels_title[i]+" vs time", fontsize=20)
        plt.gca().yaxis.grid(which='minor', alpha=0.2)
        plt.gca().xaxis.grid(True, linestyle='-', alpha=0.3)
        plt.gca().yaxis.grid(True, linestyle='-', alpha=0.3)
        plt.savefig("plt_"+plot_labels_file[i]+".png")
        plt.clf()

        plt.loglog((x)[start:], (plot_list[i][start:]))
        plt.ylabel(plot_labels_units[i], fontsize = 16)
        plt.xlabel("Time in s", fontsize = 16)
        plt.title("Run {} ".format(run_number)+plot_labels_title[i]+" vs time")
        plt.gca().yaxis.grid(which='minor', alpha=0.2)
        plt.gca().xaxis.grid(True, linestyle='-', alpha=0.3)
        plt.gca().yaxis.grid(True, linestyle='-', alpha=0.3)
        plt.savefig("plt_log_"+plot_labels_file[i]+".png")
        plt.clf()

    #    plt.plot((x), (plot_list[i]))
    #    plt.ylabel(plot_labels[i])
    #    plt.xlabel("time")
    #    plt.title(plot_labels[i]+" vs time")
    #    #plt.annotate("hrmpf", (1.26E10,0), "axes points")
    #    plt.savefig(plot_labels[i]+".png")
    #    plt.clf()

        i= i+1
#---END datplot def---
#---------------------

#--- File Plotter ---
#--------------------

def file_plotter(file_name, location, logfile):
    inF = open(location+logfile, 'r')
    for line in inF:
        if 'x_sn' in line:
            for word in line.split():
                xvalue = word
            for word in inF.next().split():
                yvalue = word
            for word in inF.next().split():
                zvalue = word
            break
    center = [float(xvalue), float(yvalue), float(zvalue)]
    #--- set values to scale axes to for all slices so they are comparable ---
    first_plt = yt.load(location+str(file_name))
    sp = first_plt.sphere(center, (10., "pc"))
    # density min/max
    dens_zmin = float(sp["density"].min()/10)
    dens_zmax = float(sp["density"].max()*1000)
    # temp min/max
    temp_zmin = float(sp["temperature"].min()/1.2)
    temp_zmax = float(sp["temperature"].max()*1.2)
    # pres min/max
    pres_zmin = float(sp["pressure"].min()/1)
    pres_zmax = float(sp["pressure"].max()*10000)
    # velx min/max
    velx_zmin = float(sp["velx"].min()/1)
    velx_zmax = float(sp["velx"].max()*1)

    for root, dirs, files in os.walk(location):
        for file in files:
            if file.startswith(file_name[:-3]):
                i=int(file[-3:])
                ds = yt.load(location+str(file))
                width_pc = int(ds.domain_width.in_units('pc')[0])
                year = round(ds.current_time.in_units('yr'), 5)

                #---density
                slc_dens = yt.SlicePlot(ds, 'z', 'density', width = (width_pc, 'pc'), center=center)
                slc_dens.set_log('density',True)
                #slc_dens.set_zlim('density', dens_zmin, dens_zmax)
                slc_dens.annotate_title('Density at  at t={:} years'.format(year))
                #slc_dens.annotate_title('Density')
                #slc_dens.annotate_timestamp()
                #eps_fig_dens = eps.single_plot(slc_dens)
                #eps_fig_dens.title_box('Density at t = '+str(year)+' yrs')

                #---temp
                slc_temp = yt.SlicePlot(ds, 'z', 'temperature', width = (width_pc, 'pc'), center=center)
                slc_temp.set_log('temperature',False)
                #slc_temp.set_zlim('temperature', temp_zmin, temp_zmax)
                #slc_temp.annotate_title('Temperature at t = '+str(year)+' years')
                slc_temp.annotate_title('Temperature at t={:} years'.format(year))
                #slc_temp.annotate_timestamp()
                #eps_fig_temp = eps.single_plot(slc_temp)
                #eps_fig_temp.title_box('Temp at t = '+str(year)+' yrs')

                #---oden
                #slc_oden = yt.SlicePlot(ds, 'z', 'oden', center=[1.119021E+20, 4.634398E+20, 3.921526E+20])
                #slc_oden.set_log('oden',False)
                #slc_dens.annotate_title('Density at t = '+str(year)+' years')
                #eps_fig_oden = eps.single_plot(slc_oden)
                #eps_fig_oden.title_box('oden at t = '+str(year)+' yrs')

                #---velx
                slc_velx = yt.SlicePlot(ds, 'z', 'velx', width = (width_pc, 'pc'), center=center)
                slc_velx.set_log('velx',False)
                #slc_velx.set_zlim('velx', velx_zmin, velx_zmax)
                slc_velx.annotate_title('Velx at  at t={:} years'.format(year))
                #slc_velx.annotate_title('Velocity (x)')
                #slc_velx.annotate_timestamp()
                #eps_fig_velx = eps.single_plot(slc_velx)
                #eps_fig_velx.title_box('Velx at t = '+str(year)+' yrs')

                #---vely
                #slc_vely = yt.SlicePlot(ds, 'z', 'vely', width = (width_pc, 'pc'), center=center)
                #slc_vely.set_log('vely',False)
                #slc_vely.annotate_title('Vely at t = '+str(year)+' years')
                #slc_vely.annotate_title('Velocity (y)')
                #slc_vely.annotate_timestamp()
                #eps_fig_vely = eps.single_plot(slc_vely)
                #eps_fig_vely.title_box('Vely at t = '+str(year)+' yrs')

                #--- output

                #eps_fig_dens.save_fig('Density_'+str(i), format='pdf')
                #eps_fig_velx.save_fig('Velx_'+str(i), format='pdf')
                #eps_fig_vely.save_fig('Vely_'+str(i), format='pdf')
                slc_dens.save('Density_{:}'.format(i))
                slc_velx.save('Velx_{:}'.format(i))
                #slc_vely.save('Vely_{:}'.format(i))
                slc_temp.save('Temp_{:}'.format(i))
                #eps_fig_temp.save_fig('Temp_'+str(i), format='pdf')
                #eps_fig_oden.save_fig('Oden_'+str(i), format='pdf')
                #slc_dens.save('density'+str(i))
                #i=i+1


### Global Variables ###
###------------------###

location = "F:\\SnDrive\\proposal runs\\384_test\\run_1\\"
chkpoint = "run_1_cube48_hdf5_plt_cnt_0001"
plt_file = chkpoint
datfile = "run_1_cube48.dat"
logfile = "run_1_cube48.log"

#--- Set working dir ---
#-----------------------
scriptdir = 'C:\\Users\\Prawn Wonton\\Dropbox\\scripts'
os.chdir(location+'plots2')

#--- file plotter ---
#file_plotter(plt_file, location, logfile)

#---

#--- dat plotter ---
datplot(location, datfile)

#--- power spectrum ---
'''
p = path(location)
for f in p.files(pattern=chkpoint[:-2]+str('*')):
    dataset = yt.load(f)
    year = round(dataset.current_time.in_units('Myr'), 5)
    j=int(f[-2:])
    if(j>0):
        if(j<20):
            get_power(j,dataset,0,0, year)
'''

#change back to script directory
os.chdir(scriptdir)
#--- npy to dat converter ---
#p = path(location)
#for f in p.files(pattern='k1d'):
#    k1d = np.load(f)
#    np.savetxt(f+str('.dat'), k1d)
#for f in p.files(pattern='*power.npy'):
#    power2 = np.load(f)
#    np.savetxt(f+str('.dat'), power2)
