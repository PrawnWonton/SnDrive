# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 14:41:01 2016

@author: prawnwonton
"""
#--- Import Gubbins ---#
#----------------------#

from fractions import Fraction
import numpy as np

#--- End Importation of Gubbins ---#
#----------------------------------#

#--- Variables ---#
#-----------------#
sigma_g = [250, 250, 150, 150, 150, 50, 50]
sigma_sn = [0.002, 0.001, 0.001, 0.0006, 0.0003, 0.0003, 0.0001]
sigma_sn = [x*10**6 for x in sigma_sn]

#Velocity disperion for c_s = 15 km/s
#vel_disp = [4.43, 3.28, 5.08, 4.07, 3.01, 7.77, 4.79] # (x^2)/3  ## old incorrect equation for density
#new and improved equations (6-6-2016)
vel_disp = [5.118, 3.785, 5.874, 4.698, 3.476, 9.011, 5.536]


#Velocity disperion for c_s = 10 km/s
#vel_disp = [3.751, 2.772, 4.308, 3.443, 2.546, 6.626, 4.059]

cs = 15
fp = 1
fg = .5
resolution = 256


#--- Output Formatting ---#
#-------------------------#

new_n0_array   = np.array([])
ndot_array     = np.array([])
R_st_array     = np.array([])
L_array        = np.array([])
LtoR_ratio     = np.array([])
rho_array      = np.array([])
box_array      = np.array([])
dt_sn_array    = np.array([])
radius_sn_array= np.array([])
velocity_array = np.array([])
dt_sn_reset_array = np.array([])
scale_height_array = np.array([])
L_strat_array  = np.array([])
H_array        = np.array([])
n0_11_array    = np.array([])
new_rho_array = np.array([]) #sigma_g / 2H

i=0
while (i<7):

    #--- Equations to crunch ---#
    #---------------------------#

    #corrected pressure fraction using corrected velocity dispersion
    f_p = 1/(1+(3*cs**2)/(vel_disp[i]**2))

    #number density [cm^-3]
    new_n0 = (9.0*10**(-4))*((f_p/fg)**(1.72))*((sigma_g[i])**4.92)*(((sigma_sn[i])**(-1.48)))

    #snupernova rate [10^-12 yr pc^-3]
    #n_dot = .013*(f_p/fg)*sigma_sn[i]*sigma_g[i]*vel_disp[i]**(-2)
    n_dot = (1.3*10**(-5))*((f_p/fg)**(1.72))*((sigma_g[i])**3.92)*((sigma_sn[i])**(-.48))

    #Cooling radius
    R_st = 19.1*((new_n0)**(Fraction(-7,17)))

    #Driving radius isotropic/stratified
    L = 45*(new_n0**Fraction(-19,119))*(n_dot**Fraction(-1,7))
    L_strat = 1440*((f_p/fg)**(-.52))*(sigma_g[i]**(-1.34))*(sigma_sn[i]**(.3))

    #scale height
    scale_height = ((3.9*10**4)*((f_p/fg)**-1.72)*(sigma_g[i]**-3.92)*(sigma_sn[i]**1.48))

    #--- Flash.par parameters ---#
    rho = new_n0*.6*(1.6726e-24)
    #Injection radius of snupernova
    R_inj = .5*R_st
    #box size
    box = R_st/10*resolution
    #supernova rate
    dt_sn = ((1e12)/n_dot)*(3.154e7)/((box)**3)
    #reset time between snupernovas
    dt_sn_reset = (box*(3.08e18))/(3e10)

    #eq. 10
    H = 74*(vel_disp[i]**2)*(sigma_g[i])**(-1)*(fg/f_p)
    #eq. 11
    n0_11 = 0.48*(f_p/fg)*(sigma_g[i]**2)*(vel_disp[i]**(-2))

    new_rho = sigma_g[i]/(2*H)

    print('---')
    print('resolution [pc/cell]: {} at {} resolution'.format(box/resolution, resolution))
    print('for sigma_g = '+str(sigma_g[i])+' M_sol pc^-2'+' and sigma_sn = '+str(sigma_sn[i])+' Myr^-1 kpc^2')
    print('n0 = {}'.format(new_n0)+' cm^-3')
    print('n_dot = '+str(n_dot)+' 10^-12 year^-1 pc^-3')
    print('Cooling radius R_st = {:.6}'.format(R_st)+' pc')
    print('Injection radius = {:.6}'.format(R_inj)+' pc' + ' [{:.5} cm]'.format(R_inj*(3.08e18)))
    print('Driving scale L (iso) = {:.6}'.format(L)+' pc')
    print('Driving scale L_strat = {:.6}'.format(L_strat)+' pc')
    print('Velocity dispersion = {:.6}'.format(vel_disp[i])+' km/s')
    print('Ratio of L/R : {:.4}'.format(L_strat/R_st))

    print('flash.par parameters:')
    print('rho_ambient = n0*.6*m_p = {:.3e}'.format(rho)+' g/cm^3')
    print('Box size (R_st)/10*{}: {:.4}'.format(resolution, box)+' pc = {:.4e}'.format(box*3.08e18)+' cm')
    print('Midpoint (box/2): {:.5} cm'.format(box*3.08e18/2))
    print('dt_sn = {:.4e}'.format(dt_sn)+' seconds. ({:.2f} years)'.format(dt_sn/(3.154e7)))
    print('dt_sn_reset = box/c = {:.4e}'.format(dt_sn_reset)+' seconds')
    print('radius_sn = R_inj = {:.6e}'.format(R_inj*(3.08e18))+' cm')
    print('Scale height = {:}'.format(scale_height)+' pc')

    print('H = {:.6e}'.format(H))
    print('eq. 11 n_0 = {:.6e}'.format(n0_11))
    print('sigma_g / 2H = {:.6e}'.format(new_rho))

    #--- Output
    new_n0_array = np.append(new_n0_array, (new_n0))
    ndot_array = np.append(ndot_array, (n_dot))
    R_st_array = np.append(R_st_array, (R_st))
    L_array = np.append(L_array, (L))
    LtoR_ratio = np.append(LtoR_ratio, (L_strat/R_st))
    rho_array = np.append(rho_array, (rho))
    box_array = np.append(box_array, (box))
    dt_sn_array = np.append(dt_sn_array, (dt_sn))
    radius_sn_array = np.append(radius_sn_array, (R_inj))
    velocity_array = np.append(velocity_array, (vel_disp[i]))
    dt_sn_reset_array = np.append(dt_sn_reset_array, (dt_sn_reset))
    scale_height_array = np.append(scale_height_array, (scale_height))
    L_strat_array = np.append(L_strat_array, (L_strat))
    H_array = np.append(H_array, (H))
    n0_11_array = np.append(n0_11_array, (n0_11))
    new_rho_array = np.append(new_rho_array, (new_rho))

    #new_n0_array = np.append(new_n0_array, str(new_n0))
    #ndot_array = np.append(ndot_array, str(n_dot))
    #R_st_array = np.append(R_st_array, str(R_st))
    #L_array = np.append(L_array, str(L_strat))
    #LtoR_ratio = np.append(LtoR_ratio, str(L_strat/R_st))
    #rho_array = np.append(rho_array, str(rho))
    #box_array = np.append(box_array, str(box))
    #dt_sn_array = np.append(dt_sn_array, str(dt_sn))
    #radius_sn_array = np.append(radius_sn_array, str(R_inj))

    i = i+1

legend = np.array(['Gas surface density [M_sol pc^-2]', 'SN rate[Myr^-1 kpc^-2]', 'number density [cm^-3]', 'ambient density [g cm^-3]', 'SN rate[10^-12year^-1pc^-3]', 'Time between SNe [secs]', 'Cooling radius [pc]', 'SN driving scale (iso) [pc]', 'SN Driving Scale (strat) [pc]', 'velocity dispersion [km s^-1]', 'ratio of L_drive/R_cool', 'Box size [pc]', 'Injection radius [pc]', 'dt_sn_reset', 'Scale height [pc]', 'H [pc]', 'eq 11 n_0', 'sigma_g/2H'])
output = np.row_stack((sigma_g, sigma_sn, new_n0_array, rho_array, ndot_array, dt_sn_array, R_st_array, L_array, velocity_array, LtoR_ratio, box_array, radius_sn_array))
output1 = np.row_stack((sigma_g, sigma_sn, new_n0_array))
output2 = np.row_stack((rho_array))
output3 = np.row_stack((ndot_array, dt_sn_array, R_st_array, L_array, velocity_array, LtoR_ratio, box_array, radius_sn_array))
#legend = np.vstack(legend)
#output = np.column_stack((legend, output))
i=0
j=0
#print(legend[i] +':   '+ str(sigma_g))
#print(legend[i+1] +':   '+ str(sigma_sn))
#print(legend[i+2] +':  '+ str(new_n0_array))
#print(legend[i+3] +':   '+ str(rho_array))
#print(legend[i+4] +':   '+ str(ndot_array))
#print(legend[i+5] +':   '+ str(dt_sn_array))
#print(legend[i+6] +':   '+ str(R_st_array))
#print(legend[i+7] +':   '+ str(L_array))
#print(legend[i+8] +':   '+ str(velocity_array))
#print(legend[i+9] +':   '+ str(LtoR_ratio))
#print(legend[i+10] +':   '+ str(box_array))
#print(legend[i+11] +':   '+ str(radius_sn_array))
print('-----')
print('{}:   {}').format(legend[i], sigma_g)
print('{}:   {}').format(legend[i+1], sigma_sn)
print legend[i+2], ['{:0.2f}'.format(x) for x in new_n0_array]
print legend[i+3], ['{:0.2e}'.format(x) for x in rho_array]
print legend[i+4], ['{:0.2f}'.format(x) for x in ndot_array]
print legend[i+5], ['{:0.2e}'.format(x) for x in dt_sn_array]
print legend[i+6], ['{:0.3f}'.format(x) for x in R_st_array]
print legend[i+7], ['{:0.3f}'.format(x) for x in L_array]
print legend[i+8], ['{:0.3f}'.format(x) for x in L_strat_array]
print legend[i+9], ['{:0.3f}'.format(x) for x in velocity_array]
print legend[i+10], ['{:0.3f}'.format(x) for x in LtoR_ratio]
print legend[i+11], ['{:0.3f}'.format(x) for x in box_array]
print legend[i+12], ['{:0.3f}'.format(x) for x in radius_sn_array]
print legend[i+13], ['{:0.2e}'.format(x) for x in dt_sn_reset_array]
print legend[i+14], ['{:0.2f}'.format(x) for x in scale_height_array]
print legend[i+15], ['{:0.2f}'.format(x) for x in H_array]
print legend[i+16], ['{:0.2f}'.format(x) for x in n0_11_array]
print legend[i+17], ['{:0.2f}'.format(x) for x in new_rho_array]

#print('{}:   {:.4e}').format(legend[i+2], new_n0_array)


#print(legend[i] + ' ' + sigma_g[i])

#f_handle=open("output.txt",'a')
#np.savetxt(f_handle, output1, delimiter='  ', fmt='%.3f')
#np.savetxt(f_handle, np.hstack(output2), delimiter='  ', fmt='%.3e')
#np.savetxt(f_handle, output3, delimiter='  ', fmt='%.3f')
#f_handle.close()
