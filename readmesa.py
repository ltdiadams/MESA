#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 11:30:23 2017

@author: clovekin
"""

import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy import constants as const
from astropy import units as u
import mesa_reader as ms

def read_profile(direc,fname,labelName):

    fname = str(direc) + '/' +  str(fname)
    file = open(fname)

    data = file.readlines()

    file.close()

    global_props = data[1]
    global_data = data[2]


    radial_props = data[5].strip().split()

    radial_data = np.loadtxt(fname,skiprows=6)

    print(radial_props)
    data_table = Table(radial_data,names = radial_props)#, meta = {'name','first table'})


    gradT = np.diff(data_table['temperature'])*u.K/(np.diff(data_table['radius'])*const.R_sun.cgs)

    F_rad = -(16*const.sigma_sb.cgs*(data_table['temperature'][1:]*u.K)**3)/(3*(10**data_table['logRho']*u.g/u.cm**3)[1:])*(1/(10**data_table['log_opacity'][1:])*u.g/u.cm**2)*gradT

    L_rad = F_rad*4*np.pi*((10**data_table['logR'][1:])*const.R_sun.cgs)**2

    # find a zero crossing for log_conv_vel

    zero_crossings = np.where(np.diff(np.sign(data_table['log_conv_vel'])))[0]
    mean_Lrad = np.mean(L_rad[:zero_crossings[0]]/(data_table['luminosity'][:zero_crossings[0]]*const.L_sun.cgs))
   # plt.plot(10**data_table['logR'],data_table['log_conv_vel'],label = "log v$_{conv}$")
   # plt.plot(10**data_table['logR']/np.max(10**data_table['logR']),data_table['log_D_conv'], label = "log D$_{conv}$")
#plt.plot(data_table['P_initial']/data_table['P_after'], data_table['P_div'], label = "Pinit vs Ptot")

    #plt.plot(10**data_table['logR'], data_table['P_div'], '-o', label = "Pinit vs Ptot")
    #plt.plot(10**data_table['logR'], data_table['P'], '-o', label = "With P_mag")
    #plt.plot(10**data_table['logR'], data_table['P_without'], '-o', label = "Without P_mag")
    plt.plot(10**data_table['logR'], data_table['log_D_mix'], '-o', label = labelName)
    
    
    
    plt.xlabel("Radius(Rsun)")
    plt.ylabel("log10 diffusion coefficient for mixing")
    
    plt.legend(loc = "upper_right")

    plt.title("log_D_mix VS radius")

    #plt.plot(10**data_table['logR'],mean_Lrad*np.ones(len(data_table['logR'])))
    #plt.plot(10**data_table['logR'][1:],L_rad/(data_table['luminosity'][1:]*const.L_sun.cgs))
    #plt.axis([0.6,1.8,-20,20])

    #print 10**data_table['logR'][zero_crossings[0]]

def make_HR(direc,numTracks):

    for i in range(1,numTracks+1):
        #profilefile = direc + "/i"+str(i).zfill(3)+"_profiles.index"
#        profilefile = direc + "profiles.index"
#        prof_info = np.loadtxt(profilefile,skiprows=1)
#        print "how many lines?", prof_info.size, profilefile
#        if prof_info.size < 4:
#            continue
#        modnum = prof_info[:,0]
#        profnum = prof_info[:,2]
#        firstpoint = np.int(modnum[profnum==2][0])-1
        history_file = "history.data"
        star = ms.MesaData(file_name=direc)
        effectiveT = star.data("log_Teff")
        logL = star.data("log_L")
    plt.plot(effectiveT,logL,'k-',linewidth = 2)
        #for profile in profnum:
#            if np.mod(profile,2) == 0:
#                index_num = np.int(modnum[profnum == profile][0]) - 1
#                print index_num
#                plt.plot(np.log10(effectiveT[index_num]),logL[index_num],'kd')
    plt.gca().invert_xaxis()


def evolution(param1,param2,filename):

    star = ms.MesaData(filename)

    plot1 = star.data(param1)
    if param2 == 'rhobar':
        mass = star.data('star_mass')
        mass = mass*1.989e33
        r = star.data('radius_cm')
        plot2 = 3 * mass / (4 * np.pi * r**3)
    else:
        plot2 = star.data(param2)
    if param1 == 'star_mass':
        plot1 = plot1 * 1.989e33
    if param2 == 'star_mass':
        plot2 = plot2 * 1.989e33


    plt.plot(plot1,plot2)
    plt.xlabel(param1)
    plt.ylabel(param2)
