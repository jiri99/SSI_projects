# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 21:25:15 2022

@author: jiri.nabelek
"""
import display as disp
from var import *
import numpy as np
import math
import matplotlib.pyplot as plt

def plot_hall(s_pedestrians, orientation_pedestrians, left_wall, right_wall, one_sided_diff, one_sided_door_diff):
    figure, axes = plt.subplots()
    plt.plot([left_wall, left_wall], [0,100], color="black")
    plt.plot([right_wall, right_wall], [0,100], color="black")
    plt.plot([left_wall, left_wall + one_sided_door_diff], [50,50], color="black")
    plt.plot([right_wall - one_sided_door_diff, right_wall], [50,50], color="black")
    plt.xlim([0, hall_length])
    plt.ylim([0, hall_length])
    plt.axis('off')
    for i in range(0,np.shape(s_pedestrians)[0]):    
        axes.add_artist(plt.Circle((left_wall+s_pedestrians[i,0], s_pedestrians[i,1]), 1, fill = False ))
        axes.add_artist(plt.Circle((left_wall+s_pedestrians[i,0], s_pedestrians[i,1]+1), 0.5, fill = True, color="red" ))
    plt.savefig('./plots/hall_plot.pdf') 
    plt.show()

s_pedestrians[0,:] = [10,5]
s_pedestrians[1,:] = [20,15]
s_pedestrians[2,:] = [10,25]
s_pedestrians[3,:] = [20,35]

disp.plot_hall(s_pedestrians, orientation_pedestrians, left_wall, right_wall, hall_length, one_sided_diff, one_sided_door_diff)

# Repulsive force to another pedestrian
def repulsive_force_pedestrian(r_i, r_j):
    r_ij = r_i + r_j
    d_ij = np.norm(r_i-r_j)
    n_ij = [(r_i-r_j)[i]/np.norm(r_i-r_j) for i in range(0,2)]
    t_ij = np.transpose([-n_ij[1], n_ij[0]])
    delta_v_ij = (v_j-v_i)*t_ij
    Fp_ij = lambda r_ij,d_ij,n_ij : (A_i*math.exp((r_ij-d_ij)/B_i)+k_1*max(r_ij-d_ij,0))*n_ij+k_2*max(r_ij-d_ij,0)*delta_v_ij*t_ij
    return Fp_ij(r_ij, d_ij, n_ij)
    

# Repulsive force to wall
def repulsive_force_wall(r_i, r_w):
    r_iw = r_i + r_w
    d_iw = np.norm(r_i-r_w)
    n_iw = [(r_i-r_w)[i]/np.norm(r_i-r_w) for i in range(0,2)]
    t_iw = np.transpose([-n_iw[1], n_iw[0]])
    delta_v_iw = (0-v_i)*t_iw
    Fw_iw = lambda r_iw,d_iw,n_iw : (A_w*math.exp((r_iw-d_iw)/B_w)+k_1*max(r_iw-d_iw,0))*n_iw+k_2*max(r_iw-d_iw,0)*delta_v_iw*t_iw
    return Fw_iw(r_iw, d_iw, n_iw)

# F_0 = np.zeros(number_of_pedestrians,2)
