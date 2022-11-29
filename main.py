# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 21:25:15 2022

@author: jiri.nabelek
"""
import display as disp
from var import *
from process import *
import numpy as np
import math
import matplotlib.pyplot as plt

def plot_hall(s_pedestrian, orientation_pedestrians, hall_length, left_wall, right_wall, one_sided_diff, one_sided_door_diff):
    figure, axes = plt.subplots()
    plt.plot([left_wall, left_wall], [0,100], color="black")
    plt.plot([right_wall, right_wall], [0,100], color="black")
    plt.plot([left_wall, left_wall + one_sided_door_diff], [hall_length/2,hall_length/2], color="black")
    plt.plot([right_wall - one_sided_door_diff, right_wall], [hall_length/2,hall_length/2], color="black")
    plt.xlim([0, hall_length])
    plt.ylim([0, hall_length])
    plt.axis('off')
    for i in range(0,np.shape(s_pedestrian)[0]):    
        axes.add_artist(plt.Circle((left_wall+s_pedestrian[i,0], s_pedestrian[i,1]), 1, fill = False ))
        axes.add_artist(plt.Circle((left_wall+s_pedestrian[i,0], s_pedestrian[i,1]+1), 0.5, fill = True, color="red" ))
    plt.savefig('./plots/hall_plot.pdf') 
    plt.show()

s_pedestrian[0,:] = [10,5]
s_pedestrian[1,:] = [12,6]
s_pedestrian[2,:] = [10,25]
s_pedestrian[3,:] = [20,35]

v_pedestrian[0,:] = [1,1]
v_pedestrian[1,:] = [-1,-1]
v_pedestrian[2,:] = [0,1]
v_pedestrian[3,:] = [1,0]

disp.plot_hall(s_pedestrian, orientation_pedestrians, hall_length, left_wall, right_wall, one_sided_diff, one_sided_door_diff)


for ped_id_i in range(0, number_of_pedestrians):
    forces_external_i = np.zeros([number_of_pedestrians,2])
    for ped_id_j in range(0, number_of_pedestrians):
        if(ped_id_i != ped_id_j):
            forces_external_i[ped_id_j, :] = repulsive_force_pedestrian(r_pedestrian[ped_id_i], r_pedestrian[ped_id_j], v_pedestrian[ped_id_i], v_pedestrian[ped_id_j], s_pedestrian[ped_id_i], s_pedestrian[ped_id_j])
    forces_external_i = np.matrix(forces_external_i)
    forces["rep_pedestrian"][ped_id_i,:] = forces_external_i.sum(axis=0)

for ped_id_i in range(0, number_of_pedestrians):
    wall_dist = [abs(s_pedestrian[ped_id_i,0]-left_wall), abs(s_pedestrian[ped_id_i,0]-right_wall)]
    if((left_wall + one_sided_door_diff) > s_pedestrian[ped_id_i,0] and (right_wall - one_sided_door_diff) < s_pedestrian[ped_id_i,0]):
        wall_dist.append(abs(s_pedestrian[ped_id_i,1]-hall_length/2))
    else:
        wall_dist.append(np.linalg.norm(np.array(s_pedestrian[ped_id_i,:]) - np.array([left_wall + one_sided_door_diff, hall_length/2])))
        wall_dist.append(np.linalg.norm(np.array(s_pedestrian[ped_id_i,:]) - np.array([right_wall - one_sided_door_diff, hall_length/2])))
    forces["rep_wall"][ped_id_i,:] = repulsive_force_wall(r_pedestrian[ped_id_i], v_pedestrian[ped_id_i], s_pedestrian[ped_id_i], min(wall_dist))

forces["external"] = forces["rep_pedestrian"] + forces["rep_wall"]

for ped_id_i in range(0, number_of_pedestrians):
    end_point_direction = end_point - s_pedestrian[ped_id_i,:]
    vd_i = vd*end_point_direction/np.linalg.norm(end_point_direction)
    forces["target"][ped_id_i,:] = m_pedestrian[ped_id_i]*(vd_i-v_pedestrian[ped_id_i])/tau_i

forces["total"] = forces["target"] + forces["external"]


# F_0 = np.zeros(number_of_pedestrians,2)
