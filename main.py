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

rep = repulsive_force_pedestrian(r_pedestrian[0], r_pedestrian[1], np.array([-1,3]), np.array([3,1]), np.array([2,1]), np.array([4,1]))

# F_0 = np.zeros(number_of_pedestrians,2)
