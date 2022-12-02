# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 14:17:44 2022

@author: jiri.nabelek
"""
import numpy as np
import math
import matplotlib.pyplot as plt

def plot_hall(s_pedestrian, q_pedestrian, hall_length, left_wall, right_wall, one_sided_diff, one_sided_door_diff):
    figure, axes = plt.subplots()
    plt.plot([left_wall, left_wall], [0,hall_length], color="black")
    plt.plot([right_wall, right_wall], [0,hall_length], color="black")
    plt.plot([left_wall, left_wall + one_sided_door_diff], [hall_length/2, hall_length/2], color="black")
    plt.plot([right_wall - one_sided_door_diff, right_wall], [hall_length/2, hall_length/2], color="black")
    plt.xlim([0, hall_length])
    plt.ylim([0, hall_length])
    plt.axis('off')
    for i in range(0,np.shape(s_pedestrian)[0]):    
        axes.add_artist(plt.Circle((left_wall+s_pedestrian[i,0], s_pedestrian[i,1]), 1, fill = False ))
        axes.add_artist(plt.Circle((left_wall+s_pedestrian[i,0] + math.cos(q_pedestrian[i,0]), s_pedestrian[i,1] + math.sin(q_pedestrian[i,0])), 0.5, fill = True, color="red" ))
    plt.savefig('./plots/hall_plot.pdf') 
    plt.show()


