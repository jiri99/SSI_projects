# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 14:17:44 2022

@author: jiri.nabelek
"""
import numpy as np
import math
import matplotlib.pyplot as plt
import time


def plot_hall(t, s_pedestrian, q_pedestrian, hall_properties, save_fig):
    figure, axes = plt.subplots()
    plt.plot([hall_properties["left_wall"], hall_properties["left_wall"]], [0,hall_properties["hall_length"]], color="black")
    plt.plot([hall_properties["right_wall"], hall_properties["right_wall"]], [0,hall_properties["hall_length"]], color="black")
    plt.plot([hall_properties["left_wall"], hall_properties["left_wall"] + hall_properties["one_sided_door_diff"]], [hall_properties["hall_length"]/2, hall_properties["hall_length"]/2], color="black")
    plt.plot([hall_properties["right_wall"] - hall_properties["one_sided_door_diff"], hall_properties["right_wall"]], [hall_properties["hall_length"]/2, hall_properties["hall_length"]/2], color="black")
    plt.xlim([0, hall_properties["hall_length"]])
    plt.ylim([0, hall_properties["hall_length"]])
    plt.axis('off')
    for i in range(0,np.shape(s_pedestrian)[0]):    
        axes.add_artist(plt.Circle((hall_properties["left_wall"]+s_pedestrian[i,0], s_pedestrian[i,1]), 1, fill = False ))
        axes.add_artist(plt.Circle((hall_properties["left_wall"]+s_pedestrian[i,0] + math.cos(q_pedestrian[i,0]), s_pedestrian[i,1] + abs(math.sin(q_pedestrian[i,0]))), 0.5, fill = True, color="red" ))
    if(save_fig):
        plt.savefig('./plots/hall_plot_' + str(t) + '.pdf')
    plt.show()
    time. sleep(0.3)


