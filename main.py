# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 21:25:15 2022

@author: jiri.nabelek
"""
import display as disp
import numpy as np
import math
import matplotlib.pyplot as plt

def plot_hall(s_peasants, orientation_peasants, left_wall, right_wall, one_sided_diff, one_sided_door_diff):
    figure, axes = plt.subplots()
    plt.plot([left_wall, left_wall], [0,100], color="black")
    plt.plot([right_wall, right_wall], [0,100], color="black")
    plt.plot([left_wall, left_wall + one_sided_door_diff], [50,50], color="black")
    plt.plot([right_wall - one_sided_door_diff, right_wall], [50,50], color="black")
    plt.xlim([0, hall_length])
    plt.ylim([0, hall_length])
    plt.axis('off')
    for i in range(0,np.shape(s_peasants)[0]):    
        axes.add_artist(plt.Circle((left_wall+s_peasants[i,0], s_peasants[i,1]), 1, fill = False ))
        axes.add_artist(plt.Circle((left_wall+s_peasants[i,0], s_peasants[i,1]+1), 0.5, fill = True, color="red" ))
    plt.savefig('./plots/hall_plot.pdf') 
    plt.show()


hall_length = 100

number_of_peasants = 4
hall_width = 30
door_width = 30

one_sided_diff = (hall_length*hall_width/100)/2
one_sided_door_diff = (2*one_sided_diff*(100-door_width)/100)/2
left_wall = hall_length/2 - one_sided_diff
right_wall = hall_length/2 + one_sided_diff

F_0 = np.zeros([number_of_peasants, 2])
s_peasants = np.zeros([number_of_peasants, 2])
orientation_peasants = np.zeros([number_of_peasants, 2])

s_peasants[0,:] = [10,5]
s_peasants[1,:] = [20,15]
s_peasants[2,:] = [10,25]
s_peasants[3,:] = [20,35]

disp.plot_hall(s_peasants, orientation_peasants, left_wall, right_wall, hall_length, one_sided_diff, one_sided_door_diff)


# F_0 = np.zeros(number_of_peasants,2)
