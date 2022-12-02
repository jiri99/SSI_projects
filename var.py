# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 15:00:25 2022

@author: jiri.nabelek
"""
import numpy as np
import math


# Hall properties 
hall_length = 100
hall_width = 30
door_width = 30

one_sided_diff = hall_width/2
one_sided_door_diff = (2*one_sided_diff*(100-door_width)/100)/2
left_wall = hall_length/2 - one_sided_diff
right_wall = hall_length/2 + one_sided_diff
end_point = np.array([hall_length/2, hall_length])


# Peasant properties
number_of_pedestrians = 4
r_pedestrian = np.random.rand(number_of_pedestrians)*0.35 + 0.25
m_pedestrian = np.random.randint(60, 90, number_of_pedestrians)
F_0 = np.zeros([number_of_pedestrians, 2])
forces = {"rep_pedestrian": np.zeros([number_of_pedestrians, 2]), 
          "rep_wall": np.zeros([number_of_pedestrians, 2]),
          "external": np.zeros([number_of_pedestrians, 2]),
          "target": np.zeros([number_of_pedestrians, 2])}
orientation_pedestrians = np.zeros([number_of_pedestrians, 2])


# Rotation matrix
R = lambda theta: np.array([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]])


# Parameters of model
k0 = 1
kd = 500
alpha = 3
klambda = 0.3
df = 2
d0 = 1
kg_1 = 200
kg_2 = 200
vd = 1.5


# Parameters used in SFM
test_start = 0
test_end = 20
tau_i = 0.5
time_discrete = np.linspace(test_start, test_end, int((test_end-test_start)/tau_i+1))
A_i = 2000
A_w = 2000
B_i = 0.08
B_w = 0.08
k_1 = 120000
k_2 = 240000



