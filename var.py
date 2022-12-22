# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 15:00:25 2022

@author: jiri.nabelek
"""
import numpy as np
import math


# Hall properties
def calculate_hall_params(hall_properties):
    hall_properties["one_sided_diff"] = hall_properties["hall_width"]/2
    hall_properties["one_sided_door_diff"] = (2*hall_properties["one_sided_diff"]*(100-hall_properties["door_width"])/100)/2
    hall_properties["left_wall"] = hall_properties["hall_length"]/2 - hall_properties["one_sided_diff"]
    hall_properties["right_wall"] = hall_properties["hall_length"]/2 + hall_properties["one_sided_diff"]
    hall_properties["end_point"] = np.array([hall_properties["hall_length"]/2 - hall_properties["left_wall"], hall_properties["hall_length"]])
    hall_properties["mid_point"] = np.array([hall_properties["hall_length"]/2 - hall_properties["left_wall"], hall_properties["hall_length"]/2])
    return hall_properties

# Peasant properties
def create_pedestrian_var(number_of_pedestrians):
    pedestrian_var = {}
    pedestrian_var["r_pedestrian"] = np.random.rand(number_of_pedestrians)*0.35 + 0.25
    pedestrian_var["m_pedestrian"] = np.random.randint(60, 90, number_of_pedestrians)
    return pedestrian_var


# Rotation matrix
R = lambda theta: np.array([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]])
A = np.array([[0,1],[0,0]])
b = lambda I: np.array([[0],[1/I]])

# Parameters of model
k0 = 1
kd = 500
alpha = 3
klambda = 0.3
df = 2
d0 = 1
kg_1 = 200
kg_2 = 200
vd = 5


# Parameters used in SFM
test_start = 0
test_end = 60
tau = 0.5
time_discrete = np.linspace(test_start, test_end, int((test_end-test_start)/tau+1))
A_i = 2000
A_w = 2000
B_i = 0.08
B_w = 0.08
k_1 = 120000
k_2 = 240000



