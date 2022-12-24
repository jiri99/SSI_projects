# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 10:26:04 2022

@author: jiri.nabelek
"""
import display as disp
from var import *
from process import *
from model import run_model
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, odeint

door_width_var = np.linspace(30, 80, 51)
number_of_pedestrians = 9

test_output = []
hall_properties = {}
hall_properties["hall_length"] = 100
hall_properties["hall_width"] = 30
flow = np.zeros([number_of_pedestrians])

# Run model
for width in door_width_var:
    hall_properties["door_width"] = width    
    hall_properties = calculate_hall_params(hall_properties)
    pedestrian_var = create_pedestrian_var(number_of_pedestrians)
    
    test_output.append(run_model(False, hall_properties, pedestrian_var, number_of_pedestrians))

# Static calculation
for test_id in range(0, len(test_output)-1):
    ped_id_through = np.zeros([number_of_pedestrians])
    for ped_id in range(0, number_of_pedestrians):
        x_ped = test_output[test_id]["s"][ped_id,1,:]
        ped_id_through[ped_id] = np.where(x_ped > hall_properties["hall_length"]/2)[0][0]
    time_diff = time_discrete[max(ped_id_through)] - time_discrete[min(ped_id_through)]
    
    flow[test_id] = number_of_pedestrians/time_diff
    
door_width_m = hall_properties["hall_width"]*(door_width_var/100)

plt.scatter(door_width_m, flow, c ="blue",
            linewidths = 2,
            marker ="s",
            edgecolor ="purple",
            s = 50)

plt.xlabel("Width of the narrowed part of the corridor [m]")
plt.ylabel("Flow [ped/s]")
plt.show()
