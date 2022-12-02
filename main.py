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
from scipy.integrate import solve_ivp, odeint

state_data = {"s": np.zeros([number_of_pedestrians, 2, len(time_discrete)]),
              "v": np.zeros([number_of_pedestrians, 2, len(time_discrete)]),
              "q": np.zeros([number_of_pedestrians, 2, len(time_discrete)]),
              "u": np.zeros([number_of_pedestrians, 2, len(time_discrete)]),
              "F": np.zeros([number_of_pedestrians, 2, len(time_discrete)])}

state_data["s"][0,:,0] = [10,5]
state_data["s"][1,:,0] = [15,8]
state_data["s"][2,:,0] = [10,25]
state_data["s"][3,:,0] = [20,35]

state_data["q"][0,:,0] = [math.pi/2,0]
state_data["q"][1,:,0] = [math.pi/3,0]
state_data["q"][2,:,0] = [math.pi/4,0]
state_data["q"][3,:,0] = [math.pi/6,0]

state_data["v"][0,:,0] = [1,1]
state_data["v"][1,:,0] = [-1,-1]
state_data["v"][2,:,0] = [0,1]
state_data["v"][3,:,0] = [1,0]

t = 0
disp.plot_hall(state_data["s"][:,:,t], state_data["q"][:,:,t], hall_length, left_wall, right_wall, one_sided_diff, one_sided_door_diff)

# Force calculation
for ped_id_i in range(0, number_of_pedestrians):
    forces_external_i = np.zeros([number_of_pedestrians,2])
    for ped_id_j in range(0, number_of_pedestrians):
        if(ped_id_i != ped_id_j):
            forces_external_i[ped_id_j, :] = repulsive_force_pedestrian(r_pedestrian[ped_id_i], r_pedestrian[ped_id_j], state_data["v"][ped_id_i,:,t], state_data["v"][ped_id_j,:,t], state_data["s"][ped_id_i,:,t], state_data["s"][ped_id_j,:,t])
    forces_external_i = np.matrix(forces_external_i)
    forces["rep_pedestrian"][ped_id_i,:] = forces_external_i.sum(axis=0)

for ped_id_i in range(0, number_of_pedestrians):
    wall_dist = [abs(state_data["s"][ped_id_i,0,t]-left_wall), abs(state_data["s"][ped_id_i,0,t]-right_wall)]
    if((left_wall + one_sided_door_diff) > state_data["s"][ped_id_i,0,t] and (right_wall - one_sided_door_diff) < state_data["s"][ped_id_i,0,t]):
        wall_dist.append(abs(state_data["s"][ped_id_i,1,t]-hall_length/2))
    else:
        wall_dist.append(np.linalg.norm(np.array(state_data["s"][ped_id_i,:,t]) - np.array([left_wall + one_sided_door_diff, hall_length/2])))
        wall_dist.append(np.linalg.norm(np.array(state_data["s"][ped_id_i,:,t]) - np.array([right_wall - one_sided_door_diff, hall_length/2])))
    forces["rep_wall"][ped_id_i,:] = repulsive_force_wall(r_pedestrian[ped_id_i], state_data["v"][ped_id_i,:,t], state_data["s"][ped_id_i,:,t], min(wall_dist))

forces["external"] = forces["rep_pedestrian"] + forces["rep_wall"]

for ped_id_i in range(0, number_of_pedestrians):
    end_point_direction = end_point - state_data["s"][ped_id_i,:,t]
    vd_i = vd*end_point_direction/np.linalg.norm(end_point_direction)
    forces["target"][ped_id_i,:] = m_pedestrian[ped_id_i]*(vd_i-state_data["v"][ped_id_i,:,t])/tau_i

forces["total"] = forces["target"] + forces["external"]

# Force projection
# for ped_id_i in range(0, number_of_pedestrians):
#     R()


# ODE solver
i = 2
F = lambda t, s: 20.0
t_eval = [time_discrete[i]]
t_span = np.linspace(time_discrete[i-1], time_discrete[i+1], 101)
v0 = [10.0]
# sol = solve_ivp(F, t_span, v0)
sol = odeint(F, v0, t_span)


# Torgue input
# theta0_i = math.atan(forces["target"][i,1]/forces["target"][i,0])
# I_i = 1/2*m_pedestrian[i]*r_pedestrian[i]
# ktheta = I_i*klambda*f0_i
# komega = I_i*(1+alpha)*math.sqrt((klambda*f0_i)/alpha)
# u0_i = -ktheta(q_pedestrian[i,0]-theta0_i)-komega*q_pedestrian[i,1]
