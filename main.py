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


s_pedestrian[0,:] = [10,5]
s_pedestrian[1,:] = [12,6]
s_pedestrian[2,:] = [10,25]
s_pedestrian[3,:] = [20,35]

q_pedestrian[0,:] = [math.pi/2,0]
q_pedestrian[1,:] = [math.pi/3,0]
q_pedestrian[2,:] = [math.pi/4,0]
q_pedestrian[3,:] = [math.pi/6,0]

v_pedestrian[0,:] = [1,1]
v_pedestrian[1,:] = [-1,-1]
v_pedestrian[2,:] = [0,1]
v_pedestrian[3,:] = [1,0]

disp.plot_hall(s_pedestrian, q_pedestrian, orientation_pedestrians, hall_length, left_wall, right_wall, one_sided_diff, one_sided_door_diff)

# Force calculation
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

# Force projection
i = 2
F = lambda t, s: 20.0
t_eval = [time_discrete[i]]
t_span = np.linspace(time_discrete[i-1], time_discrete[i+1], 101)
v0 = [10.0]
# sol = solve_ivp(F, t_span, v0)
sol = odeint(F, v0, t_span)

# Torgue input
theta0_i = math.atan(forces["target"][i,1]/forces["target"][i,0])
I_i = 1/2*m_pedestrian[i]*r_pedestrian[i]
ktheta = I_i*klambda*f0_i
komega = I_i*(1+alpha)*math.sqrt((klambda*f0_i)/alpha)
u0_i = -ktheta(q_pedestrian[i,0]-theta0_i)-komega*q_pedestrian[i,1]


