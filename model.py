# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 10:34:53 2022

@author: jiri.nabelek
"""
import display as disp
from var import *
from process import *
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, odeint

def run_model(plot_figure, hall_properties, pedestrian_var, number_of_pedestrians):
    state_data = {"s": np.zeros([number_of_pedestrians, 2, len(time_discrete)]),
                  "v": np.zeros([number_of_pedestrians, 2, len(time_discrete)]),
                  "vB": np.zeros([number_of_pedestrians, 2, len(time_discrete)]),
                  "a": np.zeros([number_of_pedestrians, 2, len(time_discrete)]),
                  "q": np.zeros([number_of_pedestrians, 2, len(time_discrete)]),
                  "uB": np.zeros([number_of_pedestrians, 2, len(time_discrete)]),
                  "uf": np.zeros([number_of_pedestrians, 1, len(time_discrete)]),
                  "u0": np.zeros([number_of_pedestrians, 1, len(time_discrete)]),
                  "utheta": np.zeros([number_of_pedestrians, 1, len(time_discrete)]),
                  "F": np.zeros([number_of_pedestrians, 2, len(time_discrete)]),
                  "F_0": np.zeros([number_of_pedestrians, 2, len(time_discrete)]),
                  "F_e": np.zeros([number_of_pedestrians, 2, len(time_discrete)])}
    
    forces = {"rep_pedestrian": np.zeros([number_of_pedestrians, 2]), 
          "rep_wall": np.zeros([number_of_pedestrians, 2])}
    
    state_data["s"][0,:,0] = [20,5]
    state_data["s"][1,:,0] = [20,10]
    state_data["s"][2,:,0] = [20,15]
    state_data["s"][3,:,0] = [10,5]
    state_data["s"][4,:,0] = [10,10]
    state_data["s"][5,:,0] = [10,15]
    state_data["s"][6,:,0] = [15,5]
    state_data["s"][7,:,0] = [15,10]
    state_data["s"][8,:,0] = [15,15]
    
    state_data["q"][0,:,0] = [0,0]
    state_data["q"][1,:,0] = [0,0]
    state_data["q"][2,:,0] = [0,0]
    
    state_data["v"][0,:,0] = [0,1]
    state_data["v"][1,:,0] = [0,1]
    state_data["v"][2,:,0] = [0,1]
    
    for t in range(0,len(time_discrete)-1):
        # Force calculation
        for ped_id_i in range(0, number_of_pedestrians):
            forces_external_i = np.zeros([number_of_pedestrians,2])
            for ped_id_j in range(0, number_of_pedestrians):
                if(ped_id_i != ped_id_j):
                    forces_external_i[ped_id_j, :] = repulsive_force_pedestrian(pedestrian_var["r_pedestrian"][ped_id_i], pedestrian_var["r_pedestrian"][ped_id_j], state_data["v"][ped_id_i,:,t], state_data["v"][ped_id_j,:,t], state_data["s"][ped_id_i,:,t], state_data["s"][ped_id_j,:,t])
            forces_external_i = np.matrix(forces_external_i)
            forces["rep_pedestrian"][ped_id_i,:] = forces_external_i.sum(axis=0)
        
        for ped_id_i in range(0, number_of_pedestrians):
            wall_dist = [np.array([0, state_data["s"][ped_id_i,1,t]]), np.array([hall_properties["hall_width"], state_data["s"][ped_id_i,1,t]])]
            if(hall_properties["one_sided_door_diff"] > state_data["s"][ped_id_i,0,t] or (hall_properties["hall_width"] - hall_properties["one_sided_door_diff"]) < state_data["s"][ped_id_i,0,t]):
                wall_dist.append(np.array([state_data["s"][ped_id_i,0,t], hall_properties["hall_length"]/2]))
            else:
                wall_dist.append(np.array([hall_properties["one_sided_door_diff"], hall_properties["hall_length"]/2]))
                wall_dist.append(np.array([hall_properties["hall_width"] - hall_properties["one_sided_door_diff"], hall_properties["hall_length"]/2]))
            min_wall_dist = [np.linalg.norm(np.array(dist) - state_data["s"][ped_id_i,:,t]) for dist in wall_dist]
            min_wall_dist_id = min_wall_dist.index(min(min_wall_dist))
            forces["rep_wall"][ped_id_i,:] = repulsive_force_wall(pedestrian_var["r_pedestrian"][ped_id_i], state_data["v"][ped_id_i,:,t], state_data["s"][ped_id_i,:,t], wall_dist[min_wall_dist_id])
        
        state_data["F_e"][:,:,t] = forces["rep_pedestrian"] + forces["rep_wall"]
        
        for ped_id_i in range(0, number_of_pedestrians):
            # if(state_data["s"][ped_id_i,1,t] < hall_properties["hall_length"]/2 - hall_properties["door_width"]/2):            
            #     end_point_direction = hall_properties["mid_point"] - state_data["s"][ped_id_i,:,t]
            # else:    
            #     end_point_direction = hall_properties["end_point"] - state_data["s"][ped_id_i,:,t]
            end_point_direction = hall_properties["end_point"] - state_data["s"][ped_id_i,:,t]
            vd_i = vd*end_point_direction/np.linalg.norm(end_point_direction)
            F0 = pedestrian_var["m_pedestrian"][ped_id_i]*(vd_i-state_data["v"][ped_id_i,:,t])/tau
            # state_data["F_0"][ped_id_i,:,t] = vd*F0/np.linalg.norm(F0)
            state_data["F_0"][ped_id_i,:,t] = F0
            
        
        state_data["F"][:,:,t] = state_data["F_0"][ped_id_i,:,t] + state_data["F_e"][ped_id_i,:,t]
        
        # Force projection
        for ped_id_i in range(0, number_of_pedestrians):
            
            uf_i, u0_i = force_input(state_data["F_0"][ped_id_i,1,t], state_data["F_e"][ped_id_i,:,t], R(state_data["q"][ped_id_i,0,t])[:,0], R(state_data["q"][ped_id_i,0,t])[:,1], state_data["vB"][ped_id_i,1,t])
        
            state_data["uB"][ped_id_i,:,t] = np.transpose([uf_i, u0_i])
            state_data["u0"][ped_id_i,:,t] = u0_i
            state_data["uf"][ped_id_i,:,t] = uf_i
        
            # Torgue input
            f0_i = np.linalg.norm(state_data["F_0"][ped_id_i,:,t])
            theta0_i = math.atan(state_data["F_0"][ped_id_i,1,t]/state_data["F_0"][ped_id_i,0,t])
            I_i = 1/2*pedestrian_var["m_pedestrian"][ped_id_i]*pedestrian_var["r_pedestrian"][ped_id_i]**2
            ktheta = I_i*klambda*f0_i
            komega = I_i*(1+alpha)*math.sqrt((klambda*f0_i)/alpha)
            state_data["utheta"][ped_id_i,:,t] = -ktheta*(state_data["q"][ped_id_i,0,t]-theta0_i)-komega*state_data["q"][ped_id_i,1,t]
        
        # ODE solver
        for ped_id_i in range(0, number_of_pedestrians):
            vB_dot = state_data["uB"][ped_id_i,:,t]/pedestrian_var["m_pedestrian"][ped_id_i]
            state_data["vB"][ped_id_i,:,t+1] = state_data["vB"][ped_id_i,:,t] + tau * vB_dot
            
            state_data["v"][ped_id_i,:,t+1] = np.dot(R(state_data["q"][ped_id_i,0,t]), state_data["vB"][ped_id_i,:,t])
            
            state_data["s"][ped_id_i,:,t+1] = state_data["s"][ped_id_i,:,t] + tau * state_data["v"][ped_id_i,:,t]
            
            q_dot = np.dot(A, state_data["q"][ped_id_i,:,t]) + np.dot(b(1/2*pedestrian_var["m_pedestrian"][ped_id_i]*pedestrian_var["r_pedestrian"][ped_id_i]**2), state_data["utheta"][ped_id_i,:,t])
            
            state_data["q"][ped_id_i,:,t+1] = state_data["q"][ped_id_i,:,t] + tau * q_dot
        
        if(plot_figure):
            disp.plot_hall(t, state_data["s"][:,:,t], state_data["q"][:,:,t], hall_properties, False)
    
    if(plot_figure):
        for t in [0, 40, 60, 90]:
            disp.plot_hall(t, state_data["s"][:,:,t], state_data["q"][:,:,t], hall_properties, True)

    return state_data
