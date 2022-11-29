# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 16:25:24 2022

@author: jiri.nabelek
"""

import display as disp
from var import *
import numpy as np
import math
import matplotlib.pyplot as plt


# Input force
def force_input(f0_i, fe_i, rf_i, r0_i, v0_i):
    uf_i = np.transpose(f0_i+fe_i)*rf_i
    u0_i = k0*np.transpose(fe_i)*r0_i-kd*v0_i
    return np.array([uf_i, u0_i])

# Input torgue
def torgue_input(r_i, r_j, v_i, v_j, s_i, s_j):
    
    uf_i = np.transpose(f0_i+fe_i)*rf_i
    u0_i = k0*np.transpose(fe_i)*r0_i-kd*v0_i
    return np.array([uf_i, u0_i])

# Group cohesion
def group_cohesion(f0_i, fe_i, rf_i, r0_i, v0_i):
    
    if(np.norm(np.dot(np.transpose(p_i),rf_i)) > df):
        h = 1   
    else:
        h = 0
    uf_i = np.transpose(f0_i+fe_i)*rf_i + kg_1*h
    u0_i = k0*np.transpose(fe_i)*r0_i - kd*v0_i + kg_2*h
    return np.array([uf_i, u0_i])

# Repulsive force to another pedestrian
def repulsive_force_pedestrian(r_i, r_j, v_i, v_j, s_i, s_j):
    r_ij = r_i + r_j
    d_ij = np.linalg.norm(s_i-s_j)
    n_ij = np.array([(s_i-s_j)[i]/np.linalg.norm(s_i-s_j) for i in range(0,2)])
    t_ij = np.transpose([-n_ij[1], n_ij[0]])
    delta_v_ij = np.transpose(v_j-v_i)*t_ij
    Fp_ij = lambda r_ij,d_ij,n_ij : np.array([A_i*math.exp((r_ij-d_ij)/B_i)+k_1*max(r_ij-d_ij,0)])*n_ij+k_2*max(r_ij-d_ij,0)*delta_v_ij*t_ij
    return Fp_ij(r_ij, d_ij, n_ij)

# Repulsive force to wall
def repulsive_force_wall(r_i, v_i, s_i, s_w):
    r_iw = r_i + 0
    d_iw = np.linalg.norm(s_i-s_w)
    n_iw = [(s_i-s_w)[i]/np.linalg.norm(r_i-0) for i in range(0,2)]
    t_iw = np.transpose([-n_iw[1], n_iw[0]])
    delta_v_iw = np.transpose(np.array([0,0])-v_i)*t_iw
    Fw_iw = lambda r_iw,d_iw,n_iw : np.array([A_w*math.exp((r_iw-d_iw)/B_w)+k_1*max(r_iw-d_iw,0)])*n_iw+k_2*max(r_iw-d_iw,0)*delta_v_iw*t_iw
    return Fw_iw(r_iw, d_iw, n_iw)
