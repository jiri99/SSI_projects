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

# Repulsive force to another pedestrian
def repulsive_force_pedestrian(r_i, r_j, v_i, v_j):
    r_ij = r_i + r_j
    d_ij = np.linalg.norm(r_i-r_j)
    n_ij = [(r_i-r_j)[i]/np.linalg.norm(r_i-r_j) for i in range(0,2)]
    t_ij = np.transpose([-n_ij[1], n_ij[0]])
    delta_v_ij =  np.transpose(v_j-v_i)*t_ij
    Fp_ij = lambda r_ij,d_ij,n_ij : (A_i*math.exp((r_ij-d_ij)/B_i)+k_1*max(r_ij-d_ij,0))*n_ij+k_2*max(r_ij-d_ij,0)*delta_v_ij*t_ij
    return Fp_ij(r_ij, d_ij, n_ij)
    

# Repulsive force to wall
def repulsive_force_wall(r_i, r_w, v_i):
    r_iw = r_i + r_w
    d_iw = np.linalg.norm(r_i-r_w)
    n_iw = [(r_i-r_w)[i]/np.linalg.norm(r_i-r_w) for i in range(0,2)]
    t_iw = np.transpose([-n_iw[1], n_iw[0]])
    delta_v_iw =  np.transpose(0-v_i)*t_iw
    Fw_iw = lambda r_iw,d_iw,n_iw : (A_w*math.exp((r_iw-d_iw)/B_w)+k_1*max(r_iw-d_iw,0))*n_iw+k_2*max(r_iw-d_iw,0)*delta_v_iw*t_iw
    return Fw_iw(r_iw, d_iw, n_iw)
