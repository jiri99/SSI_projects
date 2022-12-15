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

door_width_var = np.linspace(20, 60, 5)

test_output = []
hall_properties = {}
hall_properties["hall_length"] = 100
hall_properties["hall_width"] = 30

# Run model
for width in door_width_var:
    hall_properties["door_width"] = width    
    hall_properties = calculate_hall_params(hall_properties)
    
    test_output.append(run_model(False, door_width))

for test_id in test_output:
    test_output[test_id]
    
    # Static calculation

