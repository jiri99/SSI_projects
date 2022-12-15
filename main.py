# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 21:25:15 2022

@author: jiri.nabelek
"""
from var import *
from model import run_model
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, odeint

hall_properties = {}
hall_properties["hall_length"] = 100
hall_properties["hall_width"] = 30
hall_properties["door_width"] = 30
hall_properties = calculate_hall_params(hall_properties)

state_data = run_model(True, hall_properties)