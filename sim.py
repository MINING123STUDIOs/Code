print("Setting up simulation. . .") 
# This is a demonstration of the API.
#importing libs.
import numpy as np
from     sim_API import *
from special_API import *

scope = globals()

#sim params
timeframe        = 5 
burntime         = 0 
delta_time       = 2e-3 
ode_solver       = "iE" 
eq_solver        = "fS" 
save_data        = False
save_format      = ".csv" 
save_filename    = "Recording"
enable_console   = True
plot_type        = "Graph"
supp_hash        = "fceb32a7a49ab54130b70bffbf89880c99dcaa31d8ee1334e090fbfa3d0ee383"

set_const(scope)
#---
m1, m2, l1, l2 = 1, 1, 1, 1 #current example: double pendulum
g = 9.81

#working variables
d_state = np.array( [0.2, 0, 0, 0] ) # (Theta1, Theta2, w1, w2) state depending on ODEs 
state = np.array( [0] ) #state not depending on ODEs 

def df(t, x, s): return double_pendulum_ode(t, x, s, m1, m2, l1, l2, g)    

UI(scope, enable_console)

#rec = run_sim(delta_time, state, d_state, timeframe, burntime, no_f, df, ode_solver, eq_solver, False)
rec = adt_run_sim(delta_time, state, d_state, timeframe, burntime, no_f, df, ode_solver, 1, 0, 500)

# post processing
theta1, theta2 = rec[1,:], rec[2,:]

x, y = + l1 * np.sin(theta1) + l2 * np.sin(theta2), - l1 * np.cos(theta1) - l2 * np.cos(theta2)

#plotting & saving data
save_file(save_data, save_format, save_filename, np.array([x,y]))

plot(plot_type, np.array([x,y]), ["a"])