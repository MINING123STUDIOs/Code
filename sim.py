print("Setting up simulation. . .") 
# This is a universal ODE simulation framework.
#importing libs.
import numpy as np
from sim_API import *

scope = globals()

#sim params
timeframe        = 5 
burntime         = 0 
delta_time       = 2e-3 
ode_solver       = "GLRK2" 
eq_solver        = "fS" 
save_data        = False
save_format      = ".csv" 
save_filename    = "Recording"
enable_console   = True
confirm_num_len  = 8
plot_type        = "Graph"
supp_hash        = "fceb32a7a49ab54130b70bffbf89880c99dcaa31d8ee1334e090fbfa3d0ee383"

config, fc, file_hash, start_time = set_const(supp_hash)
#---
m1, m2, l1, l2 = 1, 1, 1, 1 #current example: double pendulum
g = 9.81

#working variables
d_state = np.array( [0.2, 0, 0, 0], dtype = np.float64 ) # (Theta1, Theta2, w1, w2) state depending on ODEs 
state = np.array( [0], dtype = np.float64 ) # (Rec) state not depending on ODEs 

def df(t, x, s):
    theta1, theta2, w1, w2 = x[0], x[1], x[2], x[3]
    
    delta = theta1 - theta2
    
    dtheta1 = w1
    dtheta2 = w2
    
    Div = ( 2 * m1 + m2 - m2 * np.cos( 2 * delta ) )
    
    dw1 = ( - g * ( 2 * m2 + m1 ) * np.sin( theta1 ) - m2 * g * np.sin( theta1 - 2 * theta2 ) - 2 * np.sin( delta ) * m2 * ( w2 ** 2 * l2 + w1 ** 2 * l1 * np.cos( delta ) ) ) / ( l1 * Div )
    dw2 = ( 2 * np.sin( delta ) * ( w1 ** 2 * l1 * ( m1 + m2 ) + g * ( m1 + m2 ) * np.cos( theta1 ) + w2 ** 2 * l2 * m2 * np.cos( delta ) ) ) / ( l2 * Div )
    
    dx = np.zeros_like(x)
    dx[0], dx[1], dx[2], dx[3] = dtheta1, dtheta2, dw1, dw2
    return dx 

UI(scope, enable_console, confirm_num_len)

rec = run_sim(delta_time, state, d_state, timeframe, burntime, no_f, df, ode_solver, eq_solver, False)

# post processing
theta1, theta2 = rec[1,:], rec[2,:]

x, y = + l1 * np.sin(theta1) + l2 * np.sin(theta2), - l1 * np.cos(theta1) - l2 * np.cos(theta2)

#plotting & saving data
save_file(save_data, save_format, save_filename, np.array([x,y]))

plot(plot_type, np.array([x,y]), ["a","b"])
