print("Setting up simulation. . .") 
# This is a universal ODE simulation framework.
#importing libs.
import numpy as np
from sim_API import *

scope = globals()

#sim params
Timeframe       = 5 
Burntime        = 0 
DeltaTime       = 2e-3 
ODEsolver       = "GLRK4" 
EQsolver        = "fS" 
Save_Data       = False
Save_Format     = ".csv" 
Save_Filename   = "Recording"
Enable_console  = True
Confirm_num_len = 8
Plot            = "Graph"
SuppHash        = "fceb32a7a49ab54130b70bffbf89880c99dcaa31d8ee1334e090fbfa3d0ee383"

config, FC, FileHash, StartTime = set_const(SuppHash)
#---
m1, m2, l1, l2 = 1, 1, 1, 1 #current example: double pendulum
g = 9.81

#working variables
dState = np.array( [0.2, 0, 0, 0], dtype = np.float64 ) # (Theta1, Theta2, w1, w2) state depending on ODEs 
State = np.array( [0], dtype = np.float64 ) # (Rec) state not depending on ODEs 

def df(t, x, s):
    Theta1, Theta2, w1, w2 = x[0], x[1], x[2], x[3]
    
    delta = Theta1 - Theta2
    
    dTheta1 = w1
    dTheta2 = w2
    
    Div = ( 2 * m1 + m2 - m2 * np.cos( 2 * delta ) )
    
    dw1 = ( - g * ( 2 * m2 + m1 ) * np.sin( Theta1 ) - m2 * g * np.sin( Theta1 - 2 * Theta2 ) - 2 * np.sin( delta ) * m2 * ( w2 ** 2 * l2 + w1 ** 2 * l1 * np.cos( delta ) ) ) / ( l1 * Div )
    dw2 = ( 2 * np.sin( delta ) * ( w1 ** 2 * l1 * ( m1 + m2 ) + g * ( m1 + m2 ) * np.cos( Theta1 ) + w2 ** 2 * l2 * m2 * np.cos( delta ) ) ) / ( l2 * Div )
    
    y = np.zeros(4)
    y[0], y[1], y[2], y[3] = dTheta1, dTheta2, dw1, dw2
    return y 

UI(DeltaTime, Burntime, Timeframe, Enable_console, Confirm_num_len, scope)

Rec, Error = run_sim(DeltaTime, State, dState, Timeframe, Burntime, no_f, df, ODEsolver, EQsolver, False)

# post processing
Theta1, Theta2 = Rec[1,:], Rec[2,:]

y, x = - l1 * np.cos(Theta1) - l2 * np.cos(Theta2), + l1 * np.sin(Theta1) + l2 * np.sin(Theta2)

#plotting & saving data
save_file(Save_Data, Save_Format, Save_Filename, np.array([x,y]))

plot(Plot, x, y)
