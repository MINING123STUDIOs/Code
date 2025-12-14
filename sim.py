print("Setting up simulation. . . (advanced simulation model 2025)")
# This is a universal ODE simulation framework.
#importing standart libs.
import numpy as np
from sim_API import *

scope = globals()

#sim params
Timeframe       = 5 
Burntime        = 0 
DeltaTime       = 2e-3 
Num_Damp        = 1
ODEsolver       = "GLRK4" 
EQsolver        = "fS" 
Save_Data       = False
Save_Format     = ".csv" 
Save_Filename   = "Recording"
Enable_console  = True
Confirm_num_len = 8
Plot            = "Graph"
SuppHash        = "6162a3a585ad85aa7096f0b64862e0be"
config, FC, FileHash, StartTime = set_const(SuppHash)
#---
m1, m2, l1, l2 = 1, 1, 1, 1 #current example: double pendulum
g = 9.81

#working variables
dState = np.array( [0.2, 0, 0, 0], dtype = np.float64 ) # (Theta1, Theta2, w1, w2) state depending on ODEs 
State = np.array( [ 0 ], dtype = np.float64 ) # (Rec) state not depending on ODEs 

def df(t, x, s):
    Theta1, Theta2, w1, w2 = x[0], x[1], x[2], x[3]
    
    delta = Theta1 - Theta2
    
    dTheta1 = w1
    dTheta2 = w2
    
    Div = ( 2 * m1 + m2 - m2 * np.cos( 2 * delta ) )
    
    dw1 = ( - g * ( 2 * m2 + m1 ) * np.sin( Theta1 ) - m2 * g * np.sin( Theta1 - 2 * Theta2 ) - 2 * np.sin( delta ) * m2 * ( w2 ** 2 * l2 + w1 ** 2 * l1 * np.cos( delta ) ) ) / ( l1 * Div )
    dw2 = ( 2 * np.sin( delta ) * ( w1 ** 2 * l1 * ( m1 + m2 ) + g * ( m1 + m2 ) * np.cos( Theta1 ) + w2 ** 2 * l2 * m2 * np.cos( delta ) ) ) / ( l2 * Div )
    
    x[0], x[1], x[2], x[3] = dTheta1, dTheta2, dw1, dw2
    return x #read note 1.

def Rec_fun(State, dState):
    Rec  = dState[0]
    TIME = dState[1]
    return Rec, TIME

UI(ODEsolver, EQsolver, Save_Data, Save_Format, Save_Filename, DeltaTime, Burntime, Timeframe, FileHash, SuppHash, Enable_console, Confirm_num_len, scope)

TIME, Rec = run_sim(DeltaTime, State, dState, Timeframe, Burntime, no_f, df, ODEsolver, EQsolver, Rec_fun, True)

# post processing
Rec, TIME = - l1 * np.cos(Rec) - l2 * np.cos(TIME), + l1 * np.sin(Rec) + l2 * np.sin(TIME)

#plotting data
save_file(Save_Data, Save_Format, Save_Filename, np.array([TIME,Rec]))

plot(Plot, TIME, Rec)

#Note 1:  modifies input array instead of making a new one to improve performance. Due to this being at the end and working one the local copy of dState it does NOT mutate the simulation.
