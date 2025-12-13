print("Setting up simulation. . . (advanced simulation model 2025)")
# This is a universal ODE simulation framework.
#importing standart libs.
import matplotlib.pyplot as plt
import numpy as np
import string
import time
from functions import *

scope = globals()
dlb = "\n\n"

#sim params
Timeframe = 5 #s
Burntime = 0 #s
DeltaTime = 2e-3 #s
Num_Damp = 1
ODEsolver = "GLRK4" # implicit_Euler / iE / explicit_Euler / eE / Runge_Kutta_4 / RK4 / Gauss_Legendre_Runge_Kutta_4 / GL4
EQsolver = "fS" #custom_Newton / cN / fSolve / fS
Save_Data = False
Save_Format = ".csv" # .csv, .txt, .npz
Save_Filename = "Recording"
Enable_console = True
Confirm_num_len = 8
Plot = "Graph"
SuppHash = "def65f0076310ce873c9a248ab99474b"
set_const(scope, SuppHash)
#---
m1 = 1 #current example: double pendulum
m2 = 1
l1 = 1
l2 = 1
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
    #np.nan_to_num(x)
    return x #READ NOTE 1!!!

def f(t, x, s):
    
    return s # (U_A, Rec)

def Rec_fun(State, dState):
    Rec  = dState[0]
    TIME = dState[1]
    return Rec, TIME

UI(ODEsolver, EQsolver, Save_Data, Save_Format, Save_Filename, DeltaTime, Burntime, Timeframe, FileHash, SuppHash, Enable_console, Confirm_num_len, scope)

TIME, Rec = run_sim(DeltaTime, State, dState, Timeframe, Burntime, f, df, ODEsolver, EQsolver, Rec_fun, True)

# post processing
Rec, TIME = - l1 * np.cos(Rec) - l2 * np.cos(TIME), + l1 * np.sin(Rec) + l2 * np.sin(TIME)

#plotting data

t = TIME
s = np.nan_to_num(Rec)
D = np.array([t,s])

save_file(Save_Data, Save_Format, Save_Filename, D)

plot(Plot, t, s)

input("The end of the programm was reached. Press enter to exit.")

"""
NOTES:
Note 1:  modifies input array instead of making a new one to improve performance. Due to this being at the end and working one the local copy of dState it does NOT mutate the simulation.
"""
