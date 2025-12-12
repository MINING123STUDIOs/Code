print("Setting up simulation. . . (advanced simulation model 2025)")

# TODO: -/-
# This is a universal ODE simulation framework.
#importing standart libs.
import matplotlib.pyplot as plt
import numpy as np
import scipy as sci
import scipy.optimize as op
import random as rng
import string
import sys
import hashlib
from cryptography.fernet import Fernet
from functions import *

lb()

set_const()

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
SuppHash = "b1fe0a2f2ac08455db72eba8a45d7999"
#---
m1 = 1
m2 = 1

l1 = 1
l2 = 1

g = 9.81

#current example: double pendulum

#calculating secondary params
burnsteps =  -int(- Burntime / DeltaTime) #number of required timesteps to burn
steps = -int(-( Timeframe )/DeltaTime) #number of required timesteps
progress_bar_update_time = max( 1, int( ( steps + burnsteps ) / 100 ) )
ODEsolver = ODEsolver.replace("eE", "explicit_Euler").replace("iE", "implicit_Euler").replace("GLRK", "Gauss_Legendre_Runge_Kutta_").replace("RK", "Runge_Kutta_").replace("_", " ")
FILECONTENT = "\n\n#sim.py: \n\n" + readfile("sim.py").replace(SuppHash, "") + "\n\n#functions.py: \n\n" + readfile("functions.py") + "\n\n#config.ini: \n\n" + readfile("config.ini")
FC = FILECONTENT
config = readfile("config.ini")
FileHash = MD5(FILECONTENT)
#---



#working variables
TIME = np.linspace(Burntime, Burntime + Timeframe, steps)
Rec = np.zeros( steps )
#---
Time = 0
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

scope = globals()

UI(steps, burnsteps, ODEsolver, EQsolver, Save_Data, Save_Format, Save_Filename, DeltaTime, Burntime, Timeframe, FileHash, SuppHash, Enable_console, Confirm_num_len, scope)

for x1 in range( steps + burnsteps ):
    Time += DeltaTime #keeping time
    #--- stuff VVV
    State = f(Time, dState, State)
    dState = step(df, Time, dState, State, ODEsolver, DeltaTime, EQsolver)
    if x1 >= burnsteps: #recording data
        Rec[x1 - burnsteps ] = dState[0]
        TIME[x1 - burnsteps ] = dState[1]
    if x1 % progress_bar_update_time == 0:
        progress(100 * x1 / ( steps + burnsteps))
        pass

mini_UI()

# post processing

x1 = + l1 * np.sin(Rec)
y1 = - l1 * np.cos(Rec)

x2 = x1 + l2 * np.sin(TIME)
y2 = y1 - l2 * np.cos(TIME)

Rec = y2
TIME = x2

#plotting data

t = TIME
s = np.nan_to_num(Rec)

for ___ in range(3): #basically redundant
    if len(s) > len(t):
        t = np.append(t, t[-1] + DeltaTime)
    
    if len(t) > len(s):
        s = np.append(s, s[-1])

if len(s) != len(t):
     lb() 
     print("\nData Error!")
     input("press Enter to exit.")
     exit()

D = np.array([t,s])

save_file(Save_Data, Save_Format, Save_Filename, D)

plot(Plot, t, s)

input("The end of the programm was reached. Press enter to exit.")

"""NOTES:
        Note 1:  modifies input array instead of making a new one to improve performance. Due to this being at the end and working one the local copy of dState it does NOT mutate the simulation.
        ⚠ *SnsymSP
"""
