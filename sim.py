print("Setting up simulation. . . (advanced simulation model 2025)")
def lb():
    print("\n")

lb()
# TODO: -/-
#importing standart libs.
import matplotlib.pyplot as plt
import numpy as np
import scipy as sci
import random as rng
from scipy.special import lambertw
import sys

#setting constants 
sigma = 5.670374419e-8 #W/m^2*K^4
pi = np.pi

#sim params
Timeframe = 0.2 #s
Burntime = 0 #s
DeltaTime = 1e-4 #s
Num_Damp = 1
ODEsolver = "eE" # iE / implicit_Euler / eE / explicit_Euler
EQsolver = "custom_Newton" #custom_Newton / cN / fSolve / fS
#---
R = 1e3
C = 1e-7
#current example: RC filter with white noise.

#calculating secondary params
burnsteps =  -int(- Burntime / DeltaTime) #number of required timesteps to burn
steps = -int(-( Timeframe )/DeltaTime) #number of required timesteps
#---


#working variables
TIME = np.arange(Burntime, Burntime + Timeframe, DeltaTime)
Rec = np.zeros( steps )
#---
Time = 0
dState = np.array( [ 0 ], dtype = np.float64 ) # (U_c) state depending on ODEs 
State = np.array( [ 0, 0 ], dtype = np.float64 ) # (U_A, Rec) state not depending on ODEs 

def UI():
    INP = "Y"
    y = True
    while y == True:
        INP = input(f"Confirm simulating {steps + burnsteps} samples? [Y]/[N]/[i]").casefold().strip()
        if INP in [ "y", "n" ]:
            y = False
            break
        if INP == "i":
            intvar0001 = ODEsolver.replace("eE", "explicit_Euler").replace("iE", "implicit_Euler").replace("_", " ")
            intvar0002 = EQsolver.replace("cN", "custom_Newton").replace("fS", "fSolve").replace("_", " ")
            print(f"Simulating {steps + burnsteps} samples. {burnsteps} samples will be discarded ({Burntime} seconds), {steps} samples will be recorded ({Timeframe} seconds).")
            print(f"Using {intvar0001} with {intvar0002}.")
        else:
            print("Wrong input, please try again.")

    if INP == "y":
        pass
    else: # INP == "n"
        exit()

def newton_solve(F, x0, tol=1e-9, max_iter=20):
    x = x0.astype(float).copy()
    n = len(x)

    for _ in range(max_iter):
        Fx = F(x)
        if np.linalg.norm(Fx, ord=2) < tol:
            return x

        J = np.zeros((n, n), dtype=float)
        eps = 1e-8

        for i in range(n):
            x_pert = x.copy()
            x_pert[i] += eps
            J[:, i] = (F(x_pert) - Fx) / eps

        delta = np.linalg.solve(J, -Fx)
        x += delta

        if np.linalg.norm(delta, ord=2) < tol:
            return x

    return x

def linint(hi, lo, s):
    if s < 0.0:
        s = 0.0
    elif s > 1.0:
        s = 1.0
    return hi * s + lo * ( 1 - s )

def thermrad(A, emiss, T, T_amb):
    return A * emiss * sigma * ( ( T + 273.15 ) ** 4 - ( T_amb + 273.15 ) ** 4 )

def clock(t, type, hi=1, lo=-1, f=1, duty=0.5, phase=0):
    if type == "square" or type == "sq":
        return hi if ( ( f * t + phase ) % 1 ) > duty else lo
    if type == "sine" or type == "sn":
        return linint( hi, lo, 0.5 + 0.5 * np.sin( 2 * pi * ( f * t + phase ) ) )
    if type == "triangle" or type == "tg":
        return linint( hi, lo, 2 * np.abs( ( f * t + phase ) % 1 - 0.5 ) )
    if type == "saw" or type == "sw":
        return linint( hi, lo, ( f * t + phase ) % 1 )
    if type == "white_noise" or type == "wn":
        return linint( hi, lo, rng.random() )
    if type == "0":
        return linint( hi, lo, 0 )

def progress(percent):
    percent = np.clip(percent, a_min=0, a_max=100)
    bar_len = 50
    filled = int(bar_len * percent / 100)
    bar = "█" * filled + " " * (bar_len - filled)
    sys.stdout.write(f"\r|{bar}| {percent:.0f}%")
    sys.stdout.flush()

def safeexp(x):
     return np.exp(np.clip(x, a_max=700, a_min=-5e18))

def df(t, x, s):
    U = x[0]
    U_A = s[1]
    dU = ( U_A - U ) / ( R * C )
    return np.array( [ dU ], dtype = np.float64 )

def f(t, d, s):
    U_A = clock(t, "wn")
    Rec = U_A
    return np.array( [ Rec, U_A ], dtype = np.float64 ) # (U_A, Rec)

def step(t, dState, State):
    if ODEsolver == "implicit_Euler" or ODEsolver == "iE":
        def F(xn):
            return xn - dState - DeltaTime * df(t, xn, State)
        if EQsolver == "custom_Newton" or EQsolver == "cN":
            return newton_solve(F, dState)
        if EQsolver == "fSolve" or EQsolver == "fS":
            return sci.optimize.fsolve(F, State)
    if ODEsolver == "explicit_Euler" or ODEsolver == "eE":
        return State + DeltaTime * df(t, dState, State)

dState = step(Time, dState, State)
State = f(Time, dState, State)

print("Done.")
UI()
print("0/1 Complete.")

for x1 in range( steps + burnsteps ):
    Time += DeltaTime #keeping time
    #--- stuff VVV
    State = f(Time, dState, State)
    dState = step(Time, dState, State)
    if x1 >= burnsteps: #recording data
        Rec[x1 - burnsteps ] = dState[0]
    progress(100 * x1 / ( steps + burnsteps))
lb()
print("1/1 Complete.")
print("Please wait. . .")
# post processing



#plotting data

t = TIME
s = np.nan_to_num(Rec)
for ___ in range(3):
    if len(s) > len(t):
        t = np.append(t, t[-1] + DeltaTime)
    
    if len(t) > len(s):
        s = np.append(s, s[-1])

if len(s) != len(t):
     lb()
     print("\nData Error!")
     input("press Enter to exit.")
     exit()

fig, ax = plt.subplots()
ax.plot(t, s)

ax.set(xlabel='Time in s', ylabel='Y-axis',
       title='Diagram')
#ax.set_yscale("log")
plt.tick_params(axis="both", which="both")
ax.grid()

fig.savefig("test.png")
print("Done.")
plt.show()
