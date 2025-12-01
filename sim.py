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

#sim params
Timeframe = 0.2 #s
Burntime = 0 #s
DeltaTime = 1e-5 #s
Num_Damp = 1
#---
R = 1e4
C = 10e-6

#calculating secondary params
burnsteps =  -int(- Burntime / DeltaTime) #number of required timesteps to burn
steps = -int(-( Timeframe )/DeltaTime) #number of required timesteps
#---


#working variables
TIME = np.arange(Burntime, Burntime + Timeframe, DeltaTime)
Rec = np.zeros( steps - 1 )
#---
Time = 0
State = np.array([ 0, 0 ]) # (U, Rec)

def linint(hi, lo, s):
    s = np.clip(s, a_max = 1, a_min = 0)
    return hi * s + lo * ( 1 - s )

def thermrad(A, emiss, T, T_amb):
    return A * emiss * sigma * ( ( T + 273.15 ) ** 4 - ( T_amb + 273.15 ) ** 4 )

def clock(type, hi, lo, t, f, duty, phase):
    if type == "square":
        return hi if ( ( f * t + phase ) % 1 ) > duty else lo
    if type == "sine":
        return linint( hi, lo, 0.5 + 0.5 * np.sin( np.pi( f * t + phase ) ) )
    if type == "triangle":
        return linint( hi, lo,  )
    if type == "saw":
        return linint( hi, lo, ( f * t + phase ) % 1 )
    if type == "white_noise":
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

def f(t, x):
    U = x[0]
    U_A = clock("white_noise", 1, 0, t, 0, 0, 0)
    U = ( U_A - U ) / ( R * C )
    rec = U
    return np.array([U, rec])

def Euler_step(t, State):
    def F(xn):
        return xn - State - DeltaTime * f(t, xn)
    return sci.optimize.fsolve(F, State)

for i in range(1): #ask
    INP = "Y"
    y = True
    while y == True:
        INP = input(f"Confirm simulating {steps + burnsteps} samples? [Y]/[N]/[i]")
        if INP == "Y" or INP == "N" or INP == "i":
            y = False
        else:
            print("Wrong input, please try again.")

    if INP == "Y":
        pass
    if INP == "N":
        exit()
    if INP == "i":
        print(f"simulating {steps + burnsteps} samples. {burnsteps} samples will be discarded ({Burntime} seconds), {steps} samples will be recorded ({Timeframe} seconds),")

State = Euler_step(Time, State)

print("Done.")
print("0/1 Complete.")

for x1 in range( steps + burnsteps ):
    N_Time = Time + DeltaTime #keeping time
    #--- stuff VVV
    State = Euler_step(Time, State)
    if x1 > burnsteps: #recording data
        Rec[x1 - burnsteps - 1] = State[ 1 ]
    progress(100 * x1 / ( steps + burnsteps))
lb()
print("1/2 Complete.")

# post processing



#plotting data

t = TIME
s = np.nan_to_num(Rec)

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
plt.show()
