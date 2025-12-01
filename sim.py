print("Setting up simulation. . . (advanced simulation model 2025)")
def lb():
    print("\n")

lb()
# TODO: -/-
#importing standart libs.
import matplotlib.pyplot as plt
import numpy as np
import scipy as sci
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
L_s = 10e-9
R_Al = 1e-2
R_EWET = 1e-2
R_EDRY = 50
C = 47e-3
R_P = 1e9
U_breakdown = 70
I_breakdown = 1e-5
U_beta = 3
V_s = 2

T_amb = 30
T_r = 20
R_ic = 0.5
C_i = 2
R_ca = 2
C_c = 0.7
alpha = 0.006

#calculating secondary params
burnsteps =  -int(- Burntime / DeltaTime) #number of required timesteps to burn
steps = -int(-( Timeframe )/DeltaTime) #number of required timesteps
#---


#working variables
TIME = np.arange(Burntime, Burntime + Timeframe, DeltaTime)
Rec = np.zeros( steps - 1 )
#---
Time = 0
State = np.array([ 0, 0, T_amb, T_amb, T_amb, V_s, 0 ]) # (I_L, U_int, T_i, T_c, T_A, V_c)

def thermrad(A, emiss, T, T_amb):
    return A * emiss * sigma * ( ( T + 273.15 ) ** 4 - ( T_amb + 273.15 ) ** 4 )

def clock(f, hi, lo, t, duty):
     return lo-(np.heaviside(duty-0.5+(1 / np.pi)*np.arcsin(np.sin(np.pi * t * f)), 1))*(lo-hi)

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
    U_A = clock(61, 50, 0, t, 0.5)
    I_L, U_int, T_i, T_c, T_A, V_c = x[0], x[1], x[2], x[3], x[4], x[5]
    R_s = ( 1 + alpha * ( T_i - T_r ) ) * ( 2 * R_Al + 1.1 * R_EDRY * ( R_EWET / R_EDRY ) ** ( V_c ** 2 / V_s ** 2 ))
    I_L = ( U_A - U_int - R_s * I_L) / L_s
    I_R = 1e-12 * ( safeexp( - U_int ) - 1 )
    I_B = I_breakdown / U_beta ** U_breakdown * ( U_beta ** ( U_int) - U_beta ** ( - U_int ))
    U_int = ( I_L + ( - U_int ) / R_P + I_R + I_B) / C
    T_i = 0
    T_c = 0
    T_A = 0
    V_c = 0
    
    rec = U_int
    return np.array([I_L, U_int, T_i, T_c, T_A, V_c, rec])

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
        Rec[x1 - burnsteps - 1] = State[ 6 ]
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