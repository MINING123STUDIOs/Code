print("Setting up simulation. . . (advanced simulation model 2025)")
def lb():
    print("")

lb()
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
#from scipy.special import lambertw
from matplotlib.animation import FuncAnimation

#setting constants 
for _ in range(1):
    sigma = 5.670374419e-8 #W/m^2*K^4
    Grav = 1#6.6743e-11
    pi = np.pi
    ###
    GL4_c1 = 0.5 - ( 3 ** 0.5) / 6
    GL4_c2 = 0.5 + ( 3 ** 0.5) / 6
    GL4_a11 = GL4_a22 = 0.25
    GL4_a12 = 0.25 - ( 3 ** 0.5) / 6
    GL4_a21 = 0.25 + ( 3 ** 0.5) / 6
    GL4_b1 = GL4_b2 = 0.5
    ###
    GL6_c1 = 0.5 - ( 15 ** 0.5 ) / 10
    GL6_c2 = 0.5
    GL6_c3 = 0.5 + ( 15 ** 0.5 ) / 10
    
    GL6_b1 = 5 / 18 
    GL6_b2 = 4 / 9 
    GL6_b3 = 5 / 18 
    
    GL6_a11 = 5 / 36 
    GL6_a12 = 2 / 9  - ( 15 ** 0.5 ) / 15
    GL6_a13 = 5 / 36 - ( 15 ** 0.5 ) / 30
    GL6_a21 = 5 / 36 + ( 15 ** 0.5 ) / 24
    GL6_a22 = 2 / 9 
    GL6_a23 = 5 / 36 - ( 15 ** 0.5 ) / 24
    GL6_a31 = 5 / 36 + ( 15 ** 0.5 ) / 30
    GL6_a32 = 2 / 9  + ( 15 ** 0.5 ) / 15
    GL6_a33 = 5 / 36 
    ###
    eps = 1e-8
    eps2 = 1e-1

#sim params
Timeframe = 5 #s
Burntime = 0 #s
DeltaTime = 2e-3 #s
Num_Damp = 1
ODEsolver = "GLRK2" # implicit_Euler / iE / explicit_Euler / eE / Runge_Kutta_4 / RK4 / Gauss_Legendre_Runge_Kutta_4 / GL4
EQsolver = "fS" #custom_Newton / cN / fSolve / fS
Save_Data = False
Save_Format = ".csv" # .csv, .txt, .npz
Save_Filename = "Recording"
Enable_console = True
Confirm_num_len = 8
Plot = "Graph"
#---
m1 = 1
m2 = 1

l1 = 1
l2 = 1

g = 9.81

#current example: double pendulum

def imp():
    print("This feature will be implemented soon.")

def MD5(strg):
    return hashlib.md5(strg.encode("utf8")).hexdigest()
    
def readfile(name):
    with open(name) as f: tmp = f.read()
    return tmp

#calculating secondary params
burnsteps =  -int(- Burntime / DeltaTime) #number of required timesteps to burn
steps = -int(-( Timeframe )/DeltaTime) #number of required timesteps
progress_bar_update_time = max( 1, int( ( steps + burnsteps ) / 100 ) )
ODEsolver = ODEsolver.replace("eE", "explicit_Euler").replace("iE", "implicit_Euler").replace("GLRK", "Gauss_Legendre_Runge_Kutta_").replace("RK", "Runge_Kutta_").replace("_", " ")

#---  
FILECONTENT = readfile("sim.py").replace("ec1d455a9d09faa4eb20722255eba31e", "")
config = readfile("config.ini")
FileHash = MD5(FILECONTENT)
#print(FILECONTENT)

#working variables
TIME = np.linspace(Burntime, Burntime + Timeframe, steps)
Rec = np.zeros( steps )
#---
Time = 0
dState = np.array( [0.2, 0, 0, 0], dtype = np.float64 ) # (Theta1, Theta2, w1, w2) state depending on ODEs 
State = np.array( [ 0 ], dtype = np.float64 ) # (Rec) state not depending on ODEs 

dSl = len(dState)

def UI():
    print("Done.")
    lb()
    INP = "Y"
    y = True
    while y == True:
        INP = input(f"Confirm simulating {steps + burnsteps} samples? [Y]/[N]/[i]").casefold().strip()
        if INP in [ "y", "n" ]:
            y = False
            break
        if INP == "i":
            lb()
            intvar0001 = ODEsolver
            intvar0002 = EQsolver.replace("cN", "custom_Newton").replace("fS", "fSolve").replace("_", " ")
            intvar0002a = f" with {intvar0002}" if ODEsolver in ["implicit Euler", "Gauss Legendre Runge Kutta 2", "Gauss Legendre Runge Kutta 4", "Gauss Legendre Runge Kutta 6"] else ""
            intvar0003 = "" if Save_Data == True else "not"
            intvar0004 = f" in a {Save_Format} file named {Save_Filename}{Save_Format}" if Save_Data == True else ""
            print(f"The Timestep in the simulation is set to {DeltaTime} seconds.")
            print(f"Simulating {steps + burnsteps} samples. {burnsteps} samples will be discarded ({Burntime} seconds), {steps} samples will be recorded ({Timeframe} seconds).")
            print(f"Using {intvar0001}{intvar0002a}.")
            print(f"The resulting data will {intvar0003} be saved to disk{intvar0004}.")
            print(f"The MD5 of the current file is                {FileHash}.")
            print(f"The MD5 of the current file is supposed to be ec1d455a9d09faa4eb20722255eba31e.")
            lb()
        elif INP == "console" and Enable_console == True:
            #RNGVAR = f"{rng.randint( 10 ** Confirm_num_len / 10, 10 ** Confirm_num_len - 1 )}"
            RNGVAR = rngstr(Confirm_num_len)
            lb()
            print("Warning: usage of this function may break the softwear! Usage may also pose a security risk due to the execution of bad code!")
            cinp = input(f"To confirm entering the console please enter the following key: \n{RNGVAR} \n")
            if cinp == RNGVAR:
                lb()
                print("Console:")
                lb()
                while True:
                    cinp = input(">>>").replace("UI()", "print(\"This function is not available.\")") #
                    if cinp == "exit":
                        print("Exited console.")
                        break
                    else:
                        exec(cinp, globals()) #sketchy, but disableable
                        lb()
            else:
                print("Incorrect numer!")
                pass
        else:
            print("Wrong input, please try again.")

    if INP == "y":
        pass
    else: 
        exit()
    print("0/1 Complete.")

def rngstr(length):
    chars = string.ascii_letters * 2 + string.digits + "+-*/=()!%&?#_;:.,$"
    return "".join( rng.choice( chars ) for _ in range( length ) )

def newton_solve(F, x0, tol=1e-9, max_iter=20):
    x = x0.astype(float).copy()
    n = len(x)
    for _ in range(max_iter):
        Fx = F(x)
        if np.linalg.norm(Fx, ord=2) < tol:
            return x
        J = np.zeros((n, n), dtype=float) #needs to be numerical for universality.
        for i in range(n):
            x_pert = x.copy()
            x_pert[i] += eps
            J[:, i] = (F(x_pert) - Fx) / eps
        delta = np.linalg.solve(J, -Fx)
        x += delta
        if np.linalg.norm(delta, ord=2) < tol:
            return x
    return x

def scalar_solve(F, x0, tol=1e-9, max_iter=20):
    x = x0
    for i in range(max_iter):
        x = x - eps * F(x) / ( F( x + eps ) - F( x ) )
        if abs(F(x)) < tol:
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

def eqsolve(F, G0):
        if EQsolver == "custom_Newton" or EQsolver == "cN":
            return newton_solve(F, G0)
        if EQsolver == "fSolve" or EQsolver == "fS":
            return op.fsolve(F, G0)
        if EQsolver == "hybrid" or EQsolver == "hybr":
            return op.root(F, G0, method="hybr")

def step(t, dState, State):
    #solve for new state
    G0 = dState
    
    if ODEsolver == "explicit Euler":
        return dState + DeltaTime * df(t, dState, State)

    if ODEsolver == "implicit Euler":
        def F(xn):
            return xn - dState - DeltaTime * df(t, xn, State)
        return eqsolve(F, G0) 
    
    if ODEsolver == "Runge Kutta 2":
        k1 = df(t, dState, State)
        k2 = df(t + DeltaTime / 2, dState + DeltaTime / 2 * k1, State)
        return dState + DeltaTime * k2
    
    if ODEsolver == "Runge Kutta 4":
        k1 = df(t, dState, State)
        k2 = df(t + DeltaTime / 2, dState + DeltaTime * k1 / 2, State)
        k3 = df(t + DeltaTime / 2, dState + DeltaTime * k2 / 2, State)
        k4 = df(t + DeltaTime, dState + DeltaTime * k3, State)
        return dState + ( DeltaTime / 6 ) * ( k1 + 2 * k2 + 2 * k3 + k4 )
    
    if ODEsolver == "Runge Kutta 6":
         k1 = df(t, dState, State)
         k2 = df(t + DeltaTime / 3, dState + DeltaTime / 3 * k1, State)
         k3 = df(t + DeltaTime / 3, dState + DeltaTime / 6 * ( k1 + k2 ), State)
         k4 = df(t + DeltaTime, dState + DeltaTime * ( k1 + k2 + k3 ), State)
         k5 = df(t + DeltaTime, dState + DeltaTime / 2 * ( k1 + k4 ), State)
         k6 = df(t + DeltaTime / 2, dState + DeltaTime / 8 * ( -3 * k1 + 9 * k4 ), State)
         k7 = df(t + DeltaTime, dState + DeltaTime * ( k1 / 2 - 3 * k3 / 2 + 2* k4 ), State)
         return dState + DeltaTime * ( k1 / 12 + k3 / 4 + k5 / 3 + k7 / 4 )
    
    if ODEsolver == "Gauss Legendre Runge Kutta 2":
        def F(k):
            return k - df(t + DeltaTime / 2, dState + DeltaTime / 2 * k, State)
        return dState + DeltaTime * eqsolve(F, G0)

    if ODEsolver == "Gauss Legendre Runge Kutta 4": 
        G0 = np.concatenate((G0,G0))
        def F(k):
            k1 = k[:dSl]
            k2 = k[dSl:]
            F1 = k1 - df(t + GL4_c1 * DeltaTime, dState + DeltaTime * ( GL4_a11 * k1 + GL4_a12 * k2 ), State)
            F2 = k2 - df(t + GL4_c2 * DeltaTime, dState + DeltaTime * ( GL4_a21 * k1 + GL4_a22 * k2 ), State)
            return np.concatenate((F1,F2))
        k = eqsolve(F, G0)
        k1 = k[:dSl]
        k2 = k[dSl:]
        return dState + DeltaTime / 2 * ( k1 + k2 )
    
    if ODEsolver == "Gauss Legendre Runge Kutta 6":
        G0 = np.concatenate((G0,G0,G0))
        def F(k):
            k1 = k[:dSl]
            k2 = k[dSl:2*dSl]
            k3 = k[2*dSl:]
            
            F1 = k1 - df(t + GL6_c1 * DeltaTime, dState + DeltaTime * ( GL6_a11 * k1 + GL6_a12 * k2 + GL6_a13 * k3 ), State)
            F2 = k2 - df(t + GL6_c2 * DeltaTime, dState + DeltaTime * ( GL6_a21 * k1 + GL6_a22 * k2 + GL6_a23 * k3 ), State)
            F3 = k3 - df(t + GL6_c3 * DeltaTime, dState + DeltaTime * ( GL6_a31 * k1 + GL6_a32 * k2 + GL6_a33 * k3 ), State)
            
            return np.concatenate((F1,F2,F3))
        k = eqsolve(F, G0)
        
        k1 = k[:dSl]
        k2 = k[dSl:2*dSl]
        k3 = k[2*dSl:]
        
        return dState + DeltaTime * ( GL6_b1 * k1 + GL6_b2 * k2 + GL6_b3 * k3 )

UI()

for x1 in range( steps + burnsteps ):
    Time += DeltaTime #keeping time
    #--- stuff VVV
    State = f(Time, dState, State)
    dState = step(Time, dState, State)
    if x1 >= burnsteps: #recording data
        Rec[x1 - burnsteps ] = dState[0]
        TIME[x1 - burnsteps ] = dState[1]
    if x1 % progress_bar_update_time == 0:
        progress(100 * x1 / ( steps + burnsteps))
progress(100)
lb()
print("1/1 Complete.")
print("Please wait. . .")
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

if Save_Data == True:
    D = np.array([t,s])
    if Save_Format == ".npz":
            np.savez(f"{Save_Filename}.npz", D)
    if Save_Format == ".txt":
        np.savetxt(f"{Save_Filename}.txt", D)
    if Save_Format == ".csv":
        np.savetxt(f"{Save_Filename}.csv", D, delimiter=",")

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

if Plot == "Graph":
    fig, ax = plt.subplots()
    ax.plot(t, s)

    ax.set(xlabel="Time in s", ylabel="Y-axis",
        title="Diagram")
    #ax.set_yscale("log")
    plt.tick_params(axis="both", which="both")
    ax.grid()

    fig.savefig("test.png")
    print("Done.")
    plt.show()

elif Plot == "Animation":
    imp()
    
    """NOTES:
        Note 1:  modifies input array instead of making a new one to improve performance. Due to this being at the end and working one the local copy of dState it does NOT mutate the simulation.
        ⚠
    """
