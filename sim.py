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
from scipy.special import lambertw
import sys

#setting constants 
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
eps = 1e-8

#sim params
Timeframe = 200 #s
Burntime = 0 #s
DeltaTime = 1e-1 #s
Num_Damp = 1
ODEsolver = "GLRK4" # implicit_Euler / iE / explicit_Euler / eE / Runge_Kutta_4 / RK4 / Gauss_Legendre_Runge_Kutta_4 / GL4
EQsolver = "fS" #custom_Newton / cN / fSolve / fS
Save_Data = False
Save_Format = ".csv" # .csv, .txt, .npz
Save_Filename = "Recording"
Enable_console = True
Confirm_num_len = 8
#---
m1 = 1
m2 = 1
m3 = 1
#current example: 3 body problem.

#calculating secondary params
burnsteps =  -int(- Burntime / DeltaTime) #number of required timesteps to burn
steps = -int(-( Timeframe )/DeltaTime) #number of required timesteps
progress_bar_update_time = max( 1, int( ( steps + burnsteps ) / 100 ) )
#---


#working variables
TIME = np.linspace(Burntime, Burntime + Timeframe, steps)
Rec = np.zeros( steps )
#---
Time = 0
dState = np.array( [ -0.97000436, 0.2438753, 0.4662036850, 0.4323657300, 0.97000436, -0.24308753, 0.4662036850, 0.4323657300, 0, 0, -0.93240737, -0.86473146], dtype = np.float64 ) # (x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3) state depending on ODEs 
State = np.array( [ 0 ], dtype = np.float64 ) # (Rec) state not depending on ODEs 

dSl = len(dState)

#figure 8: [ -0.97000436, 0.2438753, 0.4662036850, 0.4323657300, 0.97000436, -0.24308753, 0.4662036850, 0.4323657300, 0, 0, -0.93240737, -0.86473146]


def UI():
    INP = "Y"
    y = True
    while y == True:
        INP = input(f"Confirm simulating {steps + burnsteps} samples? [Y]/[N]/[i]").casefold().strip()
        if INP in [ "y", "n" ]:
            y = False
            break
        if INP == "i":
            lb()
            intvar0001 = ODEsolver.replace("eE", "explicit_Euler").replace("iE", "implicit_Euler").replace("GLRK4", "Gauss_Legendre_Runge_Kutta_4").replace("RK4", "Runge_Kutta_4").replace("_", " ")
            intvar0002 = EQsolver.replace("cN", "custom_Newton").replace("fS", "fSolve").replace("_", " ")
            intvar0002a = f" with {intvar0002}" if ODEsolver in ["iE", "implicit_Euler", "GLRK4", "Gauss_Legendre_Runge_Kutta_4"] else ""
            intvar0003 = "" if Save_Data == True else "not"
            intvar0004 = f" in a {Save_Format} file named {Save_Filename}{Save_Format}" if Save_Data == True else ""
            print(f"Simulating {steps + burnsteps} samples. {burnsteps} samples will be discarded ({Burntime} seconds), {steps} samples will be recorded ({Timeframe} seconds).")
            print(f"Using {intvar0001}{intvar0002a}.")
            print(f"The resulting data will {intvar0003} be saved to disk{intvar0004}.")
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
    x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3 = x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11]
    
    dx1 = vx1
    dy1 = vy1
    dx2 = vx2
    dy2 = vy2
    dx3 = vx3
    dy3 = vy3
    
    r12 = np.sqrt( ( x2 - x1 ) ** 2 + ( y2 - y1 ) ** 2 ) ** 3
    r13 = np.sqrt( ( x3 - x1 ) ** 2 + ( y3 - y1 ) ** 2 ) ** 3
    r23 = np.sqrt( ( x3 - x2 ) ** 2 + ( y3 - y2 ) ** 2 ) ** 3 #numerically fragile, but fine
    
    dvx1 = Grav * ( m2 * ( x2 - x1 ) / r12 + m3 * ( x3 - x1 ) / r13 )
    dvy1 = Grav * ( m2 * ( y2 - y1 ) / r12 + m3 * ( y3 - y1 ) / r13 )
    dvx2 = Grav * ( m1 * ( x1 - x2 ) / r12 + m3 * ( x3 - x2 ) / r23 )
    dvy2 = Grav * ( m1 * ( y1 - y2 ) / r12 + m3 * ( y3 - y2 ) / r23 )
    dvx3 = Grav * ( m1 * ( x1 - x3 ) / r13 + m2 * ( x2 - x3 ) / r23 )
    dvy3 = Grav * ( m1 * ( y1 - y3 ) / r13 + m2 * ( y2 - y3 ) / r23 )
    
    x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11] = dx1, dy1, dvx1, dvy1, dx2, dy2, dvx2, dvy2, dx3, dy3, dvx3, dvy3
    return x # modifies input array instead of making a new one to improve performance

def f(t, x, s):
    x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3 = x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11]
    
    r12 = np.sqrt( ( x2 - x1 ) ** 2 + ( y2 - y1 ) ** 2 )
    r13 = np.sqrt( ( x3 - x1 ) ** 2 + ( y3 - y1 ) ** 2 )
    r23 = np.sqrt( ( x3 - x2 ) ** 2 + ( y3 - y2 ) ** 2 )
    
    E = 0.5 * ( m1 * ( vx1 ** 2 + vy1 ** 2 ) + m2 * ( vx2 ** 2 + vy2 ** 2 ) + m3 * ( vx3 ** 2 + vy3 ** 2 ) ) - Grav * ( m1 * m2 / r12 + m1 * m3 / r13 + m2 * m3 / r23 )
    Rec = E
    s[0] =  Rec
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
    
    if ODEsolver == "implicit_Euler" or ODEsolver == "iE":
        def F(xn):
            return xn - dState - DeltaTime * df(t, xn, State)
        return eqsolve(F, G0) 
    
    if ODEsolver == "explicit_Euler" or ODEsolver == "eE":
        return dState + DeltaTime * df(t, dState, State)
    
    if ODEsolver == "Runge_Kutta_4" or ODEsolver == "RK4":
        k1 = df(t, dState, State)
        k2 = df(t + DeltaTime / 2, dState + DeltaTime * k1 / 2, State)
        k3 = df(t + DeltaTime / 2, dState + DeltaTime * k2 / 2, State)
        k4 = df(t + DeltaTime, dState + DeltaTime * k3, State)
        return dState + ( DeltaTime / 6 ) * ( k1 + 2 * k2 + 2 * k3 + k4 )
    
    if ODEsolver == "Gauss_Legendre_Runge_Kutta_4" or ODEsolver == "GLRK4": 
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
    
    if ODEsolver == "Runge_Kutta_2" or ODEsolver = "RK2":
        k1 = df(t, dState, State)
        k2 = df(t + DeltaTime / 2, dState + DeltaTime / 2 * k1, State)
        return dState + DeltaTime * k2
        
    if ODEsolver == "Gauss_Legendre_Runge_Kutta_2" or ODEsolver == "GLRK2":
        def F(k):
            return k - df(t + DeltaTime / 2, dState + DeltaTime / 2 * k, State)
        return dState + DeltaTime * eqsolve(F, G0)
     
    if ODEsolver == "Runge_Kutta_6" or ODEsolver == "RK6":
         k1 = df(t, dState, State)
         k2 = df(t + DeltaTime / 3, dState + DeltaTime / 3 * k1, State)
         k3 = df(t + DeltaTime / 3, dState + DeltaTime / 6 * ( k1 + k2 ), State)
         k4 = df(t + DeltaTime, dState + DeltaTime * ( k1 + k2 + k3 ), State)
         k5 = df(t + DeltaTime, dState + DeltaTime / 2 * ( k1 + k4 ), State)
         k6 = df(t + DeltaTime / 2, dState + DeltaTime / 8 * ( -3 * k1 + 9 * k4 ), State)
         k7 = df(t + DeltaTime, dState + DeltaTime * ( k1 / 2 - 3 * k3 / 2 + 2* k4 ), State)
         return dState + DeltaTime * ( k1 / 12 + k3 / 4 + k5 / 3 + k7 / 4 )
         
    if ODEsolver == "Gauss_Legendre_Runge_Kutta_6" or ODEsolver == "GLRK6":
        
         return dState

#dState = step(Time, dState, State)
#State = f(Time, dState, State)

print("Done.")
lb()
UI()
print("0/1 Complete.")

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
