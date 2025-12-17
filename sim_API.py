import matplotlib.pyplot as plt
import numpy as np
import scipy as sci
import scipy.optimize as op
import random as rng
import string
import sys
import time
import hashlib
from cryptography.fernet import Fernet

def c_encrypt(key, massage):
    f = Fernet(key)
    return f.encrypt(massage.encode("utf-8"))
    
def c_decrypt(key, digest):
    f = Fernet(key)
    return f.decrypt(digest).decode("utf-8")

def readfile(name):
    with open(name) as f: tmp = f.read()
    return tmp

def imp():
    print("This feature will be implemented soon.")

def MD5(strg):
    return hashlib.md5(strg.encode("utf8")).hexdigest()

def SHA256(strg):
    return hashlib.sha256(strg.encode("utf8")).hexdigest()

def clean(text):
    allowed = string.ascii_letters + string.digits + "_-."
    return "" if not all ( ch in allowed for ch in text ) else text

def rngstr(length):
    chars = string.ascii_letters * 2 + string.digits + "+-*/=()!%&?#_;:.,$"
    return "".join( rng.choice( chars ) for _ in range( length ) )

def UI(DeltaTime, Burntime, Timeframe, Enable_console, Confirm_num_len, scope):
    lb()
    print("Done.")
    lb()
    INP = "Y"
    
    imp_s = ["implicit Euler", "Gauss Legendre Runge Kutta 2", "Gauss Legendre Runge Kutta 4", "Gauss Legendre Runge Kutta 6"]
    
    burnsteps =  int( Burntime / DeltaTime)
    steps = int(( Timeframe )/DeltaTime)
    
    y = True
    while y == True:
        INP = input(f"Confirm simulating {steps + burnsteps} samples? [Y]/[N]/[i]").casefold().strip()
        if INP in [ "y", "n" ]:
            y = False
            break
        if INP == "i":
            exec("""def query(): return ODEsolver.replace("eE", "explicit_Euler").replace("iE", "implicit_Euler").replace("GLRK", "Gauss_Legendre_Runge_Kutta_").replace("RK", "Runge_Kutta_").replace("_", " "), Timeframe, DeltaTime, Burntime, EQsolver.replace("cN", "custom_Newton").replace("fS", "fSolve").replace("_", " "), Save_Data, Save_Format, Save_Filename, SuppHash, FileHash""" ,scope)
            ODEsolver, Timeframe, DeltaTime, Burntime, EQsolver, Save_Data, Save_Format, Save_Filename, SuppHash, FileHash = scope["query"]()
            steps, burnsteps = int( Timeframe / DeltaTime ), int( Burntime  / DeltaTime )
            info = f"""

The Timestep in the simulation is set to {DeltaTime} seconds.
Simulating {steps + burnsteps} samples. {burnsteps} samples will be discarded ({Burntime} seconds), {steps} samples will be recorded ({Timeframe} seconds).
Using {ODEsolver}{ f" with {EQsolver}" if ODEsolver in imp_s else "" }.
The resulting data will {"" if Save_Data == True else "not"}be saved to disk {f" in a {Save_Format} file named {Save_Filename}{Save_Format}" if Save_Data == True else ""}.
The SHA256 of the current file is                {FileHash} .
The SHA256 of the current file is supposed to be {SuppHash} .
{"The current hash does NOT match the supposed hash. This indicates that the file has been modified since the last update of the supposed hash." if FileHash != SuppHash else ""}
"""
            print(info)
        elif ( INP == "console" or INP == "con" ) and Enable_console == True:
            RNGVAR = rngstr(Confirm_num_len)
            lb()
            print("Warning: usage of this function may break the softwear!")
            cinp = input(f"To confirm entering the console please enter the following key: \n{RNGVAR} \n")
            if cinp == RNGVAR:
                lb()
                print("Console:") # The console is so restrictive, that its safe.
                lb()
                while True:
                    cinp = input(">>> ").casefold().strip()
                    lb()
                    if cinp == "exit":
                        print("Exited console.")
                        lb()
                        break
                    elif cinp == "set":
                        var = input("Enter the name of the variable you want to set: ")
                        val = input("Enter the value you want to set the variable to: ")
                        lb()
                        val = clean(val)
                        var = clean(var)
                        exec(f"{var} = {val}", scope)
                        print(f"Set {var} to {val}.")
                        lb()
                    elif cinp == "inspect" or cinp == "ins":
                        var = input("Enter the name of the variable you want to inspect: ")
                        lb()
                        exec(f"print({var})", scope)
                        lb()
                    elif cinp == "stringset" or cinp == "strset":
                        var = input("Enter the name of the variable you want to set: ")
                        val = input("Enter the string you want to set the variable to: ")
                        lb()
                        val = clean(val)
                        var = clean(var)
                        exec(f"{var} = '{val}'", scope)
                        print(f"Set {var} to {val}.")
                        lb()
                    elif cinp == "help":
                        lb()
                        print("Help menu:")
                        print("Available commands: set, string set, ins, exit")
                    else:
                        print("Invalid command!")
                        lb()
            else:
                lb()
                print("Incorrect numer!")
                lb()
                pass
        elif INP == "" :
            pass
        else:
            print("Wrong input, please try again.")

    if INP == "y":
        pass
    else: 
        exit()

def newton_solve(F, x0, tol=1e-9, max_iter=20):
    x = x0.astype(float).copy()
    n = len(x)
    eps = 1e-8
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
    eps = 1e-8
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
    sigma = 5.670374419e-8
    return A * emiss * sigma * ( ( T + 273.15 ) ** 4 - ( T_amb + 273.15 ) ** 4 )

def clock(t, type, hi=1, lo=-1, f=1, duty=0.5, phase=0):
    if type == "square" or type == "sq":
        return hi if ( ( f * t + phase ) % 1 ) > duty else lo
    if type == "sine" or type == "sn":
        return linint( hi, lo, 0.5 + 0.5 * np.sin( 2 * np.pi * ( f * t + phase ) ) )
    if type == "triangle" or type == "tg":
        return linint( hi, lo, 2 * np.abs( ( f * t + phase ) % 1 - 0.5 ) )
    if type == "saw" or type == "sw":
        return linint( hi, lo, ( f * t + phase ) % 1 )
    if type == "white_noise" or type == "wn":
        return linint( hi, lo, rng.random() )
    if type == "0":
        return linint( hi, lo, 0 )
    else:
        return 0

def progress(percent):
    percent = np.clip(percent, a_min=0, a_max=100)
    bar_len = 50
    filled = int(bar_len * percent / 100)
    bar = "█" * filled + " " * (bar_len - filled)
    sys.stdout.write(f"\r|{bar}| {percent:.0f}%")
    sys.stdout.flush()

def safeexp(x):
     return np.exp(np.clip(x, a_max=700, a_min=-750))

def lb():
    print("")

def set_const(SuppHash, sim_name="sim.py", config_name="config.ini"):
    config = readfile(f"{config_name}")
    FC = f"\n\n#{sim_name}: \n\n" + readfile(f"{sim_name}").replace(SuppHash, "") + "\n\n#sim_API.py: \n\n" + readfile("sim_API.py") + f"\n\n#{config_name}: \n\n" + config
    FileHash = SHA256(FC)
    StartTime = time.perf_counter()
    return config, FC, FileHash, StartTime
 
def save_file(Save_Data, Save_Format, Save_Filename, D):
    if Save_Format == ".npz" and Save_Data == True:
        np.savez(f"{Save_Filename}.npz", D)
    if Save_Format == ".txt" and Save_Data == True:
        np.savetxt(f"{Save_Filename}.txt", D)
    if Save_Format == ".csv" and Save_Data == True:
        np.savetxt(f"{Save_Filename}.csv", D, delimiter=",")

def eqsolve(F, G0, EQsolver):
        if EQsolver == "custom_Newton" or EQsolver == "cN":
            return newton_solve(F, G0)
        if EQsolver == "fSolve" or EQsolver == "fS":
            return op.fsolve(F, G0)
        if EQsolver == "hybrid" or EQsolver == "hybr":
            return op.root(F, G0, method="hybr")
 
def step(df, t, dState, State, ODEsolver, DeltaTime, EQsolver):
    #solve for new state
    G0  = dState
    dSl = len(dState)
    
    if ODEsolver == "explicit Euler":
        return dState + DeltaTime * df(t, dState, State)

    if ODEsolver == "implicit Euler":
        def F(xn):
            return xn - dState - DeltaTime * df(t, xn, State)
        return eqsolve(F, G0, EQsolver) 
    
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
        return dState + DeltaTime * eqsolve(F, G0, EQsolver)

    if ODEsolver == "Gauss Legendre Runge Kutta 4": 
        G0 = np.concatenate((G0,G0))
        def F(k):
            k1 = k[:dSl]
            k2 = k[dSl:]
            F1 = k1 - df(t + ( 0.5 - ( 3 ** 0.5) / 6 ) * DeltaTime, dState + DeltaTime * ( ( 1 / 4                  ) * k1 + ( 1 / 4 - ( 3 ** 0.5) / 6 ) * k2 ), State)
            F2 = k2 - df(t + ( 0.5 + ( 3 ** 0.5) / 6 ) * DeltaTime, dState + DeltaTime * ( ( 1 / 4 + ( 3 ** 0.5) / 6) * k1 + ( 1 / 4                   ) * k2 ), State)
            return np.concatenate((F1,F2))
        k = eqsolve(F, G0, EQsolver)
        k1 = k[:dSl]
        k2 = k[dSl:]
        return dState + DeltaTime / 2 * ( k1 + k2 )
    
    if ODEsolver == "Gauss Legendre Runge Kutta 6":
        G0 = np.concatenate((G0,G0,G0))
        def F(k):
            k1 = k[:dSl]
            k2 = k[dSl:2*dSl]
            k3 = k[2*dSl:]
            
            F1 = k1 - df(t + ( 1 / 2 - ( 15 ** 0.5 ) / 10 ) * DeltaTime, dState + DeltaTime * ( ( 5 / 36                      ) * k1 + ( 2 / 9 - ( 15 ** 0.5 ) / 15 ) * k2 + ( 5 / 36 - ( 15 ** 0.5 ) / 30 ) * k3 ), State)
            F2 = k2 - df(t + ( 1 / 2                      ) * DeltaTime, dState + DeltaTime * ( ( 5 / 36 - ( 15 ** 0.5 ) / 24 ) * k1 + ( 2 / 9                      ) * k2 + ( 5 / 36 - ( 15 ** 0.5 ) / 24 ) * k3 ), State)
            F3 = k3 - df(t + ( 1 / 2 - ( 15 ** 0.5 ) / 10 ) * DeltaTime, dState + DeltaTime * ( ( 5 / 36 + ( 15 ** 0.5 ) / 30 ) * k1 + ( 2 / 9 + ( 15 ** 0.5 ) / 15 ) * k2 + ( 5 / 36                      ) * k3 ), State)
            
            return np.concatenate((F1,F2,F3))
        k = eqsolve(F, G0, EQsolver)
        
        k1 = k[:dSl]
        k2 = k[dSl:2*dSl]
        k3 = k[2*dSl:]
        
        return dState + DeltaTime * ( 5 / 18 * k1 + 4 / 9 * k2 + 5 / 18 * k3 )

def plot(Plot, t, s, xlab="Time in s", ylab="Y-axis", dlab="Diagram"):
    if Plot == "Graph":
        fig, ax = plt.subplots()
        ax.plot(t, s)

        ax.set(xlabel=xlab, ylabel=ylab,
            title=dlab)
        #ax.set_yscale("log")
        plt.tick_params(axis="both", which="both")
        ax.grid()

        fig.savefig("test.png")
        print("Done.")
        plt.show()

    elif Plot == "Animation":
        imp()
    
    input("The end of the programm was reached. Press enter to exit.")

def run_sim(DeltaTime, State, dState, Timeframe, Burntime, f, df, ODEsolver, EQsolver, Show_bar = True, Enable_UI = True):
    
    if Enable_UI == True:
        print("0/1 Complete.")
    
    ODEsolver = ODEsolver.replace("eE", "explicit_Euler").replace("iE", "implicit_Euler").replace("GLRK", "Gauss_Legendre_Runge_Kutta_").replace("RK", "Runge_Kutta_").replace("_", " ")
    
    steps = int(( Timeframe )/DeltaTime)
    burnsteps =  int( Burntime / DeltaTime)
    
    Time = 0
    
    TIME = np.linspace(Burntime, Burntime + Timeframe, steps)
    Rec = np.zeros( (len(State) + len(dState) + 1, steps), dtype = np.float64)
    
    t_arr = np.array([0])
    
    progress_bar_update_time = max( 1, int( ( steps + burnsteps ) / 100 ) )
    
    for x1 in range( steps + burnsteps ):
        Time += DeltaTime #keeping time
        #--- stuff VVV
        State = f(Time, dState, State)
        dState = step(df, Time, dState, State, ODEsolver, DeltaTime, EQsolver)
        
        if x1 >= burnsteps: #recording data
            t_arr[0] = Time
            Sim_State = np.concatenate( (State, dState, t_arr) )
            Rec[:, x1 - burnsteps] = Sim_State
        if x1 % progress_bar_update_time == 0 and Show_bar == True:
            progress(100 * ( x1 + 1 ) / ( steps + burnsteps))
    ODEsolver = ODEsolver.replace(" ", "_")
    if ODEsolver in ["explicit Euler", "implicit_Euler"]:
        E_exp = 1
    elif ODEsolver in ["Gauss_Legendre_Runge_Kutta_2", "Runge_Kutta_2"]:
        E_exp = 2
    elif ODEsolver in ["Gauss_Legendre_Runge_Kutta_4", "Runge_Kutta_4"]:
        E_exp = 4
    elif ODEsolver in ["Gauss_Legendre_Runge_Kutta_6", "Runge_Kutta_6"]:
        E_exp = 6
    Error = f"The global simulation error is on the order of O({DeltaTime ** E_exp:.2e})"
    
    if Enable_UI == True:
        lb()
        print("1/1 Complete.")
        print("Please wait. . .")
    return Rec, Error

def integrate(f, a, b, s): # \int_{a}^{b}f(x)dx
    dx = ( b - a ) / s
    x = a
    c = 0
    fc = f( a )
    fn = f( a + dx )
    for _ in range(s):
        c += dx * 0.5 * ( fc + fn )
        x += dx
        fc = fn
        if _ < s - 1:
            fn = f( x + dx )
    return c
 
def drichelt_function(x):
    return 1 # this is a joke,because computers can onlyrepresent and store rationals

def no_f(t, x, s):
    return s
    
"""
Recommended Project configuration:
Project (Folder)
│
├──sim.py
├──sim_API.py
└──config.ini
"""
