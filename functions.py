import matplotlib.pyplot as plt
import numpy as np
import scipy as sci
import scipy.optimize as op
import random as rng
import string
import sys
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

def clean(text):
    text = text.replace("UI()", "print(\"This function is not available.\")").replace("FILECONTENT", "placeholder_var")
    text = text.replace("FileHash", "placeholder_var").replace("config", "placeholder_var").replace("session_key", "placeholder_var")
    
    
    return text

def rngstr(length):
    chars = string.ascii_letters * 2 + string.digits + "+-*/=()!%&?#_;:.,$"
    return "".join( rng.choice( chars ) for _ in range( length ) )

def UI(steps, burnsteps, ODEsolver, EQsolver, Save_Data, Save_Format, Save_Filename, DeltaTime, Burntime, Timeframe, FileHash, SuppHash, Enable_console, Confirm_num_len):
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
            print(f"The MD5 of the current file is                {FileHash} .")
            print(f"The MD5 of the current file is supposed to be {SuppHash} .")
            if FileHash != SuppHash:
                print(f"The current hash does NOT match the supposed hash. This indicates that the file has been modified sonce the last update of the supposed hash.")
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
                    cinp = input(">>> ").replace("UI()", "print(\"This function is not available.\")") 
                    cinp = clean(cinp)
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

def mini_UI():
    progress(100)
    lb()
    print("1/1 Complete.")
    print("Please wait. . .")

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

def lb():
    print("")

def set_const():
    consts = """
sigma = 5.670374419e-8 #W/m^2*K^4
Grav = 1#6.6743e-11
pi = np.pi
###
eps = 1e-8
eps2 = 1e-1
placeholder_var = 0
    """
    exec(consts, globals())
 
def save_file(Save_Data, Save_Format, Save_Filename, D):
    if Save_Format == ".npz" and Save_Data:
        np.savez(f"{Save_Filename}.npz", D)
    if Save_Format == ".txt" and Save_Data:
        np.savetxt(f"{Save_Filename}.txt", D)
    if Save_Format == ".csv" and Save_Data:
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

def plot(Plot, t, s):
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