"""
Contents:
    imports
    string & print functions
    numerical analysis functions
    misl. helper functions
    main simulation functions
    UI functions
    Notes
"""

#imports:
import matplotlib.pyplot as plt
import numpy as np
import scipy as sci
import scipy.optimize as op
import random as rng
import string
import sys
import time
from            math_API import *
from         special_API import *

#string & print functions:

def imp():
    print("This feature will be implemented soon.")

def lb():
    print("")

def progress(percent):
    percent = np.clip(percent, a_min=0, a_max=100)
    bar_len = 50
    filled = int(bar_len * percent / 100)
    bar = "█" * filled + " " * (bar_len - filled)
    sys.stdout.write(f"\r|{bar}| {percent:.0f}%")
    sys.stdout.flush()

#numerical analysis functions:

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

def eqsolve(F, G0, eq_solver):
        if eq_solver == "custom_Newton" or eq_solver == "cN":
            return newton_solve(F, G0)
        if eq_solver == "fSolve" or eq_solver == "fS":
            return op.fsolve(F, G0)
        if eq_solver == "hybrid" or eq_solver == "hybr":
            return op.root(F, G0, method="hybr")
 
def step(df, t, d_state, state, ode_solver, delta_time, eq_solver):
    #solve for new state
    g0  = df(t, d_state, state)
    dsl = len(d_state)
    
    if ode_solver == "explicit Euler":
        return d_state + delta_time * df(t, d_state, state)

    if ode_solver == "implicit Euler":
        def F(xn):
            return xn - d_state - delta_time * df(t, xn, state)
        return eqsolve(F, g0, eq_solver) 
    
    if ode_solver == "Runge Kutta 2":
        k1 = df(t, d_state, state)
        k2 = df(t + delta_time / 2, d_state + delta_time / 2 * k1, state)
        return d_state + delta_time * k2
    
    if ode_solver == "Runge Kutta 4":
        k1 = df(t, d_state, state)
        k2 = df(t + delta_time / 2, d_state + delta_time * k1 / 2, state)
        k3 = df(t + delta_time / 2, d_state + delta_time * k2 / 2, state)
        k4 = df(t + delta_time, d_state + delta_time * k3, state)
        return d_state + ( delta_time / 6 ) * ( k1 + 2 * k2 + 2 * k3 + k4 )
    
    if ode_solver == "Runge Kutta 6":
         k1 = df(t, d_state, state)
         k2 = df(t + delta_time / 3, d_state + delta_time / 3 * k1, state)
         k3 = df(t + delta_time / 3, d_state + delta_time / 6 * ( k1 + k2 ), state)
         k4 = df(t + delta_time, d_state + delta_time * ( k1 + k2 + k3 ), state)
         k5 = df(t + delta_time, d_state + delta_time / 2 * ( k1 + k4 ), state)
         k6 = df(t + delta_time / 2, d_state + delta_time / 8 * ( -3 * k1 + 9 * k4 ), state)
         k7 = df(t + delta_time, d_state + delta_time * ( k1 / 2 - 3 * k3 / 2 + 2* k4 ), state)
         return d_state + delta_time * ( k1 / 12 + k3 / 4 + k5 / 3 + k7 / 4 )
    
    if ode_solver == "Gauss Legendre Runge Kutta 2":
        def F(k):
            return k - df(t + delta_time / 2, d_state + delta_time / 2 * k, state)
        return d_state + delta_time * eqsolve(F, g0, eq_solver)

    if ode_solver == "Gauss Legendre Runge Kutta 4": 
        g0 = np.concatenate((g0,g0))
        def F(k):
            k1 = k[:dsl]
            k2 = k[dsl:]
            F1 = k1 - df(t + ( 0.5 - ( 3 ** 0.5) / 6 ) * delta_time, d_state + delta_time * ( ( 1 / 4                  ) * k1 + ( 1 / 4 - ( 3 ** 0.5) / 6 ) * k2 ), state)
            F2 = k2 - df(t + ( 0.5 + ( 3 ** 0.5) / 6 ) * delta_time, d_state + delta_time * ( ( 1 / 4 + ( 3 ** 0.5) / 6) * k1 + ( 1 / 4                   ) * k2 ), state)
            return np.concatenate((F1,F2))
        k = eqsolve(F, g0, eq_solver)
        k1 = k[:dsl]
        k2 = k[dsl:]
        return d_state + delta_time / 2 * ( k1 + k2 )
    
    if ode_solver == "Gauss Legendre Runge Kutta 6":
        g0 = np.concatenate((g0,g0,g0))
        def F(k):
            k1 = k[:dsl]
            k2 = k[dsl:2*dsl]
            k3 = k[2*dsl:]
            
            F1 = k1 - df(t + ( 1 / 2 - ( 15 ** 0.5 ) / 10 ) * delta_time, d_state + delta_time * ( ( 5 / 36                      ) * k1 + ( 2 / 9 - ( 15 ** 0.5 ) / 15 ) * k2 + ( 5 / 36 - ( 15 ** 0.5 ) / 30 ) * k3 ), state)
            F2 = k2 - df(t + ( 1 / 2                      ) * delta_time, d_state + delta_time * ( ( 5 / 36 - ( 15 ** 0.5 ) / 24 ) * k1 + ( 2 / 9                      ) * k2 + ( 5 / 36 - ( 15 ** 0.5 ) / 24 ) * k3 ), state)
            F3 = k3 - df(t + ( 1 / 2 + ( 15 ** 0.5 ) / 10 ) * delta_time, d_state + delta_time * ( ( 5 / 36 + ( 15 ** 0.5 ) / 30 ) * k1 + ( 2 / 9 + ( 15 ** 0.5 ) / 15 ) * k2 + ( 5 / 36                      ) * k3 ), state)
            
            return np.concatenate((F1,F2,F3))
        k = eqsolve(F, g0, eq_solver)
        
        k1 = k[:dsl]
        k2 = k[dsl:2*dsl]
        k3 = k[2*dsl:]
        
        return d_state + delta_time * ( 5 / 18 * k1 + 4 / 9 * k2 + 5 / 18 * k3 )

#misl. helper functions:

def save_file(save_data, save_format, save_filename, d):
    if save_format == ".npz" and save_data:
        np.savez(f"{save_filename}.npz", d)
    if save_format == ".txt" and save_data:
        np.savetxt(f"{save_filename}.txt", d)
    if save_format == ".csv" and save_data:
        np.savetxt(f"{save_filename}.csv", d, delimiter=",")

def no_f(t, x, s):
    return s

def set_const(supp_hash, sim_name="sim.py", config_name="config.ini"):
    config = readfile(f"{config_name}")
    fc = f"\n\n#{sim_name}: \n\n" + readfile(f"{sim_name}").replace(supp_hash, "") + f"\n\n#{config_name}: \n\n" + config
    file_hash = SHA256(fc)
    start_time = time.perf_counter()
    return config, fc, file_hash, start_time

#main simulation functions:

def run_sim(delta_time, state, d_state, time_frame, burntime, f, df, ode_solver, eq_solver="fS", show_bar = True, enable_ui = True):
    
    if enable_ui:
        print("0/1 Complete.")
    
    ode_solver = norm_solver_name(ode_solver)
    
    steps = int(( time_frame )/delta_time)
    burnsteps =  int( burntime / delta_time)
    
    time = 0
    
    rec = np.zeros( (len(state) + len(d_state) + 1, steps), dtype = np.float64)
    
    t_arr = np.array([0])
    
    progress_bar_update_time = max( 1, int( ( steps + burnsteps ) / 100 ) )
    
    for x1 in range( steps + burnsteps ):
        time += delta_time #keeping time
        #--- stuff VVV
        state = f(time, d_state, state)
        d_state = step(df, time, d_state, state, ode_solver, delta_time, eq_solver)
        
        if x1 >= burnsteps: #recording data
            t_arr[0] = time
            sim_state = np.concatenate( (state, d_state, t_arr) )
            rec[:, x1 - burnsteps] = sim_state
        if x1 % progress_bar_update_time == 0 and show_bar and enable_ui:
            progress(100 * ( x1 + 1 ) / ( steps + burnsteps))
    
    if enable_ui:
        lb()
        print("1/1 Complete.")
        print("Please wait. . .")
    return rec

def plot(plot_type, d, datalab, xlab="Time in s", ylab="Y-axis", dlab="Diagram", logy=False, show_plot=True, ui=True, save_plot=False, save_plot_name="test.png" ):
    if plot_type == "Graph":
        
        w = np.array(np.shape(d))
        
        fig, ax = plt.subplots(layout='constrained')
        
        for m in range( int( w[0] / 2 ) ):
            x, y = d[ 2 * m - 2, : ], d[ 2 * m - 1, : ]
            ax.plot(x, y, label=datalab[int( m - 1 ) % len( datalab ) ])
        
        ax.set(xlabel=xlab, ylabel=ylab,
            title=dlab)
            
        if logy: ax.set_yscale("log")
         
        plt.tick_params(axis="both", which="both")
        ax.grid()
        ax.legend()

        if save_plot: fig.savefig(save_plot_name)
        if ui: print("Done.")
        if show_plot: plt.show()

    elif plot_type == "Animation":
        imp()
    
    if ui: input("The end of the programm was reached. Press enter to exit.")

#UI functions:

def console(scope):
    while True:
        cinp = input(">>> ").casefold().strip()
        lb()
        if   cinp == "exit":
            print("Exited console.")
            lb()
            break
        elif cinp == "set":
            var = input("Enter the name of the variable you want to set: ")
            val = input("Enter the value you want to set the variable to: ")
            lb()
            val = clean(val)
            var = clean(var)
            scope[var] = float(val)
            print(f"Set {var} to {val}.")
            lb()
        elif cinp == "inspect" or cinp == "ins":
            var = input("Enter the name of the variable you want to inspect: ")
            lb()
            if var in scope:
                print(scope[var])
            else:
                print(f"No variable with the name \"{var}\" is defined.")
            lb()
        elif cinp == "stringset" or cinp == "strset":
            var = input("Enter the name of the variable you want to set: ")
            val = input("Enter the string you want to set the variable to: ")
            lb()
            val = clean(val)
            var = clean(var)
            scope[var] = f"{val}"
            print(f"Set {var} to {val}.")
            lb()
        elif cinp == "help":
            lb()
            print("Help menu:")
            print("Available commands: set, string set, ins, exit")
        elif cinp == "kill":
            exit()
        else:
            print("Invalid command!")
            lb()

def info(scope):
    imp_s = ["implicit Euler", "Gauss Legendre Runge Kutta 2", "Gauss Legendre Runge Kutta 4", "Gauss Legendre Runge Kutta 6"]
    ode_solver, eq_solver = norm_solver_name(scope["ode_solver"]), scope["eq_solver"]
    timeframe, delta_time, burntime = scope["timeframe"], scope["delta_time"], scope["burntime"]
    save_data, save_format, save_filename = scope["save_data"], scope["save_format"], scope["save_filename"]
    supp_hash, file_hash = scope["supp_hash"], scope["file_hash"]
    steps, burnsteps = int( timeframe / delta_time ), int( burntime  / delta_time )
    info = f"""

The Timestep in the simulation is set to {delta_time} seconds.
Simulating {steps + burnsteps} samples. {burnsteps} samples will be discarded ({burntime} seconds), {steps} samples will be recorded ({timeframe} seconds).
Using {ode_solver}{ f" with {eq_solver}" if ode_solver in imp_s else ""}.
The resulting data will {"" if save_data else "not "}be saved to disk{f" in a {save_format} file named {save_filename}{save_format}" if save_data else ""}.
The SHA256 of the current file is                {file_hash} .
The SHA256 of the current file is supposed to be {supp_hash} .
{"The current hash does NOT match the supposed hash. This indicates that the file has been modified since the last update of the supposed hash." if file_hash != supp_hash else ""}
"""
    print(info)

def UI(scope, enable_console=False, confirm_num_len=8):
    print("\nDone.\n")
    while True:
        s = int( ( scope["burntime"] + scope["timeframe"] ) / scope["delta_time"] )
        inp = input(f"Confirm simulating {s} samples? [Y]/[N]/[i]").casefold().strip()
        rng_var = rngstr(confirm_num_len)
        
        if   inp == "i": info(scope)
        elif inp == "y": break
        elif inp == "n": exit()
        elif ( inp == "console" or inp == "con" ) and enable_console:
            if input(f"\nWarning: usage of this function may break the softwear!\nTo confirm entering the console please enter the following key: \n{rng_var} \n") == rng_var:
                print("\nConsole:\n") # The console is so restrictive, that its safe
                console(scope)
            else: print("\nIncorrect numer!\n")
        elif inp !=  "": print("Wrong input, please try again.")

#Notes:
"""
    Recommended Project configuration:
    Project (Folder)
    │
    ├──sim.py
    ├──sim_API.py
    └──config.ini
"""
