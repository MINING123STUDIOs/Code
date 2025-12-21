"""
Contents:
    imports
    functions
    Notes
"""

#imports:
import numpy as np
import scipy.optimize as op
import random as rng
import hashlib
from cryptography.fernet import Fernet
from math_API import *

def c_encrypt(key, massage):
    f = Fernet(key)
    return f.encrypt(massage.encode("utf-8"))
    
def c_decrypt(key, digest):
    f = Fernet(key)
    return f.decrypt(digest).decode("utf-8")

def MD5(strg):
    return hashlib.md5(strg.encode("utf8")).hexdigest()

def SHA256(strg):
    return hashlib.sha256(strg.encode("utf8")).hexdigest()

def clean(text):
    allowed = string.ascii_letters + string.digits + "_-."
    return "" if ( not all ( ch in allowed for ch in text ) ) or text.startswith("__") else text

def rngstr(length):
    chars = string.ascii_letters * 2 + string.digits + "+-*/=()!%&?#_;:.,$"
    return "".join( rng.choice( chars ) for _ in range( length ) )

def norm_solver_name(name): 
    return name.replace("eE", "explicit_Euler").replace("iE", "implicit_Euler").replace("GLRK", "Gauss_Legendre_Runge_Kutta_").replace("RK", "Runge_Kutta_").replace("_", " ")

def safeexp(x):
     return np.exp(np.clip(x, a_max=700, a_min=-750))

def readfile(name):
    with open(name) as f: tmp = f.read()
    return tmp

def thermrad(a, emiss, t, t_amb):
    sigma = 5.670374419e-8
    return a * emiss * sigma * ( ( t + 273.15 ) ** 4 - ( t_amb + 273.15 ) ** 4 )

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

def double_pendulum_ode(t, x, s, m1, m2, l1, l2, g):
    theta1, theta2, w1, w2 = x[0], x[1], x[2], x[3]
    
    delta = theta1 - theta2
    
    dtheta1 = w1
    dtheta2 = w2
    
    Div = ( 2 * m1 + m2 - m2 * np.cos( 2 * delta ) )
    
    dw1 = ( - g * ( 2 * m2 + m1 ) * np.sin( theta1 ) - m2 * g * np.sin( theta1 - 2 * theta2 ) - 2 * np.sin( delta ) * m2 * ( w2 ** 2 * l2 + w1 ** 2 * l1 * np.cos( delta ) ) ) / ( l1 * Div )
    dw2 = ( 2 * np.sin( delta ) * ( w1 ** 2 * l1 * ( m1 + m2 ) + g * ( m1 + m2 ) * np.cos( theta1 ) + w2 ** 2 * l2 * m2 * np.cos( delta ) ) ) / ( l2 * Div )
    
    dx = np.zeros_like(x)
    dx[0], dx[1], dx[2], dx[3] = dtheta1, dtheta2, dw1, dw2
    return dx 

def crash():
    try:
        crash()
    except:
        crash()


#Notes:
"""
    Recommended Project configuration:
    Project (Folder)
    │
    ├──sim.py
    ├──sim_API.py
    └──config.ini
"""
