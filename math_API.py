"""
Contents:
    imports
    function
    Notes
"""

#imports:
#import numpy as np
#import scipy as sci
#import scipy.optimize as op

#functions:

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

def sigma_sum(f, a, b): #\sum_{n=a}^{b}f(n)
    c = 0
    for n in range(a,b+1):
        c += f(n)
    return c

def pi_prod(f, a, b): #\prod_{n=a}^{b}f(n)
    c = 1
    for n in range(a,b+1):
        c = c * f(n)
    return c





#Notes:
"""
    currently none
"""