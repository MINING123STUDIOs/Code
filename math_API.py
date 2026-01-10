"""
Contents:
    function
    Notes
"""

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
        c = c + f(n)
    return c

def pi_prod(f, a, b): #\prod_{n=a}^{b}f(n)
    c = 1
    for n in range(a,b+1):
        c = c * f(n)
    return c

def scalar_solve(F, x0, tol=1e-9, max_iter=20):
    x = x0
    eps = 1e-8
    for i in range(max_iter):
        x = x - eps * F(x) / ( F( x + eps ) - F( x ) )
        if abs(F(x)) < tol:
            return x
    return x

def deriv(f, w, eps=1e-8): #\frac{df}{dx}|x=w
    return ( f( w + eps / 2 ) - f( w - eps / 2 ) ) / eps

def s_deriv(f, w, eps=1e-8): #\frac{d^2f}{dx^2}|x=w
    return ( f( w + eps) - 2 * f( w ) + f( w - eps) ) / eps ** 2

def sign(x):
    if   x > 0: w = +1
    elif x < 0: w = -1
    else      : w = 0
    return w

def hstep(x, z=0.5):
    if   x > 0: w = 1
    elif x < 0: w = 0
    else      : w = z
    return w

def clamp(x, hi, lo):
    if   x > hi: w = hi
    elif x < lo: w = lo
    else       : w = x
    return w

def linint(hi, lo, s):
    if s < 0.0:
        s = 0.0
    elif s > 1.0:
        s = 1.0
    return hi * s + lo * ( 1 - s )
"""
def is_prime(x):
    if x % 1 > 0: return False
    
    else: return 
"""

#Notes:
"""
    none.
"""