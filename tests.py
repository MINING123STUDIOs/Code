#tests
from     sim_API import *
from    math_API import *
from special_API import *

t1 = abs(integrate(lambda x:x**2,0,1,100)-1/3)<1e-4
print(t1)

t2 = abs(sigma_sum(lambda x:x, 1, 5)-15)<1e-4
print(t2)

t3 = abs(pi_prod(lambda x:x, 1, 5)-120)<1e-4
print(t3)

t4 = abs(scalar_solve(lambda x:x**2-2, 1)-2**0.5)<1e-4
print(t4)

t5 = abs(deriv(lambda x:x**2, 2)-4)<1e-4
print(t5)

t6 = abs(s_deriv(lambda x:x**2, 2))<1e-4
print(t6)

t7 = abs(sign(-2)+1)<1e-4
print(t7)

t8 = MD5("") == "d41d8cd98f00b204e9800998ecf8427e"
print(t8)

t9 = SHA256("") == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
print(t9)

t10 = norm_solver_name("RK2") == "Runge Kutta 2"
print(t10)

t11 = clean("text") == "text" and clean("&") == ""
print(t11)

t12 = len(rngstr(5)) == 5
print(t12)

t13 = abs(safeexp(900)-safeexp(700))<1e-4
print(t13)

t14 = abs(thermrad(5, 1, 30, 30))<1e-4
print(t14)

t15 = abs(dielectric_breakdown(0, 1, 25, 3))<1e-4
print(t15)

t16 = abs(diode(-50, 1e-12, 0.025))<1e-4
print(t16)

t17 = abs(linint(1, 0, 0.5)-0.5)<1e-4
print(t17)

t18 = abs(no_f(5, 2, 6)-6)<1e-4
print(t18)

input()