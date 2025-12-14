import numpy as np
from sim_API import *
m1, m2, l1, l2, g = 1, 1, 1, 1, 9.81 #current example: double pendulum
def df(t, x, s):
    Theta1, Theta2, w1, w2 = x[0], x[1], x[2], x[3]
    x[0], x[1], x[2], x[3] = w1, w2, ( - g * ( 2 * m2 + m1 ) * np.sin( Theta1 ) - m2 * g * np.sin( Theta1 - 2 * Theta2 ) - 2 * np.sin( Theta1 - Theta2 ) * m2 * ( w2 ** 2 * l2 + w1 ** 2 * l1 * np.cos( Theta1 - Theta2 ) ) ) / ( l1 * ( 2 * m1 + m2 - m2 * np.cos( 2 * Theta1 - Theta2 ) ) ), ( 2 * np.sin( Theta1 - Theta2 ) * ( w1 ** 2 * l1 * ( m1 + m2 ) + g * ( m1 + m2 ) * np.cos( Theta1 ) + w2 ** 2 * l2 * m2 * np.cos( Theta1 - Theta2 ) ) ) / ( l2 * ( 2 * m1 + m2 - m2 * np.cos( 2 * Theta1 - Theta2 ) ) )
    return x 
def Rec_f(State, dState): return dState[0], dState[1]
TIME, Rec = run_sim(2e-3, np.array( [0], dtype = np.float64 ), np.array( [0.2, 0, 0, 0], dtype = np.float64 ), 5, 0, no_f, df, "GLRK4", "fS", Rec_f, True)
plot("Graph", + l1 * np.sin(Rec) + l2 * np.sin(TIME), - l1 * np.cos(Rec) - l2 * np.cos(TIME))
