import numpy as np
from sim_API import *
from special_API import *
r=run_sim(2e-3,np.array([0]),np.array([0.2,0,0,0]),5,0,no_f,lambda t,x,s:double_pendulum_ode(t,x,s,1,1,1,1,9.81),"GLRK2","fS")
plot("Graph",np.array([np.sin(r[1,:])+np.sin(r[2,:]),-np.cos(r[1,:])-np.cos(r[2,:])]),"a")