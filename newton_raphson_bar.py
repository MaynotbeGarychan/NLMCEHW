"""
Nonlinear method in civil engineering_ hw1
by CHEN Jiawei, master 2nd year in Mechanical engineering
Here is the .py for newton raphson algorithm
"""
import matplotlib.pyplot as plt
import numpy as np
import math

"""
    define function
"""
def formula_Kt(k_bar, u):
    Kt = (1 + 3 * u + 0.5 * pow(u, 2)) + k_bar
    return Kt

def formula_F(k_bar, u):
    F = (u + 1.5 * pow(u, 2) + 0.5 * pow(u, 3)) + k_bar * u
    return F

"""
     init param
"""
k_bar = 0.8
df = -0.02
Fmax = -2
num_steps = math.ceil(Fmax/df)
tol = 1e-6

"""
    Make new array to store the step info
    Notation for the step is "i"
    for initial state, u_i, F_i = 0
"""
u_i = np.zeros(num_steps)
F_i = np.zeros(num_steps)

"""
    First step
"""
for i in range(1,num_steps):
    u_i_j = []
    F_i_j = []
    u_i_j.append(u_i[i-1])
    F_i_j.append(F_i[i-1])
    F_i[i] = F_i[i-1]+df
    iter = 1
    r = 1
    du = 0
    Kt_i_j = formula_Kt(k_bar, u_i[i-1])
    while r > tol:
        #Kt_i_j = formula_Kt(k_bar,u_i_j[iter-1])
        R = F_i[i] - F_i_j[iter-1]
        if R < 0:
            Kt_i_j = 2
        du_j = R/Kt_i_j
        du = du + du_j
        u_i_j.append(u_i_j[iter-1]+du_j)
        F_i_j.append(formula_F(k_bar,u_i_j[iter]))
        r = abs(F_i[i] - F_i_j[iter])
        iter = iter + 1
    u_i[i] = u_i[i-1] + du
    print("step forward")
    print(F_i[i])

fig,axs = plt.subplots(sharex=False,figsize=(10,8))
axs.scatter(u_i,F_i)
plt.show()