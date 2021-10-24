"""
Nonlinear method in civil engineering_ hw1
by CHEN Jiawei, master 2nd year in Mechanical engineering
Here is the .py for forward euler algorithm
"""

import math
import matplotlib.pyplot as plt

"""
    define function
"""
def formula_Kt(k_bar, u):
    Kt = 1 + 3*u + 1.5*pow(u,2) + k_bar
    return Kt


"""
     init param
"""
k_bar = 0.8
df = -0.05
Fmax = -2
num_step = Fmax/df
Kt = []
ui = []
Fi = []

"""
    First step
"""
Kt.append(0);
ui.append(0);
Fi.append(0);

"""
    start forward euler
"""
for i in range(1,int(num_step)+1):
    Kt.append(formula_Kt(k_bar,ui[i-1]))
    du = df / Kt[i]
    ui.append(ui[i-1]+du)
    Fi.append(Fi[i-1]+df)


plt.plot(ui,Fi)
plt.show()