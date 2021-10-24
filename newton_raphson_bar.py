"""
Nonlinear method in civil engineering_ hw1
by CHEN Jiawei, master 2nd year in Mechanical engineering
Here is the .py for newton raphson algorithm
"""
import matplotlib.pyplot as plt

"""
    define function
"""
def formula_Kt(k_bar, u):
    Kt = 1 + 3*u + 1.5*pow(u,2) + k_bar
    return Kt

def formula_F(k_bar, u):
    F = k_bar*u + 0.5*pow(u,3) + 1.5*pow(u,2) + u
    return F

def formula_err(current,target):
    err = abs((target-current)/target)
    return err

"""
     init param
"""
k_bar = 0.8
df = -0.02
Fmax = -2
num_step = Fmax/df
tol = 1e-6
Kt = []
ui = []
Fi = []

"""
    First step
"""
Kt.append(0)
ui.append(0)
Fi.append(0)

"""
    start 
"""
for i in range(1,int(num_step)+1):
    ui_j = []
    Fi_j = []
    ui_j.append(ui[-1])
    Fi_j.append(Fi[-1])
    Fmax_j = Fi[-1] + df
    iter = 0
    err = 1
    while err > tol:
        Kt_1 = formula_Kt(k_bar, ui_j[-1])
        residual = Fmax_j - Fi_j[-1]
        du_j = residual/Kt_1
        ui_j.append(ui_j[-1]+du_j)

        iter = iter + 1

        Fi_j.append(formula_F(k_bar, ui_j[-1]))

        err = formula_err(Fi_j[-1],Fmax_j)
        #print(err)
    ui.append(ui[-1]+ui_j[-1])
    Fi.append(Fi[-1]+Fi_j[-1])
    print(Fi[-1])

plt.plot(ui,Fi)
plt.show()
