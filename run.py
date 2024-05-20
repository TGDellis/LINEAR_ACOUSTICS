import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import os

n = 201
m = 201
dt = 1e-3
t_total = 0.8
t_total_1 = 0.15
t_total_2 = 0.4
DL = 130
N = n + 2 * DL
M = m + 2 * DL
x = np.linspace(-100 - DL, 100 + DL, N)
y = np.linspace(-100 - DL, 100 + DL, M)

r_m = np.empty((len(x) * len(y), 2))
i = 0
for x_1 in x:
    for y_1 in y:
        r_m[i, :] = np.array([x_1, y_1])
        i = i + 1
np.savetxt("grid.txt", r_m, fmt="%lf")
np.savetxt("x.txt", x, fmt="%lf")
np.savetxt("y.txt", y, fmt="%lf")
ins = np.empty((1, 6))
ins[0] = np.array([N, M, dt, t_total, t_total_1, t_total_2])
np.savetxt("input.txt", ins, fmt="%d %d %lf %lf %lf %lf")
Mach = -0.3
f = 10 + 1  # Hertz
Th_0 = len("dellis") + 10
a = 1 / len("dellis")
T_0 = 273 + Th_0
pars = np.empty((1, 4))
pars[0] = np.array([Mach, T_0, f, a])
np.savetxt("params.txt", pars, fmt="%lf %lf %lf %lf")
os.system("./run > log")
