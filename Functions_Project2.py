import numpy as np
import matplotlib.pyplot as plt
from Parameters import *


def LotkaVolterraModel (x, params):
    alpha = params["alpha"]
    beta = params["beta"]
    gamma = params["gamma"]
    delta = params["delta"]

    xdot = np.array([alpha*x[0] - beta*x[0]*x[1],
                     delta*x[0]*x[1] - gamma*x[1]])

    return xdot


# RK-4 method:
def RK4(f, x0, t0, tf, h):
    # Function to solve
    # x0 : Initial conditions
    # t0 : initial time
    # tf : final time
    # h : time step

    t = np.arange(t0, tf, h)
    nt = t.size  # size of time

    nx = x0.size
    x = np.zeros((nx, nt))

    x[:, 0] = x0

    for i in range(nt-1):
        k1 = h * f(t[i], x[:, i])
        k2 = h * f(t[i] + h/2, x[:, i] + k1/2)
        k3 = h * f(t[i] + h/2, x[:, i] + k2/2)
        k4 = h * f(t[i] + h, x[:, i] + k3)

        dx = (k1 + 2*k2 + 2*k3 + k4) / 6

        x[:, i+1] = x[:, i] + dx

    return x, t


# Define Problem
params = {"alpha": 1.1, "beta": 0.4, "gamma": 0.4, "delta": 0.1}

f = lambda t, x: LotkaVolterraModel(x, params)

x0 = np.array([20,5])

# solve the Diff. Eq

t0 = 0
tf = 100
h = 0.01

x, t = RK4(f,x0, t0, tf, h)

# Plot the results
#plt.subplot(1, 2)
plt.plot(t, x[0, :], 'r', label="Preys")
plt.plot(t, x[1, :], 'b', label="Predators")
plt.xlabel("Time (t)")
plt.grid()
plt.legend()
plt.show()