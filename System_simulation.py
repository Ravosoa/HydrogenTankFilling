import matplotlib.pyplot as plt
from Initial_conditions import *
import sys


def temperature_root(p):
    t_array = np.roots(p)
    #print(t_array)
    t = T_amb
    for i in range(len(t_array)):
        if isinstance(t_array[i], complex):
            continue
        elif t_array[i] <= 0:
            continue
        else:
            t = t_array[0]
    return t


def LotkaVolterraModel(x, u):
    # State space representation
    v = np.array([(temperature_root([C_v, 0, -x[0], (3*a / (2*b)) * np.log(1 / (1 + b*x[16]))]))**2,
                  (temperature_root([C_v, 0, -x[1], (3*a / (2*b)) * np.log(1 / (1 + b*x[17]))]))**2,
                  (temperature_root([C_v, 0, -x[2], (3*a / (2*b)) * np.log(1 / (1 + b*x[18]))]))**2,
                  (temperature_root([C_v, 0, -x[3], (3*a / (2*b)) * np.log(1 / (1 + b*x[19]))]))**2])

    y = np.array([r * v[0] / (1 / x[16] - b) - a / (np.sqrt(v[0])*(1/x[16])*(b + 1/x[16])) if (x[16] != 0) and (v[0] != 0) else P_atm,
                  r * v[1] / (1 / x[17] - b) - a / (np.sqrt(v[1])*(1/x[17])*(b + 1/x[17])) if (x[17] != 0) and (v[1] != 0) else P_atm,
                  r * v[2] / (1 / x[18] - b) - a / (np.sqrt(v[2])*(1/x[18])*(b + 1/x[18])) if (x[18] != 0) and (v[2] != 0)else P_atm,
                  r * v[3] / (1 / x[19] - b) - a / (np.sqrt(v[3])*(1/x[19])*(b + 1/x[19])) if (x[19] != 0) and (v[3] != 0) else P_atm])

    d = np.array([(u[0] + u[1] + u[2]) / x[16] if x[16] != 0 else 0,
                  C_D10 * A_1 * (y[1] - y[0]) / u[0] if u[0] != 0 else 0,
                  C_D20 * A_2 * (y[2] - y[0]) / u[1] if u[1] != 0 else 0,
                  C_D30 * A_3 * (y[3] - y[0]) / u[2] if u[2] != 0 else 0])

    w = np.array([u[0] / (A_1 * d[1]) if d[1] != 0 else 0,
                  u[1] / (A_2 * d[2]) if d[2] != 0 else 0,
                  u[2] / (A_3 * d[3]) if d[3] != 0 else 0,
                  y[1] - u[0] * d[1] / (C_D11 * A_1),
                  y[2] - u[1] * d[2] / (C_D11 * A_2),
                  y[3] - u[2] * d[3] / (C_D11 * A_3),
                  ])

    e = np.array([(temperature_root([r / (1 / w[0] - b), -w[3], -a / ((1 / w[0]) * (b + 1 / w[0]))])) ** 2 if w[0] != 0 else T_amb,
                  (temperature_root([r / (1 / w[1] - b), -w[4], -a / ((1 / w[1]) * (b + 1 / w[1]))])) ** 2 if w[1] != 0 else T_amb,
                  (temperature_root([r / (1 / w[2] - b), -w[5], -a / ((1 / w[2]) * (b + 1 / w[2]))])) ** 2 if w[2] != 0 else T_amb])
    if(v[0] == 0):
        sys.breakpointhook()
    s = np.array([C_v * e[0] + 3 * a / (2 * b * np.sqrt(e[0])) * np.log(1 / (1 + b * w[0])) + r * e[0] / (
                              1 - b * w[0]) - a / (np.sqrt(e[0]) * (b + 1 / w[0])) if (e[0] != 0) and (w[0] != 0) else 0,
                  C_v * e[1] + 3 * a / (2 * b * np.sqrt(e[1])) * np.log(1 / (1 + b * w[1])) + r * e[1] / (
                              1 - b * w[1]) - a / (np.sqrt(e[1]) * (b + 1 / w[1])) if (e[1] != 0) and (w[1] != 0) else 0,
                  C_v * e[2] + 3 * a / (2 * b * np.sqrt(e[2])) * np.log(1 / (1 + b * w[2])) + r * e[2] / (
                          1 - b * w[2]) - a / (np.sqrt(e[2]) * (b + 1 / w[2])) if (e[2] != 0) and (w[2] != 0) else 0,
                  x[20] / (x[23] * cp_H2O) + u[3] if x[23] != 0 else T_amb,
                  x[21] / (x[24] * cp_H2O) + u[4] if x[24] != 0 else T_amb,
                  x[22] / (x[25] * cp_H2O) + u[5] if x[25] != 0 else T_amb,
                  x[20] / (C_p * x[5]) + e[0] if x[5] != 0 else 0,
                  x[21] / (C_p * x[6]) + e[1] if x[6] != 0 else 0,
                  x[22] / (C_p * x[7]) + e[2] if x[7] != 0 else 0,
                  x[0] + r * v[0] / (1 - b * x[16]) - a / (np.sqrt(v[0]) * (b + 1 / x[16]))])

    z = np.concatenate((d, w, e, s, v, y), axis=0)


    xdot = np.array([-1/x[4] * ((u[0]+u[1]+u[2])*(z[22] + (((u[0]+u[1]+u[2])/(A_0*x[16]))**2)/2) - UA0_in*(x[12]-z[23])
                               - x[0]*(u[0]+u[1]+u[2])) if (x[4] != 0) else z[22] - z[27]/x[16],                   # x0
                     1/x[5] * (u[0]*(z[13] + x[20]/x[5] + ((u[0]/(A_1*x[17]))**2)/2) - UA_in*(x[13]-z[24])
                               - x[1]*u[0]) if x[5] != 0 else 0,                                                   # x1
                     1/x[6] * (u[1]*(z[14] + x[21]/x[6] + ((u[1]/(A_1*x[18]))**2)/2) - UA_in*(x[14]-z[25])
                               - x[2]*u[1]) if x[6] != 0 else 0,                                                   # x2
                     1/x[7] * (u[2]*(z[15] + x[22]/x[7] + ((u[2]/(A_1*x[19]))**2)/2) - UA_in*(x[15]-z[26])
                               - x[3]*u[2]) if x[7] != 0 else 0,                                                   # x3
                     u[0] + u[1] + u[2],                                                                           # x4
                     u[0],                                                                                         # x5
                     u[1],                                                                                        # x6
                     u[2],                                                                                         # x7
                     1 / (m0_CFRP * c_CFRP) * (UA0_int * (x[12] - x[8]) - UA0_out * (x[8] - T_amb)),               # x8
                     1 / (m_CFRP * c_CFRP) * (UA_int * (x[13] - x[9]) - UA0_out * (x[9] - T_amb)),               # x9
                     1 / (m_CFRP * c_CFRP) * (UA_int * (x[14] - x[10]) - UA0_out * (x[10] - T_amb)),               # x10
                     1 / (m_CFRP * c_CFRP) * (UA_int * (x[15] - x[11]) - UA0_out * (x[11] - T_amb)),               # x11
                     1 / (m0_poly * c_poly) * (UA0_in * (z[23] - x[12]) - UA0_int * (x[12] - x[8])),               # x12
                     1 / (m_poly * c_poly) * (UA_in * (z[24] - x[13]) - UA_int * (x[13] - x[9])),               # x13
                     1 / (m_poly * c_poly) * (UA_in * (z[25] - x[14]) - UA_int * (x[14] - x[10])),               # x14
                     1 / (m_poly * c_poly) * (UA_in * (z[26] - x[15]) - UA_int * (x[15] - x[11])),               # x15
                     -(u[0] + u[1] + u[2]) / V0,                                                              # x16
                     u[0] / V1,                                                                            # x17
                     u[1] / V2,                                                                              # x18
                     u[2] / V3,                                                                                # x19
                     UA_hex * ((z[19] - u[3]) - (z[10] - z[16])) / (np.log((z[19] - u[3]) / (z[10] - z[16])))
                     if ((z[10] - z[16])*(z[19] - u[3]) > 0) and
                        ((np.log((z[19] - u[3]) / (z[10] - z[16]))) != 0) else 0,                              # x20
                     UA_hex * ((z[20] - u[4]) - (z[11] - z[17])) / (np.log((z[20] - u[4]) / (z[11] - z[17])))
                     if ((z[11] - z[17])*(z[20] - u[4]) > 0) and
                        ((np.log((z[20] - u[4]) / (z[11] - z[17]))) != 0) else 0,                          # x21
                     UA_hex * ((z[21] - u[5]) - (z[12] - z[18])) / (np.log((z[21] - u[5]) / (z[12] - z[18])))
                     if ((z[12] - z[18])*(z[21] - u[5]) > 0) and
                        (np.log((z[21] - u[5]) / (z[12] - z[18])) != 0) else 0,                          # x22
                     u[6],                                                                             # x23
                     u[7],                                                                                 # x24
                     u[8]])                                                                                # x25

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


# Define the inputs
u_0 = 0.005  # Mass flow rate m_dot_1 in [kg/s]
u_1 = 0  # Mass flow rate m_dot_2 in [kg/s]
u_2 = 0  # Mass flow rate m_dot_3 in [kg/s]
u_3 = 25 + 273  # Temperature of incoming water in the heat exchanger 1 [K]
u_4 = 25 + 273  # Temperature of incoming water in the heat exchanger 2 [K]
u_5 = 25 + 273  # Temperature of incoming water in the heat exchanger 3 [K]
u_6 = 0  # Mass flow rate of water entering the heat exchanger 1 [kg/s]
u_7 = 0  # Mass flow rate of water entering the heat exchanger 2 [kg/s]
u_8 = 0  # Mass flow rate of water entering the heat exchanger 2 [kg/s]

u = np.array([u_0, u_1, u_2, u_3, u_4, u_5, u_6, u_7, u_8])

f = lambda t, x: LotkaVolterraModel(x, u)

# solve the Diff. Eq

t0 = 0  # [s]
tf = 300  # [s]
h = 0.01

x, t = RK4(f, x0, t0, tf, h)

# Plot the results

#plt.plot(t, x[0, :], 'r', label="internal energy in Tank 0")
#plt.plot(t, x[13, :], 'b', label="Temperature of the liner of Tank 1 ")
#plt.xlabel("Time (t)")
#plt.grid()
#plt.legend()
#plt.show()


fig, axs = plt.subplots(4, 2)
axs[0, 0].plot(t, x[0])
axs[0, 0].set_title('u_0')
axs[0, 1].plot(t, x[1], 'tab:orange')
axs[0, 1].set_title('u_1')
axs[1, 0].plot(t, x[12], 'tab:green')
axs[1, 0].set_title('T_poly_0')
axs[1, 1].plot(t, x[13], 'tab:red')
axs[1, 1].set_title('T_poly_1')

axs[2, 0].plot(t, x[4])
axs[2, 0].set_title('m_0')
axs[2, 1].plot(t, x[5], 'tab:orange')
axs[2, 1].set_title('m_1')
axs[3, 0].plot(t, x[16], 'tab:green')
axs[3, 0].set_title('Rho_0')
axs[3, 1].plot(t, x[17], 'tab:red')
axs[3, 1].set_title('Rho_1')


plt.show()

