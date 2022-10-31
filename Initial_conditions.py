import math
import numpy as np
from CoolProp.CoolProp import PropsSI
from Parameters import *

def printx(x):
    print(f"u0 = {x[0]}, u1 = {x[1]}, u2 = {x[2]}, u3 = {x[3]}")
    print(f"m0 = {x[4]}, m1 = {x[5]}, m2 = {x[6]}, m3 = {x[7]}")
    print(f"T_CFRP0 = {x[8]}, T_CFRP1 = {x[9]}, T_CFRP2 = {x[10]}, T_CFRP3 = {x[11]}")
    print(f"T_poly0 = {x[12]}, T_poly1 = {x[13]}, T_poly2 = {x[14]}, T_poly3 = {x[15]}")
    print(f"rho_0 = {x[16]}, rho_1 = {x[17]}, rho_2 = {x[18]}, rho_3 = {x[19]}")
    print(f"Q_hex1 = {x[20]}, Q_hex2 = {x[21]}, Q_hex3 = {x[22]}")
    print(f"m_H20_1 = {x[23]}, m_H20_2 = {x[24]}, m_H20_3 = {x[25]}")

# Initial conditions :
#e = 0.000001 # Density close to 0
P_0i = 77  # [MPa]
T_0i = 30 + 273  # [K]
rho_0i = PropsSI('D', 'P', P_0i*(10**5), 'T', T_0i, 'H2')
u_0i = C_v * T_0i + 3 * a * np.log(1/(1 + b * rho_0i)) / (2 * b * np.sqrt(T_0i))

# We assume that initialy the tank 1, 2, 3 are empty but still remaining H2 inside.
# Thus, the density of H2 is computed at atmospheric pressure and ambiant temperature
rho_1i = PropsSI('D', 'P', P_atm, 'T', T_amb, 'H2')
#print(rho_0i)
u_1i = C_v * T_amb + 3 * a * np.log(1/(1 + b * rho_1i)) / (2 * b * np.sqrt(T_amb))

x0 = np.array([u_0i,  # x0 - u0
               u_1i,  # x1 - u1
               u_1i,   # x2 - u2
               u_1i,   # x3 - u3
               rho_0i*V0,   # x4 - m0
               rho_1i*V1,   # x5 - m1
               rho_1i*V2,   # x6 - m2
               rho_1i*V3,   # x7 - m3
               T_0i,   # x8 - T_CFRP0
               T_amb,   # x9 - T_CFRP1
               T_amb,   # x10 - T_CFRP2
               T_amb,   # x11 - T_CFRP3
               T_0i,     # x12 - T_poly0
               T_amb,   # x13 - T_poly1
               T_amb,   # x14 - T_poly2
               T_amb,   # x15 - T_poly3
               rho_0i,   # x16 - rho0
               rho_1i,   # x17 - rho1
               rho_1i,   # x18 - rho2
               rho_1i,   # x19 - rho3
               0,   # x20 - Q_hex1
               0,   # x21 - Q_hex2
               0,   # x22 - Q_hex3
               0,   # x23 - m_H20_1
               0,   # x24 - m_H20_2
               0])   # x25 - m_H20_3

#printx(x0)