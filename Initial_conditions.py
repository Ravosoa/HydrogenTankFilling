from Parameters import parameters as p
from scipy.optimize import fsolve
import numpy as np

def density(P, T):
    """Compute the density with Redlich-Kwong Equation of State.

            Args:
                P: Pressure [Pa]
                T: Temperature [K]

            Returns:
                Density [kg/m^3]
    """
    def func(dens):
        return p['r'] * T / (1 / dens - p['b']) - p['a'] / (np.sqrt(T) * (1 / dens) * ((1 / dens) + p['b'])) - P

    z = fsolve(func, 100)
    np.isclose(func(z), 0)  # func(root) should be almost 0.0
    return z[0]


def internal_energy(T, rho):
    """Compute the internal energy from a derivation of Redlich-Kwong Equation of State.

            Args:
                rho: Density [kg/m^3]
                T: Temperature [K]

            Returns:
                Internal energy [J/kg]
    """
    return p['C_v'] * T + 3 * p['a'] * np.log(1 / (1 + p['b'] * rho)) / (2 * p['b'] * np.sqrt(T))

# Define the inputs
u_0 = 0.05  # Mass flow rate m_dot_1 in [kg/s]
u_1 = 0.05  # Mass flow rate m_dot_2 in [kg/s]
u_2 = 0.05  # Mass flow rate m_dot_3 in [kg/s]
u_3 = 70000  # Q_hex_1_dot [J/s]
u_4 = 70000  # Q_hex_2_dot [J/s]
u_5 = 70000  # Q_hex_3_dot [J/s]

u0 = np.array([u_0, u_1, u_2, u_3, u_4, u_5])

timestep_0 = 0.1

# Initial conditions :

# For the reservoir
P_0i = 300*10**5 # [Pa]
T_0i = 300 + 25  # [K]
rho_0i = density(P_0i, T_0i)
u_0i = internal_energy(T_0i, rho_0i)

# We assume that initially the tank 1, 2, 3 are empty but still remaining H2 inside.
# Thus, the density of H2 is computed at atmospheric pressure and ambient temperature
T_1i = 273 + 25  # [K]
P_1i = 2*10**6  # [Pa]
rho_1i = density(P_1i, T_1i)
u_1i = internal_energy(T_1i, rho_1i)
h_0i = u_0i + P_0i / rho_0i
h_1i = u_1i + p['P_atm'] / rho_1i


masse_0 = rho_1i*p['V1']
x0 = np.array([u_0i,  # x0 - u0
               u_1i,  # x1 - u1
               u_1i,   # x2 - u2
               u_1i,   # x3 - u3
               rho_0i*(p['V0']),  # x4 - m0
               rho_1i*p['V1'],    # x5 - m1
               rho_1i*p['V2'],    # x6 - m2
               rho_1i*p['V3'],    # x7 - m3
               T_0i,   # x8 - T_CFRP0
               T_1i,   # x9 - T_CFRP1
               T_1i,   # x10 - T_CFRP2
               T_1i,   # x11 - T_CFRP3
               T_0i,     # x12 - T_poly0
               T_1i,   # x13 - T_poly1
               T_1i,   # x14 - T_poly2
               T_1i,   # x15 - T_poly3
               rho_0i,   # x16 - rho0
               rho_1i,   # x17 - rho1
               rho_1i,   # x18 - rho2
               rho_1i,  # x19 - rho3
               ])


# initial guess for the research of zeros for z

z0 = [(u0[0] + u0[1] + u0[2]) / (x0[16] * p['A_0']),  # z0
      u0[0] / (x0[16] * p['A_1']),  # z1
      u0[1] / (x0[16] * p['A_2']),  # z2
      u0[2] / (x0[16] * p['A_3']),  # z3
      P_0i,  # z4 - P_valve_1
      P_0i,  # z5 - P_valve_2
      P_0i,  # z6 - P_valve_3
      T_0i,  # z7 - T_valve_1
      T_0i,  # z8 - T_valve_2
      T_0i,  # z9 - T_valve_3
      h_0i,  # z10 - h_valve_1
      h_0i,  # z11 - h_valve_2
      h_0i,  # z12 - h_valve_3
      T_0i,  # z13 - T_hex_1
      T_0i,  # z14 - T_hex_2
      T_0i,  # z15 - T_hex_3
      h_0i,  # z16 - h_0
      T_1i,  # z17 - T_0
      T_1i,  # z18 - T_1
      T_1i,  # z19 - T_2
      T_1i,  # z20 - T_3
      P_0i,  # z21 - P_0
      P_1i,  # z22 - P_1
      P_1i,  # z23 - P_2
      P_1i,  # z24 - P_3
      ]

# Creating bounds
T_min = 273 - 40
T_max = 273 + 300
P_min = 0
P_max = 70 * 10**7
P_max_0 = 300 * 10**5
v_min = 0
v_max = 500
h_min = 0
h_max = 10E13
Q_min = 0
Q_max = 100000
rho_min = 0
rho_max = 100
mass_max = 5
mass_min = 0
mass_max_0 = 50
MFR_min = 0  # Mass flow rate
MFR_max = 0.6
u_min = 0  # Internal energy
u_max = 5E6


Bounds_z = ((v_min, v_max), (v_min, v_max), (v_min, v_max), (v_min, v_max), (P_min, P_max_0), (P_min, P_max_0), (P_min, P_max_0),
            (T_min, T_max), (T_min, T_max), (T_min, T_max), (h_min, h_max), (h_min, h_max), (h_min, h_max),
            (T_min, T_max), (T_min, T_max), (T_min, T_max), (h_min, h_max), (T_min, T_max), (T_min, T_max), (T_min, T_max),
            (T_min, T_max), (P_min, P_max_0), (P_min, P_max), (P_min, P_max), (P_min, P_max))

Bounds_x = ((u_min, u_max), (u_min, u_max), (u_min, u_max), (u_min, u_max), (mass_min, mass_max_0), (mass_min, mass_max),
            (mass_min, mass_max), (mass_min, mass_max), (T_min, T_max), (T_min, T_max), (T_min, T_max), (T_min, T_max),
            (T_min, T_max), (T_min, T_max), (T_min, T_max), (T_min, T_max), (rho_min, rho_max), (rho_min, rho_max),
            (rho_min, rho_max), (rho_min, rho_max))

Bounds_u = ((MFR_min, MFR_max), (MFR_min, MFR_max), (MFR_min, MFR_max),
            (Q_min, Q_max), (Q_min, Q_max), (Q_min, Q_max))

Bound_h = (0.0001, 2)

