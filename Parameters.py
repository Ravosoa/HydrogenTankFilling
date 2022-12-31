# Module importations
import math
import numpy as np
from CoolProp.CoolProp import PropsSI

"""This file contains a dictionary of all parameters needed for the simulation """

# Type 4 tank in general:
C_v = PropsSI('Cvmass', 'P', [k*1.01325e5 for k in [1, 200, 300, 500, 800]], 'T', 25+273.15, 'H2').mean()  # [J/kg/K]
C_p = PropsSI('Cpmass', 'P', [k*1.01325e5 for k in [1, 200, 300, 500, 800]], 'T', 25+273.15, 'H2').mean()  # [J/kg/K]

Rho_70 = PropsSI('Dmass', 'P', [70*10**6], 'T', 15+273.15, 'H2')  # [J/kg/K]

# For CFRP
c_CFRP = 1120  # Specific heat capacity of CFRP in [J/kg/K]
Lambda_CFRP = 0.55  # CFRP thermal conductivity in [W/m/K]
h_conv_out = 6  # Outer forced convection coefficient in [W/m^2/k]
Rho_CFRP = 1360  # Density [kg/m^3]
e_CFRP = 0.005  # Thickness of CFRP in [m]

# For the polymer (here HDPE)
c_poly = 1880  # Specific heat capacity of HDPE in [J/kg/K]
Lambda_poly = 0.55 # HDPE thermal conductivity in [W/m/K] 0.55
h_conv_in = 230 # Inner forced convection coefficient in [W/m^2/k]
Rho_poly = 940  # Density in [kg/m^3]
e_poly = 0.020  # Thickness of HDPE in [m]

# From : Performance of Hydrogen Storage Tanks of Type IV in a Fire: Effect of the State of Charge
# h_out from : "Thermodynamic analysis of the emptying process of compressed hydrogen tanks" @ t = 200s

# For the reservoir : Tank 0
R0_in = 0.282  # Inner radius [m]
R0_out = R0_in + e_poly + e_CFRP  # External radius in [m]
L0 = 2  # Height in [m]



# For the tanks to fill : Tank 1, Tank 2, Tank 3
# From: "The role of initial tank temperature on refuelling of on-board hydrogen tanks"
R_in = 0.115  # Inner radius [m]
R_out = R_in + e_poly + e_CFRP  # External radius in [m]
L = 0.827  # Height in [m]



parameters = {'T_amb': 273 + 30,  # K
              'P_atm': 1.013E5,  # [Pa]
        # H2 properties:
              'C_v': PropsSI('Cvmass', 'P', [k*1.01325e5 for k in [1, 200, 300, 500, 800]], 'T', 25+273.15, 'H2').mean(),  # [J/kg/K]
              'C_p' : PropsSI('Cpmass', 'P', [k*1.01325e5 for k in [1, 200, 300, 500, 800]], 'T', 25+273.15, 'H2').mean(),  # [J/kg/K]
        # For CFRP:
              'c_CFRP': 1120,  # Specific heat capacity of CFRP in [J/kg/K]

        # For the polymer (here HDPE):
              'c_poly': 1880,  # Specific heat capacity of HDPE in [J/kg/K]

        # For the reservoir : Tank 0
              'V0': math.pi * pow(R0_in, 2) * L0,  # [m^3]
              'm0_poly': Rho_poly * L0 * math.pi * (pow(R0_in + e_poly, 2)-pow(R0_in, 2)),   # [kg]
              'm0_CFRP': Rho_CFRP * L0 * math.pi * (pow(R0_out, 2)-pow(R0_in + e_poly, 2)),   # [kg]
              'UA0_in': 1 / ((1/(2*math.pi*R0_in*L0*h_conv_in))+(np.log((R0_in + e_poly)/R0_in)/(2*math.pi*L0*Lambda_poly))),  # [W/K]
              'UA0_out': 1 / ((1/(2*math.pi*R0_out*L0*h_conv_out))+(np.log(R0_out/(R0_in+e_poly))/(2*math.pi*L0*Lambda_CFRP))),  # [W/K]
              'UA0_int': 1 / ((np.log((R0_in + e_poly)/R0_in)/(2*math.pi*L0*Lambda_poly)) +
                              (np.log(R0_out/(R0_in+e_poly))/(2*math.pi*L0*Lambda_CFRP))),  # [W/K]

        # For the tanks to fill : Tank 1, Tank 2, Tank 3
              'V1': math.pi * pow(R_in, 2) * L,  # [m^3]
              'V2': math.pi * pow(R_in, 2) * L,  # [m^3]
              'V3': math.pi * pow(R_in, 2) * L,  # [m^3]
              'm_poly': Rho_poly * L * math.pi * (pow(R_in + e_poly, 2)-pow(R_in, 2)),    # [kg]
              'm_CFRP': Rho_CFRP * L * math.pi * (pow(R_out, 2)-pow(R_in + e_poly, 2)),   # [kg]
              'UA_in': 1 / ((1/(2*math.pi*R_in*L*h_conv_in))+(np.log((R_in + e_poly)/R_in)/(2*math.pi*L*Lambda_poly))),  # [W/K]
              'UA_out': 1 / ((1/(2*math.pi*R_out*L*h_conv_out))+(np.log(R_out/(R_in+e_poly))/(2*math.pi*L*Lambda_CFRP))),  # [W/K]
              'UA_int': 1 / ((np.log((R_in + e_poly)/R_in)/(2*math.pi*L*Lambda_poly)) +
                             (np.log(R_out/(R_in+e_poly))/(2*math.pi*L*Lambda_CFRP))),  # [W/K]

        # Cross-section of pipes:
              'A_0': math.pi * 0.005**2,  # [m^2]
              'A_1': math.pi * 0.005**2,  # [m^2]
              'A_2': math.pi * 0.005**2,
              'A_3': math.pi * 0.005**2,

        # For the equations of states:
              'r': 4124.2,  # [J/kgK]
              'a': 3.6E4,  # [J K^(1/2)/ kg^2 / m^3]
              'b': 9.2E-3,  # [m^3 / kg]

        # Heat exchanger:
              'cp_H2O': PropsSI('Cpmass', 'P', [k*1.01325e5 for k in [1, 200, 300, 500, 800]], 'T', 25+273.15, 'Water').mean(),  # [J/kg/K]
        # Limitation at T = 15°C and P = 70MPa
              'Rho_70': PropsSI('Dmass', 'P', [70*10**6], 'T', 15+273.15, 'H2')  # [J/kg/K]
              }


