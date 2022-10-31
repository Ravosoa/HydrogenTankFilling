# Module importations
import math
import numpy as np
from ht import *
from CoolProp.CoolProp import PropsSI


T_amb = 273 + 20  # K
P_atm = 1.013E5  # [Pa]

# Type 4 tank in general:

Rho_h2 = PropsSI('D', 'P', [k*1.01325e5 for k in [1, 200, 300, 500, 800]], 'T', 25+273.15, 'H2').mean()  # [kg/m^3]
C_v = PropsSI('Cvmass', 'P', [k*1.01325e5 for k in [1, 200, 300, 500, 800]], 'T', 25+273.15, 'H2').mean()  # [J/kg/K]
C_p = PropsSI('Cpmass', 'P', [k*1.01325e5 for k in [1, 200, 300, 500, 800]], 'T', 25+273.15, 'H2').mean()  # [J/kg/K]

# For CFRP
c_CFRP = 1120  # Specific heat capacity of CFRP in [J/kg/K]
Lambda_CFRP = 0.9  # CFRP thermal conductivity in [W/m/K]
h_conv_out = 6  # Inner forced convection coefficient in [W/m^2/k]
Rho_CFRP = 1360  # Density [kg/m^3]
e_CFRP = 0.005  # Thickness of CFRP in [m]

# For the polymer (here HDPE)
c_poly = 1880  # Specific heat capacity of HDPE in [J/kg/K]
Lambda_poly = 0.4  # HDPE thermal conductivity in [W/m/K]
h_conv_in = 230 # Outer forced convection coefficient in [W/m^2/k]
Rho_poly = 940  # Density in [kg/m^3]
e_poly = 0.020  # Thickness of HDPE in [m]

# From : Performance of Hydrogen Storage Tanks of Type IV in a Fire: Effect of the State of Charge
# h_out from : "Thermodynamic analysis of the emptying process of compressed hydrogen tanks" @ t = 200s

# For the reservoir : Tank 0
R0_in = 0.282  # Inner radius [m]
R0_out = R0_in + e_poly + e_CFRP  # External radius in [m]
L0 = 2  # Height in [m]
V0 = math.pi * pow(R0_in, 2) * L0   # [m^3]
m0_poly = Rho_poly * L0 * math.pi * (pow(R0_in + e_poly, 2)-pow(R0_in, 2))   # [kg]
m0_CFRP = Rho_CFRP * L0 * math.pi * (pow(R0_out, 2)-pow(R0_in + e_poly, 2))   # [kg]


UA0_in = 1 / ((1/(2*math.pi*R0_in*L0*h_conv_in))+(np.log((R0_in + e_poly)/R0_in)/(2*math.pi*L0*Lambda_poly)))  # [W/K]
UA0_out = 1 / ((1/(2*math.pi*R0_out*L0*h_conv_out))+(np.log(R0_out/(R0_in+e_poly))/(2*math.pi*L0*Lambda_CFRP)))  # [W/K]
UA0_int = 1 / ((np.log((R0_in + e_poly)/R0_in)/(2*math.pi*L0*Lambda_poly)) +
               (np.log(R0_out/(R0_in+e_poly))/(2*math.pi*L0*Lambda_CFRP)))  # [W/K]

# For the tanks to fill : Tank 1, Tank 2, Tank 3
# From: "The role of initial tank temperature on refuelling of on-board hydrogen tanks"
R_in = 0.115  # Inner radius [m]
R_out = R_in + e_poly + e_CFRP  # External radius in [m]
L = 0.827  # Height in [m]
V = math.pi * pow(R_in, 2) * L   # [m^3]
V1 = V
V2 = V
V3 = V
m_poly = Rho_poly * L * math.pi * (pow(R_in + e_poly, 2)-pow(R_in, 2))    # [kg]
m_CFRP = Rho_CFRP * L * math.pi * (pow(R_out, 2)-pow(R_in + e_poly, 2))   # [kg]

UA_in = 1 / ((1/(2*math.pi*R_in*L*h_conv_in))+(np.log((R_in + e_poly)/R_in)/(2*math.pi*L*Lambda_poly)))  # [W/K]
UA_out = 1 / ((1/(2*math.pi*R_out*L*h_conv_out))+(np.log(R_out/(R_in+e_poly))/(2*math.pi*L*Lambda_CFRP)))  # [W/K]
UA_int = 1 / ((np.log((R_in + e_poly)/R_in)/(2*math.pi*L*Lambda_poly)) +
              (np.log(R_out/(R_in+e_poly))/(2*math.pi*L*Lambda_CFRP)))  # [W/K]


# In the pipes:
# R_in = 0.005 # internal Radius of pipe before the valve [m]
# R_in = 0.003 # internal Radius of pipe after the valve [m]
A_0 = math.pi * 0.005**2  # [m^2]
A_1 = math.pi * 0.003**2  # [m^2]
A_2 = A_1
A_3 = A_1


# In the Heat exchanger
# Material of the pipe inside = austenitic stainless steel
L_hex = 0.300  # Length [m]
e_pipe = 0.002  # thickness of tube [m]
k_pipe = 21.5  # Conduction coefficient of the pipe in [W/m/K]
R_hex_in = 0.005  # Internal radius of pipe [m]
R_hex_out = R_hex_in + e_pipe  # External radius of pipe [m]
h_hex_in = 3540  # Internal convection heat transfer coefficient from calculation in pipe loss [W/m^2K]
k_hex = 21.5  # Conduction coefficient of the pipe in [W/m/K]
m_dot_H2O = 0.0156  # Mass flow rate of water [kg/s]

# To compute the external convection heat transfer coefficient
Prandtl_H2O = PropsSI('PRANDTL', 'P', [k*1.01325e5 for k in [1,200,300,500,800]], 'T',25+273.15, 'Water').mean()
Viscosity_H2O = PropsSI('VISCOSITY', 'P', [k*1.01325e5 for k in [1, 200, 300, 500, 800]], 'T', 25+273.15, 'Water').mean()
Re_H2O = 2 * m_dot_H2O / (math.pi*R_hex_out*Viscosity_H2O)
k_H2O = PropsSI('CONDUCTIVITY', 'P', [k*1.01325e5 for k in [1, 200, 300, 500, 800]], 'T', 25+273.15, 'H2').mean()  # [W/m/K]
Nu = conv_external.Nu_cylinder_Churchill_Bernstein(Re_H2O, Prandtl_H2O)  # Nusselt number

h_hex_out = Nu * k_H2O/(2*R_hex_out)  # in [W/m^2/K]
UA_hex = 1/(1/(h_hex_in*math.pi*2*R_hex_in*L_hex) + h_hex_in*(R_hex_out/R_hex_in)/(2*k_pipe*math.pi*L_hex)
            + 1/(h_hex_out*math.pi*2*R_hex_out*L_hex))  # [W/K]
cp_H2O = PropsSI('Cpmass', 'P', [k*1.01325e5 for k in [1, 200, 300, 500, 800]], 'T', 25+273.15, 'Water').mean()  # [J/kg/K]


# Pressure drop coefficient
C_D10 = -1
C_D20 = 1
C_D30 = 1
C_D11 = 1
C_D22 = 1
C_D33 = 1


# For the equations of states:
r = 4124.2  # [J/kgK]
a = 3.6E4  # [J K^(1/2)/ kg^2 / m^3]
b = 9.2E-3  # [m^3 / kg]


