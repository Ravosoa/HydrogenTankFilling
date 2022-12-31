from Initial_conditions import *
import math
import numpy as np
from CoolProp.CoolProp import PropsSI
from ht import *
from Parameters import *

# Parameters for the heat exchanger (changes with m_dot_water:

# Material of the pipe inside = austenitic stainless steel
L_hex = 0.300  # Length [m]
e_pipe = 0.002  # thickness of tube [m]
k_pipe = 21.5  # Conduction coefficient of the pipe in [W/m/K]
R_hex_in = 0.005  # Internal radius of pipe [m]
R_hex_out = R_hex_in + e_pipe  # External radius of pipe [m]
k_hex = 21.5  # Conduction coefficient of the pipe in [W/m/K]
m_dot_H2O_1 = 0.005  # Mass flow rate of water in hex 1 [kg/s]
m_dot_H2O_2 = 0.005  # Mass flow rate of water in hex 2[kg/s]
m_dot_H2O_3 = 0.005   # Mass flow rate of water in hex 3[kg/s]
m_dot_H2_1 = 0.005  # Mass flow rate of H2 in hex 1 [kg/s]
m_dot_H2_2 = 0.005  # Mass flow rate of H2 in hex 2[kg/s]
m_dot_H2_3 = 0.005   # Mass flow rate of H2 in hex 3[kg/s]

k_H2 = PropsSI('CONDUCTIVITY','P',[k*1.01325e5 for k in [1,200,300,500,800]],'T',25+273.15,'H2').mean() # en [W/m/K]


# External heat transfer computation
# To compute the external convection heat transfer coefficient
Prandtl_H2O = PropsSI('PRANDTL', 'P', [k*1.01325e5 for k in [1,200,300,500,800]], 'T',25+273.15, 'Water').mean()
k_H2O = PropsSI('CONDUCTIVITY', 'P', [k*1.01325e5 for k in [1, 200, 300, 500, 800]], 'T', 25+273.15, 'H2').mean()  # [W/m/K]
Viscosity_H2O = PropsSI('VISCOSITY', 'P', [k*1.01325e5 for k in [1, 200, 300, 500, 800]], 'T', 25+273.15, 'Water').mean()
cp_H2O = PropsSI('Cpmass', 'P', [k*1.01325e5 for k in [1, 200, 300, 500, 800]], 'T', 25+273.15, 'Water').mean()  # [J/kg/K]
Viscosity = PropsSI('VISCOSITY','P',[k*1.01325e5 for k in [1,200,300,500,800]],'T',25+273.15,'H2').mean()  # en [kg/m^3]
Prandtl = PropsSI('PRANDTL','P',[k*1.01325e5 for k in [1,200,300,500,800]],'T',25+273.15,'H2').mean() # en [kg/m^3]

if m_dot_H2O_1 == 0:
    UA_hex_1 = 0
else:
    Re_H2O_1 = 2 * m_dot_H2O_1 / (math.pi * R_hex_out * Viscosity_H2O)
    Nu_1 = conv_external.Nu_cylinder_Churchill_Bernstein(Re_H2O_1, Prandtl_H2O)
    h_hex_out_1 = Nu_1 * k_H2O / (2 * R_hex_out)  # in [W/m^2/K]
    Re_H2_1 = 2 * m_dot_H2_1 / (math.pi * R_in * Viscosity)
    f_1 = pow(0.790 * np.log(Re_H2_1) - 1.64, -2)
    Nu_H2_1 = (f_1 / 8) * (Re_H2_1 - 1000) * Prandtl / (1 + 12.7 * pow(f_1 / 8, 0.5) * (pow(Prandtl, 2 / 3) - 1))
    h_hex_in_1 = k_H2 * Nu_H2_1 / (2 * R_in)  # Internal convection heat transfer
                                            # coefficient from calculation in pipe loss [W/m^2K]
    UA_hex_1 = 1 / (1 / (h_hex_in_1 * math.pi * 2 * R_hex_in * L_hex) + (R_hex_out / R_hex_in) /
                    (2 * k_pipe * math.pi * L_hex) + 1 / (h_hex_out_1 * math.pi * 2 * R_hex_out * L_hex))  # [W/K]



if m_dot_H2O_2 == 0:
    UA_hex_2 = 0
else:
    Re_H2O_2 = 2 * m_dot_H2O_2 / (math.pi * R_hex_out * Viscosity_H2O)
    Nu_2 = conv_external.Nu_cylinder_Churchill_Bernstein(Re_H2O_2, Prandtl_H2O)
    h_hex_out_2 = Nu_2 * k_H2O / (2 * R_hex_out)  # in [W/m^2/K]
    Re_H2_2 = 2 * m_dot_H2_2 / (math.pi * R_in * Viscosity)
    f_2 = pow(0.790 * np.log(Re_H2_2) - 1.64, -2)
    Nu_H2_2 = (f_2 / 8) * (Re_H2_2 - 1000) * Prandtl / (1 + 12.7 * pow(f_2 / 8, 0.5) * (pow(Prandtl, 2 / 3) - 1))
    h_hex_in_2 = k_H2 * Nu_H2_2 / (2 * R_in)  # Internal convection heat transfer
                                                # coefficient from calculation in pipe loss [W/m^2K]
    UA_hex_2 = 1 / (1 / (h_hex_in_2 * math.pi * 2 * R_hex_in * L_hex) + (R_hex_out / R_hex_in) /
                    (2 * k_pipe * math.pi * L_hex) + 1 / (h_hex_out_2 * math.pi * 2 * R_hex_out * L_hex))  # [W/K]


if m_dot_H2O_3 == 0:
    UA_hex_3 = 0
else:
    Re_H2O_3 = 2 * m_dot_H2O_3 / (math.pi*R_hex_out*Viscosity_H2O)
    Nu_3 = conv_external.Nu_cylinder_Churchill_Bernstein(Re_H2O_3, Prandtl_H2O)
    h_hex_out_3 = Nu_3 * k_H2O/(2*R_hex_out)  # in [W/m^2/K]
    Re_H2_3 = 2 * m_dot_H2_3 /(math.pi*R_in*Viscosity)
    f_3 = pow(0.790*np.log(Re_H2_3)-1.64, -2)
    Nu_H2_3 = (f_3/8)*(Re_H2_3-1000)*Prandtl / (1 + 12.7*pow(f_3/8, 0.5)*(pow(Prandtl, 2/3)-1))
    h_hex_in_3 = k_H2*Nu_H2_3/(2*R_in)  # Internal convection heat transfer
                                        # coefficient from calculation in pipe loss [W/m^2K]
    UA_hex_3 = 1/(1/(h_hex_in_3*math.pi*2*R_hex_in*L_hex) + (R_hex_out/R_hex_in)/(2*k_pipe*math.pi*L_hex)
                + 1/(h_hex_out_3*math.pi*2*R_hex_out*L_hex))  # [W/K]


# print(UA_hex_1*(300-298))
# print(m_dot_H2O_1*cp_H2O*(273-280))
# print(UA_hex_1)
# print(f'Temps de résidence du HEX 1 :\n     {L_hex*2*math.pi*R_hex_in*rho_0i/u0[0]} [s]')
