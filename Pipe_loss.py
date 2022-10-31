# Module importation
import numpy
from CoolProp.CoolProp import PropsSI
import math
import numpy as np
import matplotlib.pyplot as plt

# Constant
R_in = 0.005 # internal Radius of inlet pipe in [m]
R_out = R_in + 0.002 # external Radius of inlet pipe in [m]
            # Data from: MODELING THE TRANSIENT TEMPERATURE DISTRIBUTION WITHIN A HYDROGEN CYLINDER DURING REFUELING
k = 21.5 # Conduction coefficient of the pipe in [W/m/K]
h_out = 1.79 # [W/m^2/K]
L_max = 30 # Length of the pipe in [m]
T_amb = 25 + 273 # [K]
MFR_max = 0.005 # Mass flow rate in [kg/s]
T_max = 80 + 273 # [K]
#https://reader.elsevier.com/reader/sd/pii/S2238785422006433?token=02C23CA72A2D1285B9AC9465740E3FF5B2954E1C36EF3FC1AAF43
# 6BAE0B43EC2E01A83C6E2585E11FBD9A241998CC43B&originRegion=eu-west-1&originCreation=20221018074136

# Computing the Hydrogen density for the multiple pressure:
Rho_Tmin = PropsSI('D','P',[k*1.01325e5 for k in [1,200,300,500,800]],'T',25+273.15,'H2') # en [kg/m^3]


#Rho_Tmax = PropsSI('D','P',[k*1.01325e5 for k in [1,200,300,500,800]],'T',80+273.15,'H2')
#print ("Rho_Tmax = ", Rho_Tmax)


# Get the fluid viscosity:
Viscosity = PropsSI('VISCOSITY','P',[k*1.01325e5 for k in [1,200,300,500,800]],'T',25+273.15,'H2') # en [kg/m^3]
Viscosity_min = np.min(Viscosity)
print("Viscosity min = ", Viscosity_min) # in [Pa s]



Re_max = 2 * MFR_max /(math.pi*R_in*Viscosity_min)
print(f"the max of Reynold's number is : \n Re_max = {Re_max} ")

# Get the Prandtl number:
Prandtl = PropsSI('PRANDTL','P',[k*1.01325e5 for k in [1,200,300,500,800]],'T',25+273.15,'H2') # en [kg/m^3]
Pr_max = np.max(Prandtl)
Pr_min = np.min(Prandtl)
print(f"Prandtl_max = {Pr_max} \nPrandtl_min = {Pr_min} ") # in [Pa s]

# Computing Nusselt number:
    # The friction factor is:
f = pow(0.790*np.log(Re_max)-1.64, -2)

Nu_max = (f/8)*(Re_max-1000)*Pr_max / (1 + 12.7*pow(f/8, 0.5)*(pow(Pr_max, 2/3)-1))
print(f"The Nu_max = {Nu_max}")

# Fluid conductivity:
k_fluid = PropsSI('CONDUCTIVITY','P',[k*1.01325e5 for k in [1,200,300,500,800]],'T',25+273.15,'H2') # en [W/m/K]
k_fluid_max = np.max(k_fluid) # in [W/m/K]

print(f"k_fluid_max = {k_fluid_max} [W/m/K]")


# Getting the internal convection coefficient:
h_in_max = k_fluid_max*Nu_max/(2*R_in)

print(f"h_in_max = {h_in_max} [W/m^2/K]")


R_conv_in = 1/(2 * math.pi * R_in * L_max * h_in_max)
R_conv_out = 1/(2 * math.pi * R_out * L_max * h_out)
R_cond = np.log(R_out/R_in)/(2 * math.pi * L_max * k)

R_tot = R_conv_in + R_cond + R_conv_out
print(f"R_tot = {R_tot} [K /W]")

Q_loss_max = (T_max-T_amb)/R_tot

print(f"The heat loss through pipes bu convection and conduction is: \nQ_loss_max ={Q_loss_max} [W]")
# It is very small compared to other heat exchanges in the system that are in the order of magnitude of [kW]
# Hence, heat losses through conduction and convection are neglected
