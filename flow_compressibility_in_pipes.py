# Module importation
import numpy
from CoolProp.CoolProp import PropsSI
import math
import numpy as np
import matplotlib.pyplot as plt

# Constant
R = 0.005 # Radius of inlet pipe in [m]
A = math.pi * pow(R, 2) # Area of the pipe in [m^2]

# Computing the sound speed in hydrogen:
V_Tmin = PropsSI('A','P',[k*1.01325e5 for k in [1,200,300,500,800]],'T',25+273.15,'H2')
V_Tmax = PropsSI('A','P',[k*1.01325e5 for k in [1,200,300,500,800]],'T',80+273.15,'H2')

# Computing the Hydrogen density for the multiple pressure:
Rho_Tmin = PropsSI('D','P',[k*1.01325e5 for k in [1,200,300,500,800]],'T',25+273.15,'H2') # en [kg/m^3]
Rho_Tmax = PropsSI('D','P',[k*1.01325e5 for k in [1,200,300,500,800]],'T',80+273.15,'H2')

# Computing the velocity of Hydrogen:
MFR = numpy.linspace(0, 0.6, 20) # Mass Flow rate in [kg/s]

# We look for the Mach number at multiple presure and Temperature = 273 + 25 K
Ma_Tmin = np.zeros((len(MFR), 5))
Ma_Tmax = np.zeros((len(MFR), 5))

for c in range(0, 5):
    for l in range(0,len(MFR)):
        Ma_Tmin[l, c] = (MFR[l]/(Rho_Tmin[c]*A))/V_Tmin[c]

# print(V_Tmax, V_Tmin)
# we look for the Mach number at multiple presure and Temperature = 273 + 80 K

for c in range(0, 5):
    for l in range(0,len(MFR)):
        Ma_Tmax[l, c] = (MFR[l]/(Rho_Tmax[c]*A))/V_Tmax[c]



# Plot
l1, = plt.plot(MFR, Ma_Tmin[:, 1], 'b--')
l2, = plt.plot(MFR, Ma_Tmin[:, 2], 'g--')
l3, = plt.plot(MFR, Ma_Tmin[:, 3], 'c--')
l4, = plt.plot(MFR, Ma_Tmin[:, 4], 'm--')

l5, = plt.plot(MFR, Ma_Tmax[:, 1], 'b')
l6, = plt.plot(MFR, Ma_Tmax[:, 2], 'g')
l7, = plt.plot(MFR, Ma_Tmax[:, 3], 'c')
l8, = plt.plot(MFR, Ma_Tmax[:, 4], 'm')


legend = plt.legend((l1, l2, l3, l4, l5, l6, l7, l8), ["P = 200 atm, T= 25°C", "P = 300 atm, T= 25°C",
                    "P = 500 atm, T= 25°C", "P = 800 atm, T= 25°C", "P = 200 atm, T= 85°C",
         "P = 300 atm, T= 85°C", "P = 500 atm, T= 85°C", "P = 800 atm, T= 85°C"], loc="upper left")
plt.gca().add_artist(legend)

plt.ylabel('Mach number')
plt.xlabel('Mass flow rate [kg/s]')
plt.title('Variation of Mach number with mass flow rate')
plt.show()

# The flow is then considered as incompressible in the pipes by limiting the mass flow rate at 0.46 kg/s

