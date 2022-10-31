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
#print (V_Tmin)
#array([1315.38841602, 1491.63249042, 1581.48572567, 1756.26291713, 2000.44266183]) # en [m/s]
V_Tmax = PropsSI('A','P',[k*1.01325e5 for k in [1,200,300,500,800]],'T',80+273.15,'H2') # k [atm] changé en [Pa]
#print(V_Tmax)
#array([1428.8541971 , 1590.36217432, 1671.38859285, 1829.0351964 , 2051.07153236]) # en [m/s]

# Computing the Hydrogen density for the multiple pressure:
Rho_Tmin = PropsSI('D','P',[k*1.01325e5 for k in [1,200,300,500,800]],'T',25+273.15,'H2') # en [kg/m^3]
#print ("Rho_Tmin = ", Rho_Tmin)

Rho_Tmax = PropsSI('D','P',[k*1.01325e5 for k in [1,200,300,500,800]],'T',80+273.15,'H2')
#print ("Rho_Tmax = ", Rho_Tmax)

# Computing the velocity of Hydrogen:
MFR = numpy.linspace(0, 0.02, 20) # Mass Flow rate in [kg/s]
#print('MFR = ', MFR)

# we look for the Mach number at multiple presure and Temperature = 273 + 25 K
Ma_Tmin = np.zeros((len(MFR), 5))

for c in range(0, 5):
    for l in range(0,len(MFR)):
        Ma_Tmin[l, c] = (MFR[l]/(Rho_Tmin[c]*A))/V_Tmin[c]

#print('Ma_Tmin = ', Ma_Tmin)

# # Plot
# l1, = plt.plot(MFR, Ma_Tmin[:, 0], 'k')
# l2, = plt.plot(MFR, Ma_Tmin[:, 1], 'b')
# l3, = plt.plot(MFR, Ma_Tmin[:, 2], 'g')
# l4, = plt.plot(MFR, Ma_Tmin[:, 3], 'c')
# l5, = plt.plot(MFR, Ma_Tmin[:, 4], 'm')
# l6, = plt.plot(MFR, np.ones(len(MFR))*0.3, 'r')
#
# legend = plt.legend((l1, l2, l3, l4, l5, l6), ["P = 1 [atm]", "P = 200 [atm]", "P = 300 [atm]", "P = 500 [atm]", "P = 800 [atm]", "Ma limit for incompressibility"], loc="upper right")
# plt.gca().add_artist(legend)
#
# plt.ylabel('Mach number')
# plt.xlabel('Mass flow rate in [kg/s]')
# plt.title('Mach number Vs MFR for T = 273 + 25 K')
# plt.show()


# we look for the Mach number at multiple presure and Temperature = 273 + 80 K
Ma_Tmax = np.zeros((len(MFR), 5))

for c in range(0, 5):
    for l in range(0,len(MFR)):
        Ma_Tmax[l, c] = (MFR[l]/(Rho_Tmax[c]*A))/V_Tmax[c]

#print('Ma_Tmax = ', Ma_Tmax)

# # Plot T = 80 °
# l1, = plt.plot(MFR, Ma_Tmax[:, 0], 'k')
# l2, = plt.plot(MFR, Ma_Tmax[:, 1], 'b')
# l3, = plt.plot(MFR, Ma_Tmax[:, 2], 'g')
# l4, = plt.plot(MFR, Ma_Tmax[:, 3], 'c')
# l5, = plt.plot(MFR, Ma_Tmax[:, 4], 'm')
# l6, = plt.plot(MFR, np.ones(len(MFR))*0.3, 'r')
#
#
# legend = plt.legend((l1, l2, l3, l4, l5, l6), ["P = 1 [atm]", "P = 200 [atm]", "P = 300 [atm]", "P = 500 [atm]", "P = 800 [atm]", "Ma limit for incompressibility"], loc="upper right")
# plt.gca().add_artist(legend)
#
# plt.ylabel('Mach number')
# plt.xlabel('Mass flow rate in [kg/s]')
# plt.title('Mach number Vs MFR for T = 273 + 80 K')
# plt.show()


# Plot T = 25°C and T = 80°C

#l1, = plt.plot(MFR, Ma_Tmin[:, 0], 'k--')
l2, = plt.plot(MFR, Ma_Tmin[:, 1], 'b--')
l3, = plt.plot(MFR, Ma_Tmin[:, 2], 'g--')
l4, = plt.plot(MFR, Ma_Tmin[:, 3], 'c--')
l5, = plt.plot(MFR, Ma_Tmin[:, 4], 'm--')

#l6, = plt.plot(MFR, Ma_Tmax[:, 0], 'k')
l7, = plt.plot(MFR, Ma_Tmax[:, 1], 'b')
l8, = plt.plot(MFR, Ma_Tmax[:, 2], 'g')
l9, = plt.plot(MFR, Ma_Tmax[:, 3], 'c')
l10, = plt.plot(MFR, Ma_Tmax[:, 4], 'm')
#l11, = plt.plot(MFR, np.ones(len(MFR))*0.3, 'r')


legend = plt.legend((l2, l3, l4, l5), ["P = 200 [atm]", "P = 300 [atm]", "P = 500 [atm]", "P = 800 [atm]"], loc="upper right")
plt.gca().add_artist(legend)

plt.ylabel('Mach number')
plt.xlabel('Mass flow rate in [kg/s]')
plt.title('Mach number Vs MFR for T = 273 + 80 K and T = 273 +25 K')
plt.show()
# The flow is then considered as incompressible in the pipes

