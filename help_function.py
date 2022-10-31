# This doc will be used to define more easily the differential equation to be solved with RK4
from Parameters import *
import numpy as np


def h(T, rho):
    # Enthalpy for real gaz with Redlich-Kwong EoS
    # t : Temperature
    # rho: Density
    a = 3.6E4
    b = 9.2E-3
    r = 4124.2  # [J/kgK]
    nu = 1/rho
    return C_v*T + 3*a*np.log(nu/(nu+b))/(2*b*np.sqrt(T)) + r*T*nu/(nu-b)-a/((nu+b)*np.sqrt(T))


def Q(ua, T1, T2):
    # Heat loss through the tank
    return ua*(T2-T1)


def v(A, rho, dm):
    # compute the velocity
    return dm/(A*rho)

def P(T, rho):
    a = 3.6E4
    b = 9.2E-3
    r = 4124.2  # [J/kgK]
    nu = 1/rho
    return r*T/(nu-b) - a/(np.sqrt(T)*nu*(nu+b))

def Q_hex_H2O (dm,T_in, T_out):
    cp =
    return dm*cp*(T_out-T_in)