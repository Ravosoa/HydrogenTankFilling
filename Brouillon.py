
# For the equation of state
R = 8.3144598 # [J/mol/K]
M_H2 = 2.0159 # Molar mass [g/mol]
Rho_H2 = 0.08988 # Volumic mass [g/L]
r = R/M_H2
#Vm = (M_H2/Rho_H2)*pow(10,-3) # Molar volume [m^3/mol]
Vm = 33.2E-3
Tc_H2 = 33.15 # Temperature of critical point [K]
Pc_H2 = 13.0  # Pressure of critical point [bar]


def state_equation(T):
    # Redlich-Kwong equation of state for real gas
    a = 3.6E4  # [JK^1/2m^3/kg^2]
    b = 9.2E-3  # [m^3/kg]
    #a = (1 / (9 * (pow(2, 3 / 2) - 1))) * (pow(R, 2) * pow(Tc_H2, 2.5)) / Pc_H2
    #b = ((pow(2, 3 / 2) - 1) / 3) * (R * Tc_H2) / Pc_H2
    #print(r*T/(Vm-b))
    #print(a/(pow(T,0.5)*Vm*(Vm+b)))
    #print((Vm*(Vm+b)))
    P = r*T/(Vm-b) - a/(pow(T,0.5)*Vm*(Vm+b))
    return P*pow(10,-5)

#print(state_equation(40 +273))
# 50 MPa


# RK-4 method python program

# function to be solved
def f(x, y):
    return x + y


# or
# f = lambda x: x+y