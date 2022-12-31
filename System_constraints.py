# Module importation
from Parameters import parameters as p
import pickle
from casadi import *
from Initial_conditions import P_max, T_max, mass_max, MFR_max, Bound_h, Q_max,masse_0


def normalize(m, lb, ub):
    """Normalize the value so that it gets a value between 0 and 1.

        Args:
            m: the variable to normalize
            lb: lower bound (tuple)
            ub: upper bound (tuple)

        Returns:
            The normalized variable
    """
    return (m - lb) / (ub - lb)


def unormalize(ms, lb, ub):
    """Recover the original variable.

        Args:
            ms: normalized variable
            lb: lower bound (tuple)
            ub: upper bound (tuple)

        Returns:
            The recovered variable
    """
    return ms * (ub-lb) + lb


def equality_constraints(x, u, z):
    """Equalities to be satisfied to get the z-states of the system.
            Args:
                x: state variables (array)
                u: controllers (array)
                z: state variable to solve (array)

            Returns:
                List of equality constraints to solve
    """
    return vertcat(*[z[0] - (u[0] + u[1] + u[2]) / (x[16] * p['A_0']),  # z0 - v_0
                     z[1] - u[0] / (x[16] * p['A_1']),  # z1 - v_valve_1
                     z[2] - u[1] / (x[16] * p['A_2']),  # z2 - v_valve_2
                     z[3] - u[2] / (x[16] * p['A_3']),  # z3 - v_valve_3
                     z[4] - z[21],  # z4 P_valve_1
                     z[5] - z[21],  # z5 P_valve_2
                     z[6] - z[21],  # z6 P_valve_3
                     p['r'] * (z[7]) / (1 / x[16] - p['b']) - p['a'] / (sqrt(z[7]) * (1 / x[16]) * ((1 / x[16]) + p['b'])) - z[4],  # z7 T_valve_1
                     p['r'] * (z[8]) / (1 / x[16] - p['b']) - p['a'] / (sqrt(z[8]) * (1 / x[16]) * ((1 / x[16]) + p['b'])) - z[5],  # z8 T_valve_2
                     p['r'] * (z[9]) / (1 / x[16] - p['b']) - p['a'] / (sqrt(z[9]) * (1 / x[16]) * ((1 / x[16]) + p['b'])) - z[6],  # z9 T_valve_3
                     z[10] - p['C_v'] * z[7] - (3 * p['a'] / (2 * p['b'] * sqrt(z[7]))) * log(1 / (1 + p['b'] * x[16])) - z[4] / x[16],  # z10 h_valve_1
                     z[11] - p['C_v'] * z[8] - (3 * p['a'] / (2 * p['b'] * sqrt(z[8]))) * log(1 / (1 + p['b'] * x[16])) - z[5] / x[16],  # z11 h_valve_2
                     z[12] - p['C_v'] * z[9] - (3 * p['a'] / (2 * p['b'] * sqrt(z[9]))) * log(1 / (1 + p['b'] * x[16])) - z[6] / x[16],  # z12 h_valve_3
                     (z[13] - z[7]) * (u[0] * p['C_p']) - u[3],  # z13 T_hex_1
                     (z[14] - z[8]) * (u[1] * p['C_p']) - u[4],  # z14 T_hex_2
                     (z[15] - z[9]) * (u[2] * p['C_p']) - u[5],  # z15 T_hex_3
                     z[16] - (x[0] + p['r'] * (z[17]) / (1 - p['b'] * x[16]) - p['a'] / (sqrt(z[17]) * (1 / x[16] + p['b']))),  # z16 h_0
                     p['C_v'] * z[17] - (3 * p['a'] / (2 * p['b'] * sqrt(z[17]))) * log(1 + p['b'] * x[16]) - x[0],  # z17 T_0
                     p['C_v'] * z[18] * sqrt(z[18]) - (3 * p['a'] / (2 * p['b'])) * log(1 + p['b'] * x[17]) - x[1] * sqrt(z[18]),  # z18 T_1
                     p['C_v'] * z[19] * sqrt(z[19]) - (3 * p['a'] / (2 * p['b'])) * log(1 + p['b'] * x[18]) - x[2] * sqrt(z[19]),  # z19 T_2
                     p['C_v'] * z[20] * sqrt(z[20]) - (3 * p['a'] / (2 * p['b'])) * log(1 + p['b'] * x[19]) - x[3] * sqrt(z[20]),  # z20 T_3
                     z[21] - (p['r'] * z[17] / ((1 / x[16]) - p['b']) - p['a'] / (sqrt(z[17]) * (1 / x[16]) * (p['b'] + 1 / x[16]))),  # z21 P_0
                     z[22] - (p['r'] * z[18] / ((1 / x[17]) - p['b']) - p['a'] / (sqrt(z[18]) * (1 / x[17]) * (p['b'] + 1 / x[17]))),  # z22 P_1
                     z[23] - (p['r'] * z[19] / ((1 / x[18]) - p['b']) - p['a'] / (sqrt(z[19]) * (1 / x[18]) * (p['b'] + 1 / x[18]))),  # z23 P_2
                     z[24] - (p['r'] * z[20] / ((1 / x[19]) - p['b']) - p['a'] / (sqrt(z[20]) * (1 / x[19]) * (p['b'] + 1 / x[19]))),  # z24 P_3
                     ])


def inequality_constraints(u, z, h):
    """Inequalities to be satisfied to get a safe filling process. THe inequalities have the format: g(x) <= 0
                Args:
                    h: time step of the shooting interval (float)
                    u: controllers (array)
                    z: state variable to solve (array)

                Returns:
                    List of inequality constraints to satisfy
    """
    m_dot_maximum = 0.1
    P_lim_sup = 1.25 * 70 * 10**6  # MPa
    T_lim_sup = 273 + 85  # K
    Q_lim_sup = 60000
    return vertcat(*[(z[22] - P_lim_sup)/P_max,
                     (z[23] - P_lim_sup)/P_max,
                     (z[24] - P_lim_sup)/P_max,
                     (z[18] - T_lim_sup)/T_max,
                     (z[19] - T_lim_sup)/T_max,
                     (z[20] - T_lim_sup)/T_max,
                     (h - 0.8)/Bound_h[1],
                     (0.01 - h)/Bound_h[1],
                     # (u[3] - Q_lim_sup) / Q_max,
                     # (u[4] - Q_lim_sup) / Q_max,
                     # (u[5] - Q_lim_sup) / Q_max,
                     (u[0] - m_dot_maximum)/MFR_max,
                     (u[1] - m_dot_maximum)/MFR_max,
                     (u[2] - m_dot_maximum)/MFR_max,
                     ])


def inequality_objectif(x, z):
    """Terminal constraints to satisfy to get the final mass in the tank.
                Args:
                    x: state variables (array)
                    z: state variable  (array)

                Returns:
                    List of inequality terminal constraints
    """
    mass_lim_inf1 = 1
    mass_lim_inf2 = 0.5
    mass_lim_inf3 = 0.05
    mass_lim_sup1 = 1.002
    mass_lim_sup2 = 0.52
    mass_lim_sup3 = 0.056
    return vertcat(*[(mass_lim_inf1 - x[5])/mass_max,
                     (mass_lim_inf2 - x[6])/mass_max,
                     (mass_lim_inf3 - x[7])/mass_max,
                     (-mass_lim_sup1 + x[5]) / mass_max,
                     (-mass_lim_sup2 + x[6]) / mass_max,
                     (-mass_lim_sup3 + x[7]) / mass_max,
                     ])


def LotkaVolterraModel(x, u, z):
    """Equalites in the state-space representation: x_dot= f(x)
            Args:
                x: state variables to solve (array)
                u: controllers (array)
                z: state variables (array)

            Returns:
                List of equality constraints
    """
    xdot = vertcat(*[1/x[4] * (-(u[0]+u[1]+u[2])*(z[16] + ((-(u[0]+u[1]+u[2])/(p['A_0']*x[16]))**2)/2) + p['UA0_in']*(x[12] - (z[17]))
                               + x[0]*(u[0]+u[1]+u[2])),                                                      # x0 u_0
                     1/x[5] * (u[0]*z[10] - u[3] + u[0]*(z[1]**2)/2 + p['UA_in']*(x[13]-(z[18])) - x[1]*u[0]),    # x1 u_1
                     1/x[6] * (u[1]*z[11] - u[4] + u[1]*(z[2]**2)/2 + p['UA_in']*(x[14]-(z[19])) - x[2]*u[1]),   # x u_2
                     1/x[7] * (u[2]*z[12] - u[5] + u[2]*(z[3]**2)/2 + p['UA_in']*(x[15]-(z[20])) - x[3]*u[2]),  # x3 u_3
                     -(u[0] + u[1] + u[2]),                                                             # x4 m_0
                     u[0],                                                                              # x5 m_1
                     u[1],                                                                              # x6 m_2
                     u[2],                                                                              # x7 m_3
                     1 / (p['m0_CFRP'] * p['c_CFRP']) * (p['UA0_int'] * (x[12] - x[8]) - p['UA0_out'] * (x[8] - p['T_amb'])),    # x8 T_CFRP_0
                     1 / (p['m_CFRP'] * p['c_CFRP']) * (p['UA_int'] * (x[13] - x[9]) - p['UA_out'] * (x[9] - p['T_amb'])),      # x9 T_CFRP_1
                     1 / (p['m_CFRP'] * p['c_CFRP']) * (p['UA_int'] * (x[14] - x[10]) - p['UA_out'] * (x[10] - p['T_amb'])),    # x10 T_CFRP_2
                     1 / (p['m_CFRP'] * p['c_CFRP']) * (p['UA_int'] * (x[15] - x[11]) - p['UA_out'] * (x[11] - p['T_amb'])),    # x11 T_CFRP_3
                     1 / (p['m0_poly'] * p['c_poly']) * (p['UA0_in'] * ((z[17]) - x[12]) - p['UA_int'] * (x[12] - x[8])),  # x12 T_poly_0
                     1 / (p['m_poly'] * p['c_poly']) * (p['UA0_in'] * ((z[18]) - x[13]) - p['UA_int'] * (x[13] - x[9])),    # x13 T_poly_1
                     1 / (p['m_poly'] * p['c_poly']) * (p['UA0_in'] * ((z[19]) - x[14]) - p['UA_int'] * (x[14] - x[10])),   # x14 T_poly_2
                     1 / (p['m_poly'] * p['c_poly']) * (p['UA0_in'] * ((z[20]) - x[15]) - p['UA_int'] * (x[15] - x[11])),   # x15 T_poly_3
                     -(u[0] + u[1] + u[2]) / (p['V0']),                                                     # x16 rho_0
                     u[0] / p['V1'],                                                                        # x17 rho_1
                     u[1] / p['V2'],                                                                        # x18 rho_2
                     u[2] / p['V3'] ])                                                                       # x19 rho_3
    return xdot


def save_sol(savefilename, sol):
    """Save the solution file in ".pkl" format.
            Args:
                savefilename: Name of file to save (string)
                sol: the file to save

            Returns:
                -
    """
    with open(savefilename+'.pkl', 'wb') as handle:
        pickle.dump(sol, handle, protocol=pickle.HIGHEST_PROTOCOL)


def load_sol(savefilename):
    """Load the solution file in ".pkl" format.
            Args:
                savefilename: Name of file to open (string)

            Returns:
                -
    """
    with open(savefilename+'.pkl', 'rb') as handle:
        sol = pickle.load(handle)
    return sol


