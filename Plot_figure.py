import numpy as np
import matplotlib.pyplot as plt
from System_constraints import load_sol


def plot_solution(sol):
    """Plot the mass, temperature, pressure in tank 1, 2, 3.
       Plot the controllers: masse flow rate and power of heat exchanger.

    Note: Comment the lines for figure 4 and 5 when the controllers are not optimise. An error will happen if not.

        Args:
            sol: solution of the optimization problem to plot (".pkl" file)

        Returns:
            Plots
    """
    # Take the data from the solution
    z = sol['z']
    x = sol['x']
    u = sol['u']
    N = sol['num_iteration']
    t = np.arange(0, N * sol['h'], sol['h'])
    # Because h is not always dividing N at the right number, sometimes it is necessary to remove one step
    if len(t) != N:
        t = t[:-1]

    print(f"The optimal time step h is {sol['h']}")

    plt.figure(1)
    plt.plot(t, x[5], 'g', label=r'$m_1$')
    plt.plot(t, x[6], 'c--', label=r'$m_2$')
    plt.plot(t, x[7], 'm--', label=r'$m_3$')
    plt.axhline(1, color='r', label='Target ')
    #plt.axhline(0.5, color='r', label='Target 2')
    plt.ylabel('Mass [kg]')
    plt.xlabel('Time [s]')
    plt.title('Mass evolution in tanks with time')
    plt.legend(loc='upper left')

    plt.figure(2)
    #plt.plot(t, z[17], 'b--', label='T_0')
    plt.plot(t, z[18] - 273, 'g', label=r'$T_1$')
    plt.plot(t, z[19] - 273, 'c--', label=r'$T_2$')
    plt.plot(t, z[20] - 273, 'm--', label=r'$T_3$')
    plt.axhline(85, color='r', label='Limiting temperature')
    plt.ylabel('Temperature [°C]')
    #plt.ylim(0, 90)
    plt.xlabel('Time [s]')
    plt.title('Temperature evolution in tanks with time')
    plt.legend(loc='lower right')

    plt.figure(3)
    #plt.plot(t, z[21] / (10 ** 5), 'b--', label='P_0')
    plt.plot(t, z[22] / (10 ** 5), 'g', label=r'$P_1$')
    plt.plot(t, z[23] / (10 ** 5), 'c--', label=r'$P_2$')
    plt.plot(t, z[24] / (10 ** 5), 'm--', label=r'$P_3$')
    plt.axhline(1.25 * 70 * 10, color='r', label=r'$P_{lim}$')
    plt.title('Pressure evolution in tanks with time')
    plt.ylabel('Pressure [bar]')
    plt.xlabel('Time [s]')
    plt.legend(loc='upper right')

    # Comment the lines for figure 4 and 5 if not optimizing the controllers
    plt.figure(4)
    plt.plot(t, u[0], 'g', label=r'$\dot{m}_1$')
    plt.plot(t, u[1], 'c--', label=r'$\dot{m}_2$')
    plt.plot(t, u[2], 'm--', label=r'$\dot{m}_3$')
    plt.ylabel(r'$\dot{m}$ [kg/s]')
    plt.xlabel('Time [s]')
    plt.legend(loc='upper right')
    plt.title('Controlled mass flow rate with time')

    plt.figure(5)
    plt.plot(t, u[3], 'g', label=r'$\dot{Q}_{HEX, 1}$')
    plt.plot(t, u[4], 'c--', label=r'$\dot{Q}_{HEX, 2}$')
    plt.plot(t, u[5], 'm--', label=r'$\dot{Q}_{HEX, 3}$')
    plt.ylabel(r'$\dot{Q}_{HEX}$ [W]')
    plt.xlabel('Time [s]')
    plt.legend(loc='lower left')
    plt.title('Controlled heat exchanged with time')

    plt.show()


if __name__ == "__main__":
    sol = load_sol('sol_tmp')
    #sol = load_sol("sol_target_1_05_005")
    #sol = load_sol("sol_no_hex")
    plot_solution(sol)
