# Hydrogen fast and safe filling control

EPFL - Project in Energy I, Autumn 2022 

Supervisor:	Prof. Haussener Sophia

Mentor: Matthieu Jonin 

## Repository Overview
- `Parameters.py` - Provide the physical parameters of the system
- `Initial_conditions.py` - Initial conditions to run the optimisation problem
- `main.py` - Controller optimisation problem
- `Simulation_casadi.py` - Optimisation problem without optimising the controllers
- `System_constraints.py` - Constraints of the optimisation problem
- `Plot_figure.py` - Plot the system state and controllers

The following files are solutions of the optimisation with multiple targets. They can be used as the initial condition of the 
optimisation problem:
- `sol_no_hex.pkl` - All heat exchangers are turned off
- `sol_target_1kg.pkl` - Target mass of the 3 tanks is 1kg
- `sol_taget_1_05_005.pkl` - Three targets: 1kg, 0.5kg, 0kg

The following files are used to verify some assumptions:
- `Flow_compressibility_in_pipes.py` - Compute Mach number
- `Pipe_loss.py` - Compute the heat loss through pipes

## Abstract
In this work, we present a system for filling hydrogen tanks based on thermodynamic principles and include a control 
system to ensure safety and efficiency. The system consists of a large hydrogen reservoir tank and three smaller tanks 
that need to be filled.

We first validate our model by comparing it to a previously published model and show that it 
accurately predicts the filling process. We then investigate the use of a heat exchanger as a way to control the 
temperature of the hydrogen gas and improve the efficiency of the filling process.

Our analysis demonstrates that the proposed hydrogen storage and the filling system are feasible but demanding in terms 
of power consumption.

## Instructions

If a basic simulation of the system is needed, `Simulation_casadi.py` file is sufficient to run.

To perform an optimal control simulation of the system, `main.py` 
should be run. It is important to pay attention to the initial conditions as the IPopt solver is sensitive to them. 
Some feasible initial conditions are provided in the file.

To run a simulation using the control system:
- Set the desired goal in `System_constraints.py` in the `inequality_objective` function 
(Note that the minimum mass in the tank is 0.056kg)
- Choose the initial guess `m0_s` in `main.py`
- Run `main.py`


### Steps

## Environment
Version: `Python 3.10`

Libraries: `numpy`, `matplotlib`, `scipy`, `casadi`, `CoolProp`
