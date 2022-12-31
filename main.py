import casadi as cs
import time

# import functions
from System_constraints import equality_constraints, LotkaVolterraModel, unormalize, save_sol, load_sol, \
    inequality_constraints, inequality_objectif
from Initial_conditions import Bounds_x, Bounds_u, Bounds_z, Bound_h, z0, u0, x0, timestep_0
from Plot_figure import plot_solution

# parameters
N = 90
tol_cons_eq = 1e-6  # Tolerance on equality constraints
tol_cons_ineq = 1e-6  # Tolerance on inequalities
n_inequality = 9 + 2

nx = len(x0)
nu = len(u0)
nz = len(z0)

# extract bounds
lbz = [];
ubz = []
for i in range(len(Bounds_z)):
    lbz.append(Bounds_z[i][0])
    ubz.append(Bounds_z[i][1])
lbx = [];
ubx = []
for i in range(len(Bounds_x)):
    lbx.append(Bounds_x[i][0])
    ubx.append(Bounds_x[i][1])
lbu = [];
ubu = []
for i in range(len(Bounds_u)):
    lbu.append(Bounds_u[i][0])
    ubu.append(Bounds_u[i][1])

lb_time = [];
ub_time = []
lb_time.append(Bound_h[0])
ub_time.append(Bound_h[1])

lbz = cs.DM(lbz)
ubz = cs.DM(ubz)
lbx = cs.DM(lbx)
ubx = cs.DM(ubx)
lbu = cs.DM(lbu)
ubu = cs.DM(ubu)
lb_time = cs.DM(lb_time)
ub_time = cs.DM(ub_time)

# Concatenate DM vectors vertically
lbm = cs.vertcat(cs.vertcat(*[lbx] * N), cs.vertcat(*[lbz] * N), cs.vertcat(*[lbu] * N), lb_time)
ubm = cs.vertcat(cs.vertcat(*[ubx] * N), cs.vertcat(*[ubz] * N), cs.vertcat(*[ubu] * N), ub_time)

# symbolic variables
X = cs.SX.sym('x', nx, N)
X_vec = cs.reshape(X, nx * N, 1)

Z = cs.SX.sym('z', nz, N)
Z_vec = cs.reshape(Z, nz * N, 1)

U = cs.SX.sym('u', nu, N)
U_vec = cs.reshape(U, nu * N, 1)

h = cs.SX.sym('h', 1, 1)

m = cs.vertcat(X_vec, Z_vec, U_vec, h)

# build problem
H = [];
lbh = [];
ubh = []  # Equality constraints
G = [];
lbg = [];
ubg = []

for i in range(N):
    if i == 0:
        x_i = x0
    else:
        x_i = X[0:nx, i - 1]
    z_i = Z[0:nz, i]
    u_i = U[0:nu, i]

    # RK4 on x_i
    k1 = h * LotkaVolterraModel(x=x_i, u=u_i, z=z_i)
    k2 = h * LotkaVolterraModel(x_i + k1 / 2, u=u_i, z=z_i)
    k3 = h * LotkaVolterraModel(x_i + k2 / 2, u=u_i, z=z_i)
    k4 = h * LotkaVolterraModel(x_i + k3, u=u_i, z=z_i)
    dx = (k1 + 2 * k2 + 2 * k3 + k4) / 6
    H.append(x_i + dx - X[0:nx, i])  # For the multiple shooting
    lbh.append(-tol_cons_eq * cs.DM.ones(nx, 1))
    ubh.append(tol_cons_eq * cs.DM.ones(nx, 1))

    # Algebraic
    H.append(equality_constraints(x=x_i, u=u_i, z=z_i))  # add the equality constraints
    lbh.append(-tol_cons_eq * cs.DM.ones(nz, 1))
    ubh.append(tol_cons_eq * cs.DM.ones(nz, 1))

    G.append(inequality_constraints(u=u_i, z=z_i, h=h))  # add the inequality constraints
    lbg.append(-cs.inf * cs.DM.ones(n_inequality, 1))
    ubg.append(tol_cons_ineq * cs.DM.ones(n_inequality, 1))

    if i >= N - 3:
        G.append(inequality_objectif(x=X[0:nx, i], z=Z[0:nz, i]))  # add the terminal constraints
        lbg.append(-cs.inf * cs.DM.ones(6, 1))
        ubg.append(tol_cons_ineq * cs.DM.ones(6, 1))

# concatenate
G = cs.vertcat(*G)
H = cs.vertcat(*H)

cons = cs.Function('cons', [m], [cs.vertcat(G, H)])
lb_cons = cs.vertcat(cs.vertcat(*lbg), cs.vertcat(*lbh))
ub_cons = cs.vertcat(cs.vertcat(*ubg), cs.vertcat(*ubh))

f = cs.Function('f', [m], [h ** 2])

# def problem
m_s = cs.MX.sym('m_s', N * (nx + nz + nu) + 1, 1)  # scale variable
m_unsc = unormalize(ms=m_s, lb=lbm, ub=ubm)

# Load initial condition
#sol_init = load_sol("sol_target_1kg")
sol_init = load_sol("sol_target_1kg_05kg_005kg")
# sol_init = load_sol("sol_tmp")


# m0_s = normalize(cs.vertcat(sol_init['m_opt_s'], cs.vertcat(*[u0]*N), timestep_0), lb=lbm, ub=ubm)
m0_s = sol_init['m_opt_s']

nlp = {'x': m_s, 'f': f(m_unsc), 'g': cons(m_unsc)}
opt = {'ipopt': {'max_iter': 10000, 'acceptable_tol': tol_cons_ineq}}
S = cs.nlpsol('S', 'ipopt', nlp, opt)

print(S)

# solve the NLP
t_start = time.perf_counter()
res = S(x0=m0_s,
        lbx=cs.DM.zeros((nz + nx + nu) * N + 1, 1),
        ubx=cs.DM.ones((nz + nx + nu) * N + 1, 1),
        lbg=lb_cons,
        ubg=ub_cons
        )

t_end = time.perf_counter()

print(f"\n Solving time : {t_end - t_start} s")

# process
m_opt_s = res['x']
m_opt_unsc = unormalize(ms=m_opt_s, ub=ubm, lb=lbm)  # Need to unscale the optimal value

sol = {'x': cs.reshape(m_opt_unsc[0:nx * N], nx, N).toarray(),
       'z': cs.reshape(m_opt_unsc[nx * N:(nx + nz) * N], nz, N).toarray(),
       'u': cs.reshape(m_opt_unsc[(nx + nz) * N:(nx + nz + nu) * N], nu, N).toarray(),
       'obj': res['f'].toarray(),
       'h': cs.reshape(m_opt_unsc[(nx + nz + nu) * N:(nx + nz + nu) * N + 1], 1, 1).toarray(),
       'num_iteration': N,
       'm_opt_s': m_opt_s,
       }

time.sleep(1)

# save
save_sol('sol_tmp', sol)

# call plot
solution = load_sol('sol_tmp')
plot_solution(solution)

print('End')
