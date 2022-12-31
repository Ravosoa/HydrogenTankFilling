import casadi as cs
import time
from System_constraints import equality_constraints, LotkaVolterraModel, normalize, unormalize, save_sol, load_sol
from Initial_conditions import Bounds_x, Bounds_z, z0, u0, x0
from Plot_figure import plot_solution

# Parameters
N = 90  # Number of intervals
h = 0.3  # Time step
tol_cons_eq = 1e-6  # Tolerance
tol_cons_ineq = 1e-6
n_inequality = 9

nx = len(x0)
nu = len(u0)
nz = len(z0)

# Extract bounds
lbz = []; ubz = []
for i in range(len(Bounds_z)):
    lbz.append(Bounds_z[i][0])
    ubz.append(Bounds_z[i][1])

lbx = []; ubx = []
for i in range(len(Bounds_x)):
    lbx.append(Bounds_x[i][0])
    ubx.append(Bounds_x[i][1])


lbz = cs.DM(lbz)
ubz = cs.DM(ubz)
lbx = cs.DM(lbx)
ubx = cs.DM(ubx)

# Concatenate DM vectors vertically
lbm = cs.vertcat(cs.vertcat(*[lbx]*N), cs.vertcat(*[lbz]*N))
ubm = cs.vertcat(cs.vertcat(*[ubx]*N), cs.vertcat(*[ubz]*N))

# Symbolic variables
X = cs.SX.sym('x', nx, N)
X_vec = cs.reshape(X, nx*N, 1)

Z = cs.SX.sym('z', nz, N)
Z_vec = cs.reshape(Z, nz*N, 1)

m = cs.vertcat(X_vec, Z_vec)

# build problem
H = []; lbh = []; ubh = []  # Equality constraints
G = []; lbg = []; ubg = []  # Inequality constraints


for i in range(N):
    if i == 0:
        x_i = x0
    else:
        x_i = X[0:nx, i-1]
    z_i = Z[0:nz, i]
    #u_i = U[0:nu, i]
    u_i = u0

    # rk4 on x_i
    k1 = h * LotkaVolterraModel(x=x_i, u=u_i, z=z_i)
    k2 = h * LotkaVolterraModel(x_i + k1 / 2, u=u_i, z=z_i)
    k3 = h * LotkaVolterraModel(x_i + k2 / 2, u=u_i, z=z_i)
    k4 = h * LotkaVolterraModel(x_i + k3, u=u_i, z=z_i)
    dx = (k1 + 2 * k2 + 2 * k3 + k4) / 6

    # Equality constraints
    H.append(x_i + dx - X[0:nx, i])  # For the multiple shooting
    lbh.append(-tol_cons_eq*cs.DM.ones(nx, 1))
    ubh.append(tol_cons_eq*cs.DM.ones(nx, 1))

    # algebraic
    H.append(equality_constraints(x=x_i, u=u_i, z=z_i))  # add the equality constraints over z --> on veut égale à zéro
    lbh.append(-tol_cons_eq*cs.DM.ones(nz, 1))
    ubh.append(tol_cons_eq*cs.DM.ones(nz, 1))



# Concatenate the constraints
G = cs.vertcat(*G)
H = cs.vertcat(*H)

cons = cs.Function('cons', [m], [cs.vertcat(G, H)])
lb_cons = cs.vertcat(cs.vertcat(*lbg), cs.vertcat(*lbh))
ub_cons = cs.vertcat(cs.vertcat(*ubg), cs.vertcat(*ubh))

# f = cs.Function('f', [m], [h**2 + cs.sumsqr(S)])
f = cs.Function('f', [m], [h**2])

# Defining the problem
m_s = cs.MX.sym('m_s', N*(nx+nz), 1)  # scaled variables
m_unsc = unormalize(ms=m_s, lb=lbm, ub=ubm)  # unscaled variables
m0_s = normalize(cs.vertcat(cs.vertcat(*[x0]*N), cs.vertcat(*[z0]*N)), lb=lbm, ub=ubm)

nlp = {'x': m_s, 'f': f(m_unsc), 'g': cons(m_unsc)}
opts = {'ipopt': {'max_iter': 10000, 'acceptable_tol': tol_cons_eq}}
S = cs.nlpsol('S', 'ipopt', nlp, opts,)


print(S)
t_start = time.perf_counter()
res = S(x0=m0_s,
        lbx=cs.DM.zeros((nz+nx)*N, 1),
        ubx=cs.DM.ones((nz+nx)*N, 1),
        lbg=lb_cons,
        ubg=ub_cons
        )

t_end = time.perf_counter()
print(t_end - t_start)

# process
m_opt_s = res['x']
m_opt_unsc = unormalize(ms=m_opt_s, ub=ubm, lb=lbm) # Need to unscale the optimal value

sol = {'x': cs.reshape(m_opt_unsc[0:nx*N], nx, N).toarray(),
       'z': cs.reshape(m_opt_unsc[nx*N:(nx+nz)*N], nz, N).toarray(),
       'u': u0,
       'h': h,
       'obj': res['f'].toarray(),
       'num_iteration': N,
       'm_opt_s': m_opt_s
       }

time.sleep(1)

# save
save_sol('sol_tmp', sol)


# call plot
solution = load_sol('sol_tmp')
plot_solution(solution)



print('End')

