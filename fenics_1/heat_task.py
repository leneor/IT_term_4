from fenics import *
import math

T = 90.0            # final time
num_steps = 100     # number of time steps
dt = T / num_steps  # time step size
alpha = -1          # parameter alpha
beta = 1.5          # parameter beta

# Create mesh and define function space
nx = 20
ny = 20
mesh = UnitSquareMesh(nx, ny)
V = FunctionSpace(mesh, 'Lagrange', 2)

# Define boundary condition
u_D = Expression('3 - beta*(x[0]*x[0]-91*x[0]) - alpha*(x[1]*x[1]-100*x[1]) + beta*t', degree=1, alpha=alpha, beta=beta, t=0)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)

# Define initial value
u_n = interpolate(u_D, V)
#u_n = project(u_D, V)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(beta - 2 - 2*alpha*alpha)

F = u*v*dx + dt*dot(grad(u), grad(v))*dx - (u_n + dt*f)*v*dx
a, L = lhs(F), rhs(F)

res_file = File('heat_task/solution.pvd')

# Time-stepping
u = Function(V)
t = 0
for n in range(num_steps):

    # Update current time
    t += dt
    u_D.t = t

    # Compute solution
    solve(a == L, u, bc)

    # Save solution to VTK
    res_file << u

    # Update previous solution
    u_n.assign(u)
