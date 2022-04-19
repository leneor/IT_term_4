from __future__ import print_function
from fenics import *
from ufl import nabla_div
import matplotlib.pyplot as plt
from vtkplotter.dolfin import plot

# Параметры балки
l = 1 # Длина балки, метры, Ox
H = 0.2 # Высота балки, метры, Oy
W = 0.1 # Толщина балки, метры, Oz

# Параметры материала, из которого сделана балка
E = 2 * 10 ** 11  # Модуль Юнга, Па
rho = 7900  # Плотность материала балки
nu = 0.29  # Коэффициент Пуассона
mu = E / (2 + 2 * nu)  # Модуль сдвига
lambda_ = E * nu / ((1 + nu) * (1 - 2 * nu))  # Коэффициент Ламе
g = 9.80665  # Ускорение свободного падения

F = 500 #Сила на правом конце, Н
Q = F / (W*H*g) #Напряжение на правом конце, кгс

# Создаем модель балки
mesh = BoxMesh(Point(0, 0, 0), Point(l, H, W), 10000, 20, 1)
V = VectorFunctionSpace(mesh, 'P', 1)
# Определяем начальные условия
tol = 1E-14


def clamped_boundary(x, on_boundary):
    return on_boundary and x[0] < tol


bc = DirichletBC(V, Constant((0, 0, 0)), clamped_boundary)


# Определение деформации и напряжения
def epsilon(u):
    return 0.5 * (nabla_grad(u) + nabla_grad(u).T)


def sigma(u):
    return lambda_ * nabla_div(u) * Identity(d) + 2 * mu * epsilon(u)


# Определяем условия задачи
u = TrialFunction(V)
d = u.geometric_dimension()
v = TestFunction(V)
f = Constant((0, -rho*g, 0))  # Сила тяжести на единицу объема
T = Constant((0, -Q, 0))  # Напряжение на конце, Н/м^2
a = inner(sigma(u), epsilon(v)) * dx
L = dot(f, v) * dx + dot(T, v) * ds

# Вычисление решения
u = Function(V)
solve(a == L, u, bc)
V = FunctionSpace(mesh, 'P', 1)

# Вычисление величины смещения
u_magnitude = sqrt(dot(u, u))
u_magnitude = project(u_magnitude, V)
print('Смещение на правом конце балки:',
      round(u_magnitude.vector().get_local().max() * 1000, 3), 'mm.')

# Сохранить решение в VTK формат
File('results/magnitude.pvd') << u_magnitude
