import random
from dolfin import *
import numpy as np
import os


class InitialConditions(UserExpression):
    def __init__(self, **kwargs):
        random.seed(2 + MPI.rank(MPI.comm_world))
        super().__init__(**kwargs)

    def eval(self, values, x):
        values[0] = 0.63 + 0.02 * (0.5 - random.random())
        values[1] = 0.0

    def value_shape(self):
        return (2,)


class CahnHilliardEquation(NonlinearProblem):
    def __init__(self, a, L):
        NonlinearProblem.__init__(self)
        self.L = L
        self.a = a

    def F(self, b, x):
        assemble(self.L, tensor=b)

    def J(self, A, x):
        assemble(self.a, tensor=A)


def create_d(lmbda):
    os.makedirs("lmbda" + str(lmbda))
    dt = 1.0e-05
    theta = 0.5

    parameters["form_compiler"]["optimize"] = True
    parameters["form_compiler"]["cpp_optimize"] = True

    mesh = UnitSquareMesh.create(96, 96, CellType.Type.quadrilateral)
    P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
    ME = FunctionSpace(mesh, P1 * P1)

    du = TrialFunction(ME)
    q, v = TestFunctions(ME)

    u = Function(ME)
    u0 = Function(ME)

    dc, dmu = split(du)
    c, mu = split(u)
    c0, mu0 = split(u0)

    u_init = InitialConditions(degree=1)
    u.interpolate(u_init)
    u0.interpolate(u_init)

    c = variable(c)
    f = 100 * c ** 2 * (1 - c) ** 2
    dfdc = diff(f, c)

    mu_mid = (1.0 - theta) * mu0 + theta * mu

    L0 = c * q * dx - c0 * q * dx + dt * dot(grad(mu_mid), grad(q)) * dx
    L1 = mu * v * dx - dfdc * v * dx - lmbda * dot(grad(c), grad(v)) * dx
    L = L0 + L1

    a = derivative(L, u, du)

    problem = CahnHilliardEquation(a, L)
    solver = NewtonSolver()
    solver.parameters["linear_solver"] = "lu"
    solver.parameters["convergence_criterion"] = "incremental"
    solver.parameters["relative_tolerance"] = 1e-6

    file = File("output.pvd", "compressed")

    ind = 0
    t = 0.0
    T = 50 * dt
    while t < T:
        t += dt
        u0.vector()[:] = u.vector()
        solver.solve(problem, u.vector())
        np.save(
            "lmbda" + str(lmbda) + "/numpy_t" + str(ind + 1),
            u.split()[0].vector().get_local(),
        )

        ind += 1
        element = ME.element()
        dofmap = ME.dofmap()
        ii = 0
        for cell in cells(mesh):
            if ii == 0:
                all_coords = element.tabulate_dof_coordinates(cell)
                all_order = dofmap.cell_dofs(cell.index())
            else:
                all_coords = np.concatenate(
                    (all_coords, element.tabulate_dof_coordinates(cell))
                )
                all_order = np.concatenate((all_order, dofmap.cell_dofs(cell.index())))
            ii += 1
        np.save("lmbda" + str(lmbda) + "/numpy_dims" + str(ind), all_coords)
        np.save("lmbda" + str(lmbda) + "/numpy_order" + str(ind), all_order)


VALS = np.linspace(1.0e-02, 0.0575, 20)

for i in range(len(VALS)):
    lmbda_val = round(VALS[i], 4)
    create_d(lmbda_val)
