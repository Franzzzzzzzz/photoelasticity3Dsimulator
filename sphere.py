import pyvista
import dolfinx
from dolfinx import mesh, fem, plot, io
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI
import ufl
import numpy as np
import gmsh

L = 1
W = 0.2
mu = 1
rho = 1
delta = W / L
gamma = 0.4 * delta**2
beta = 1.25
lambda_ = beta
g = gamma

#domain = mesh.create_box(MPI.COMM_WORLD, [np.array([0, 0, 0]), np.array([L, W, W])],
#                         [20, 6, 6], cell_type=mesh.CellType.hexahedron)

def gmsh_sphere(model: gmsh.model, name: str) -> gmsh.model:
    """Create a Gmsh model of a sphere.

    Args:
        model: Gmsh model to add the mesh to.
        name: Name (identifier) of the mesh to add.

    Returns:
        Gmsh model with a sphere mesh added.

    """
    model.add(name)
    model.setCurrent(name)
    sphere = model.occ.addSphere(0, 0, 0, 1, tag=1)

    # Synchronize OpenCascade representation with gmsh model
    model.occ.synchronize()

    # Add physical marker for cells. It is important to call this
    # function after OpenCascade synchronization
    model.add_physical_group(dim=3, tags=[sphere])

    # Generate the mesh
    model.mesh.generate(dim=3)
    return model
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)

# Create model
model = gmsh.model()

model = gmsh_sphere(model, "Sphere")
#model.mesh.refine()
#model.mesh.refine()
model.setCurrent("Sphere")

import basix.ufl
domain, cell_tags, facet_tags = io.gmshio.model_to_mesh(model, MPI.COMM_WORLD, 0)
element = basix.ufl.element("Lagrange", domain.topology.cell_name(), 1, shape=(domain.geometry.dim, ))

#V = fem.functionspace(domain, element)
V = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim, )))

def clamped_boundary(x):
    return x[0]<-0.9 #np.sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])<1
def clamped_boundary(x):
    return np.isclose(x[0], 0)


fdim = domain.topology.dim - 1
boundary_facets = mesh.locate_entities_boundary(domain, fdim, clamped_boundary)
print(boundary_facets)
u_D = np.array([0., 0., 0.])
bc = fem.dirichletbc(u_D, fem.locate_dofs_topological(V, fdim, boundary_facets), V)

T = fem.Constant(domain, (0., 0., 0.))
ds = ufl.Measure("ds", domain=domain)

def epsilon(u):
    return ufl.sym(ufl.grad(u))  # Equivalent to 0.5*(ufl.nabla_grad(u) + ufl.nabla_grad(u).T)


def sigma(u):
    return lambda_ * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2 * mu * epsilon(u)


u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
f = fem.Constant(domain, (0., 0., rho*g))
a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx
L = ufl.dot(f, v) * ufl.dx + ufl.dot(T, v) * ds


problem = LinearProblem(a, L, bcs=[bc], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
uh = problem.solve()

stress =  lambda_ * ufl.nabla_div(uh) * ufl.Identity(len(uh)) + 2 * mu * epsilon(uh)
expressions = [(stress, [[0.25, 0.25, 0.25]])]

with io.XDMFFile(domain.comm, "deformation.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    uh.name = "Deformation"
    xdmf.write_function(uh)
