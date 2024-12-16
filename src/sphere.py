import pyvista
import dolfinx
from dolfinx import mesh, fem, plot, io
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI
import ufl
import numpy as np
import gmsh

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
mu = fem.Constant(domain, 1.)
lambda_ = fem.Constant(domain, 1.25)
#mu = 1
#lambda_ = 1.25
element = basix.ufl.element("Lagrange", domain.topology.cell_name(), 1, shape=(domain.geometry.dim, ))

#V = fem.functionspace(domain, element)
V = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim, )))

def clamped_boundary(x):
    return x[2]<-0.9 #np.sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])<1
#def clamped_boundary(x):
    #return np.isclose(x[0], 0)


fdim = domain.topology.dim - 1
boundary_facets = mesh.locate_entities_boundary(domain, fdim, clamped_boundary)
print(boundary_facets)
u_D = np.array([0., 0., 0.])
bc = fem.dirichletbc(u_D, fem.locate_dofs_topological(V, fdim, boundary_facets), V)

T = fem.Constant(domain, (0., 0., 0.))
#T = fem.Function(domain)
#T.x.array[:]=0  ;
#T.x.array[0,:]=(10000, 0, 0) ;
ds = ufl.Measure("ds", domain=domain)

def epsilon(u):
    return ufl.sym(ufl.grad(u))  # Equivalent to 0.5*(ufl.nabla_grad(u) + ufl.nabla_grad(u).T)


def sigma(u):
    return lambda_ * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2 * mu * epsilon(u)

#class PointLoad(UserExpression):
    #def __init__(self, pt, vl,tol,**kwargs):
        #super().__init__(**kwargs)
        #self.point = pt
        #self.value = vl
        #self.tol = tol
    #def eval(self, values, x):
        #if near (x[0], self.point[0],self.tol) and near(x[1], self.point[1],self.tol) and near(x[2], self.point[2],self.tol):
            #values[0] = self.value[0]
            #values[1] = self.value[1]
            #values[2] = self.value[2]
        #else:
            #values[0] = 0
            #values[1] = 0
            #values[2] = 0
    #def value_shape(self):
        #return (2,)

u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
f = fem.Constant(domain, (0., 0., 0.))
#f = PointLoad(pt=(2.5,0), vl=(a0,b0), tol=1e-1,degree = 1)
#f = fem.Function(V)
#f.x.array[:]=0  ;
#f.x.array[0::3]=rho*g ;

a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx
L = ufl.dot(f, v) * ufl.dx + ufl.dot(T, v) * ds

#problem = LinearProblem(a, L, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
problem = LinearProblem(a, L, bcs=[bc], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
uh = problem.solve()

stress =  lambda_ * ufl.nabla_div(uh) * ufl.Identity(len(uh)) + 2 * mu * epsilon(uh)
expressions = [(stress, [[0.25, 0.25, 0.25]])]

with io.XDMFFile(domain.comm, "deformation.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    uh.name = "Deformation"
    xdmf.write_function(uh)
