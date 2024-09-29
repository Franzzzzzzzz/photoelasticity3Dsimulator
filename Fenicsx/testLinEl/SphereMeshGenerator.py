from mpi4py import MPI
from dolfinx.io import XDMFFile, gmshio
import gmsh  # type: ignore

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
  
def create_mesh(comm: MPI.Comm, model: gmsh.model, name: str, filename: str, mode: str):
    """Create a DOLFINx from a Gmsh model and output to file.

    Args:
        comm: MPI communicator top create the mesh on.
        model: Gmsh model.
        name: Name (identifier) of the mesh to add.
        filename: XDMF filename.
        mode: XDMF file mode. "w" (write) or "a" (append).
    """
    msh, ct, ft = gmshio.model_to_mesh(model, comm, rank=0)
    msh.name = name
    ct.name = f"{msh.name}_cells"
    ft.name = f"{msh.name}_facets"
    with XDMFFile(msh.comm, filename, mode) as file:
        msh.topology.create_connectivity(2, 3)
        file.write_mesh(msh)
        file.write_meshtags(
            ct, msh.geometry, geometry_xpath=f"/Xdmf/Domain/Grid[@Name='{msh.name}']/Geometry"
        )
        file.write_meshtags(
            ft, msh.geometry, geometry_xpath=f"/Xdmf/Domain/Grid[@Name='{msh.name}']/Geometry"
        )

  
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)

# Create model
model = gmsh.model()

model = gmsh_sphere(model, "Sphere")
#model.mesh.refine()
#model.mesh.refine()
model.setCurrent("Sphere")

create_mesh(MPI.COMM_WORLD, model, "sphere", "sphere_mesh.xdmf", "a")
