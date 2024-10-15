#include "sphere.h"
#include <basix/finite-element.h>
#include <cmath>
#include <dolfinx.h>
#include <dolfinx/fem/Constant.h>
#include <dolfinx/fem/petsc.h>
#include <dolfinx/la/petsc.h>
#include <dolfinx/io/XDMFFile.h>
#include <petscmat.h>
#include <petscsys.h>
#include <petscsystypes.h>
#include <utility>
#include <vector>

using namespace dolfinx;
using T = PetscScalar;
using U = typename dolfinx::scalar_value_type_t<T>;

class FEsolver{
public:
  FEsolver(int argc, char **argv[]) 
  {
    dolfinx::init_logging(argc, *argv);
    PetscInitialize(&argc, argv, nullptr, nullptr);
  }
  //------------------------------------------------  
  void prepare_mesh (std::string filename, double scaling=1) 
  {
    io::XDMFFile meshfile(MPI_COMM_WORLD, filename.c_str(), "r");
    auto cel = fem::CoordinateElement<double>(mesh::CellType::tetrahedron, 1, basix::element::lagrange_variant::unset);
    mesh = std::make_shared<mesh::Mesh<double>>(meshfile.read_mesh(cel, mesh::GhostMode::shared_facet, "sphere"));
    
    
    // Calculate cells center of mass
    auto meshpoints=mesh->geometry().x() ;
    if (scaling != 1)
      for (size_t i=0 ; i<meshpoints.size() ; i++)
        meshpoints[i] *= scaling ;
    auto topology=mesh->topology()->connectivity(3,0)->array() ; 
    
    const int CONN = 4 ; // TODO do not assume tetrahedron
    const int DIM = 3 ; 
    ncell=topology.size()/CONN ; 
    com.resize(ncell*DIM,0) ; 
    for (int i=0 ; i<ncell ; i++)
      for (int j=0 ; j<CONN ; j++)
        for (int dd=0 ; dd<DIM ; dd++)
          com[i*DIM+dd] += meshpoints[topology[i*CONN+j]*DIM+dd]/4. ;     
  }
  //------------------------------------------------
  void get_sigma(std::vector<double> &result)
  {    
    auto element = basix::create_element<U>(
        basix::element::family::P, basix::cell::type::tetrahedron, 1,
        basix::element::lagrange_variant::unset,
        basix::element::dpc_variant::unset, false);
    auto V = std::make_shared<fem::FunctionSpace<U>>(fem::create_functionspace(mesh, element, {3}));
    
    auto f = std::make_shared<fem::Constant<T>>(std::vector<T>{0, 0, 0.016});
    auto g = std::make_shared<fem::Constant<T>>(std::vector<T>{0, 0, 0});
    
    // Define variational forms
    auto a = std::make_shared<fem::Form<T>>(fem::create_form<T>(*form_sphere_a, {V, V}, {}, {}, {}, {}));
    auto L = std::make_shared<fem::Form<T>>(fem::create_form<T>(*form_sphere_L, {V}, {}, {{"f", f}, {"T", g}}, {}, {}));
        
    auto u = std::make_shared<fem::Function<T>>(V);
  
    auto bdofs_left = fem::locate_dofs_geometrical(
        *V,
        [](auto x)
        {
          //constexpr U eps = 1.0e-6;
          std::vector<std::int8_t> marker(x.extent(1), false);
          
          /*for (std::size_t p = 0; p < x.extent(1); ++p)
          {
            if (std::abs(x(0, p)) < eps)
            {
              marker[p] = true;
            }
          }*/
          for (std::size_t p = 0; p < x.extent(1); ++p)
          {
            if (x(1, p)<-0.45)
            {
              marker[p] = true;
            }
          }          
          return marker;
        });
    
    auto bc = std::make_shared<const fem::DirichletBC<T>>(std::vector<T>{0, 0, 0},bdofs_left, V);
    
    auto A = la::petsc::Matrix(fem::petsc::create_matrix(*a), false);
    la::Vector<T> b(L->function_spaces()[0]->dofmap()->index_map,
                    L->function_spaces()[0]->dofmap()->index_map_bs());

    
    MatZeroEntries(A.mat());
    fem::assemble_matrix(la::petsc::Matrix::set_block_fn(A.mat(), ADD_VALUES),
                         *a, {bc});
    MatAssemblyBegin(A.mat(), MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(A.mat(), MAT_FLUSH_ASSEMBLY);
    fem::set_diagonal<T>(la::petsc::Matrix::set_fn(A.mat(), INSERT_VALUES), *V,
                         {bc});
    MatAssemblyBegin(A.mat(), MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A.mat(), MAT_FINAL_ASSEMBLY);

    b.set(0.0);
    fem::assemble_vector(b.mutable_array(), *L);
    fem::apply_lifting<T, U>(b.mutable_array(), {a}, {{bc}}, {}, T(1));
    b.scatter_rev(std::plus<T>());
    bc->set(b.mutable_array(), {});

    la::petsc::KrylovSolver lu(MPI_COMM_WORLD);
    la::petsc::options::set("ksp_type", "preonly");
    la::petsc::options::set("pc_type", "lu");
    lu.set_from_options();

    lu.set_operator(A.mat());
    la::petsc::Vector _u(la::petsc::create_vector_wrap(*u->x()), false);
    la::petsc::Vector _b(la::petsc::create_vector_wrap(b), false);
    lu.solve(_u.vec(), _b.vec());

    // Update ghost values before output
    u->x()->scatter_fwd();

    //  The function `u` will be modified during the call to solve. A
    //  {cpp:class}`Function` can be saved to a file. Here, we output
    //  the solution to a `VTK` file (specified using the suffix `.pvd`)
    //  for visualisation in an external program such as Paraview.
    constexpr auto family = basix::element::family::P;
        
    auto cell_type
        = mesh::cell_type_to_basix_type(mesh->topology()->cell_type());
    constexpr bool discontinuous = true;
    basix::FiniteElement S_element = basix::create_element<U>(
        family, cell_type, 0, basix::element::lagrange_variant::unset,
        basix::element::dpc_variant::unset, discontinuous);
    auto S = std::make_shared<fem::FunctionSpace<U>>(fem::create_functionspace(
        mesh, S_element, std::vector<std::size_t>{3, 3}));
    auto sigma_expression = fem::create_expression<T, U>(
        *expression_sphere_stress, {{"uh", u}}, {});

    auto sigma = fem::Function<double>(S);
    sigma.name = "cauchy_stress";
    sigma.interpolate(sigma_expression, *mesh);
    auto res = sigma.x()->mutable_array() ; 
    
    result.resize(res.size()) ; 
    std::copy(res.begin(), res.end(), result.begin()) ; 
    
    // Save solution in VTK format
    io::VTKFile file(MPI_COMM_WORLD, "u.pvd", "w");
    file.write<T>({*u}, 0.0);
    // Save Cauchy stress in XDMF format
    io::XDMFFile file_sigma(mesh->comm(), "sigma.xdmf", "w");
    file_sigma.write_mesh(*mesh);
    file_sigma.write_function(sigma, 0.0);
  }
  //------------------------------------------------
  
  std::vector<int> interpolate(std::vector<vec> & pos) 
  {
    std::vector<int> res ; 
    res.resize(pos.size()) ; 
    for (size_t i=0 ; i<pos.size() ; i++)
      res[i] = interpolate_nearestneighbour(pos[i]) ; 
    return res ; 
  }
  int interpolate_nearestneighbour(vec & pos)
  {
    double mindst ; 
    int idmin = 0 ; 
    double dst = 0 ; 
    for (int i=0 ; i<ncell ; i++)
    {
      dst = (pos[0]-com[i*3+0])*(pos[0]-com[i*3+0])+
            (pos[1]-com[i*3+1])*(pos[1]-com[i*3+1])+
            (pos[2]-com[i*3+2])*(pos[2]-com[i*3+2]) ; 
      if (i==0) mindst = dst ; 
      if (dst<mindst)
      {
        mindst=dst ; 
        idmin = i ; 
      }
    }
    return idmin ;     
  }
  
  std::shared_ptr<mesh::Mesh<double>> mesh ;
  int ncell ; 
  std::vector<double> com ; 
} ; 





