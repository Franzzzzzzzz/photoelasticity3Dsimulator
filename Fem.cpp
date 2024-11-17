#include "Fem.h"

void FEsolver::prepare_mesh (std::string filename, double scaling) 
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
  
  std::vector<double> dst ; dst.resize(6) ; 
  double dstmax=0 ; 
  for (int i=0 ; i<ncell ; i++)
  {
    auto c = [&](int pt, int dim){return meshpoints[topology[i*CONN+pt]*DIM+dim] ; };
    for (int j=0 ; j<CONN ; j++)
      for (int dd=0 ; dd<DIM ; dd++)
        com[i*DIM+dd] += meshpoints[topology[i*CONN+j]*DIM+dd]/4. ; 
    
    dst[0] = (c(0,0)-c(1,0))*(c(0,0)-c(1,0))+(c(0,1)-c(1,1))*(c(0,1)-c(1,1))+(c(0,2)-c(1,2))*(c(0,2)-c(1,2)) ; 
    dst[1] = (c(0,0)-c(2,0))*(c(0,0)-c(2,0))+(c(0,1)-c(2,1))*(c(0,1)-c(2,1))+(c(0,2)-c(2,2))*(c(0,2)-c(2,2)) ; 
    dst[2] = (c(0,0)-c(3,0))*(c(0,0)-c(3,0))+(c(0,1)-c(3,1))*(c(0,1)-c(3,1))+(c(0,2)-c(3,2))*(c(0,2)-c(3,2)) ; 
    dst[3] = (c(1,0)-c(2,0))*(c(1,0)-c(2,0))+(c(1,1)-c(2,1))*(c(1,1)-c(2,1))+(c(1,2)-c(2,2))*(c(1,2)-c(2,2)) ; 
    dst[4] = (c(1,0)-c(3,0))*(c(1,0)-c(3,0))+(c(1,1)-c(3,1))*(c(1,1)-c(3,1))+(c(1,2)-c(3,2))*(c(1,2)-c(3,2)) ; 
    dst[5] = (c(2,0)-c(3,0))*(c(2,0)-c(3,0))+(c(2,1)-c(3,1))*(c(2,1)-c(3,1))+(c(2,2)-c(3,2))*(c(2,2)-c(3,2)) ; 
    if (*std::max_element(dst.begin(), dst.end())>dstmax)
      dstmax = *std::max_element(dst.begin(), dst.end()); 
  }
  std::vector<std::pair<double,double>> boundaries = {{-1,1},{-1,1},{-1,1}} ;
  cells.init_cells(boundaries, 1.1*sqrt(dstmax)) ; 
  cells.allocate_to_cells(com) ; 
  
  // getting connectivity and tetrahedra
  filename = filename+".con" ; 
  FILE *in=fopen(filename.c_str(), "rb") ;
  if (in)
  {
    printf("%s exists, using it for the mesh connectivity.\n", filename.c_str()) ; 
    // Reading the connectivity file
    size_t tmp; size_t res ; 
    res=fread(&tmp, sizeof(size_t), 1, in) ; 
    if (res!=1) printf("Error reading size") ; 
    tetras.resize(tmp) ; 
    for (size_t i=0 ; i<tetras.size() ; i++)
      tetras[i].read(in) ;
    fclose(in); 
  }
  else
  {
    printf("Creating the connectivity file. This may take time, but if you do not change the mesh the .con file created may be reused.\n") ; 
    tetras.resize(ncell) ; 
    auto coord = [&](int cell,int pt,int dim){return meshpoints[(topology[cell*CONN+pt]*DIM+dim)];} ; 
    for (int i=0 ; i<ncell ; i++)
    {
      tetras[i].vertex_id={topology[i*CONN+0], topology[i*CONN+1], topology[i*CONN+2], topology[i*CONN+3]} ; 
      
      tetras[i].triangles[0].M = {coord(i,0,0), coord(i,0,1), coord(i,0,2)} ; 
      tetras[i].triangles[1].M = tetras[i].triangles[0].M ;
      tetras[i].triangles[2].M = tetras[i].triangles[0].M ;
      tetras[i].com = tetras[i].triangles[0].M ; 
      Eigen::Vector3d pt = {coord(i,1,0), coord(i,1,1), coord(i,1,2)} ;
      tetras[i].triangles[0].u = pt - tetras[i].triangles[0].M ; 
      tetras[i].triangles[1].u = pt - tetras[i].triangles[1].M ; 
      tetras[i].triangles[3].M = pt ;
      tetras[i].com += pt ; 
      pt = Eigen::Vector3d(coord(i,2,0), coord(i,2,1), coord(i,2,2)) ;
      tetras[i].triangles[0].v = pt - tetras[i].triangles[0].M ; 
      tetras[i].triangles[2].u = pt - tetras[i].triangles[2].M ; 
      tetras[i].triangles[3].u = pt - tetras[i].triangles[3].M ;
      tetras[i].com += pt ; 
      pt = Eigen::Vector3d(coord(i,3,0), coord(i,3,1), coord(i,3,2)) ;
      tetras[i].triangles[1].v = pt - tetras[i].triangles[1].M ; 
      tetras[i].triangles[2].v = pt - tetras[i].triangles[2].M ; 
      tetras[i].triangles[3].v = pt - tetras[i].triangles[3].M ;
      tetras[i].com += pt ; 
      tetras[i].com /= 4 ;
    }
    
    // Finding connectivity
    std::vector<int> intersection ;
    for (int i=0 ; i<ncell ; i++)
    {
      std::set<int> tetra1{topology[i*CONN+0], topology[i*CONN+1], topology[i*CONN+2], topology[i*CONN+3]} ; 
      for (size_t j=i+1 ; j<topology.size()/CONN ; j++)
      {
        std::set<int> tetra2{topology[j*CONN+0], topology[j*CONN+1], topology[j*CONN+2], topology[j*CONN+3]} ; 
        intersection.clear() ; 
        std::set_intersection(tetra1.begin(), tetra1.end(), tetra2.begin(), tetra2.end(), std::back_inserter(intersection)) ;
        if (intersection.size()==3)
        {
          char tri_i= tetras[i].identify_tri(intersection) ; 
          char tri_j= tetras[j].identify_tri(intersection) ; 
          
          tetras[i].neigh.insert(j) ; 
          tetras[i].surface |= 1<<tri_i ; 
          tetras[j].neigh.insert(i) ; 
          tetras[j].surface |= 1<<tri_j ; 
        }
        if (intersection.size()>=1)
        {
          tetras[i].neigh_extended.insert(j) ; 
          tetras[j].neigh_extended.insert(i) ; 
        }
      }
    }
    for (int i=0 ; i<ncell ; i++)
      tetras[i].surface = (~tetras[i].surface)&(0x0F) ; 
    
    // Writing the connectivity file
    FILE *out=fopen(filename.c_str(), "wb") ;
    size_t tmp = tetras.size() ; 
    fwrite(&tmp, sizeof(size_t), 1, out) ; 
    for (size_t i=0 ; i<tetras.size() ; i++)
      tetras[i].write(out) ;
    fclose(out);      
  }

  printf("Finished preparing the mesh\n") ; fflush(stdout) ; 
}
//------------------------------------------------
void FEsolver::get_sigma(std::vector<double> &result)
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
std::vector<int> FEsolver::interpolate(std::vector<vec> & pos) 
{
  std::vector<int> res ; 
  res.resize(pos.size()) ; 
  for (size_t i=0 ; i<pos.size() ; i++)
    res[i] = interpolate_nearestneighbour_cell(pos[i]) ; 
  return res ; 
}
//------------------------------------------------
int FEsolver::interpolate_nearestneighbour_cell (vec & pos)
{
  return cells.closest(com,pos) ; 
}
//------------------------------------------------
int FEsolver::interpolate_nearestneighbour(vec & pos)
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
