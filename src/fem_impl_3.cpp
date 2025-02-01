auto A = la::petsc::Matrix(fem::petsc::create_matrix(*a), false);
la::Vector<T> b(L->function_spaces()[0]->dofmap()->index_map,
                L->function_spaces()[0]->dofmap()->index_map_bs());


MatZeroEntries(A.mat());
fem::assemble_matrix(la::petsc::Matrix::set_block_fn(A.mat(), ADD_VALUES),
                      *a, THISBCS() );
MatAssemblyBegin(A.mat(), MAT_FLUSH_ASSEMBLY);
MatAssemblyEnd(A.mat(), MAT_FLUSH_ASSEMBLY);
fem::set_diagonal<T>(la::petsc::Matrix::set_fn(A.mat(), INSERT_VALUES), *V,
                      THISBCS() );
MatAssemblyBegin(A.mat(), MAT_FINAL_ASSEMBLY);
MatAssemblyEnd(A.mat(), MAT_FINAL_ASSEMBLY);

b.set(0.0);
fem::assemble_vector(b.mutable_array(), *L);
fem::apply_lifting<T, U>(b.mutable_array(), {a}, { THISBCS() }, {}, T(1));
b.scatter_rev(std::plus<T>());
//bc->set(b.mutable_array(), {});

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
    = mesh::cell_type_to_basix_type(mesh_scaled[selected_mesh].mesh->topology()->cell_type());
constexpr bool discontinuous = true;
basix::FiniteElement S_element = basix::create_element<U>(
    family, cell_type, 0, basix::element::lagrange_variant::unset,
    basix::element::dpc_variant::unset, discontinuous);
auto S = std::make_shared<fem::FunctionSpace<U>>(fem::create_functionspace(
    mesh_scaled[selected_mesh].mesh, S_element, std::vector<std::size_t>{3, 3}));
auto sigma_expression = fem::create_expression<T, U>(
    *expression_sphere_stress, {{"uh", u}}, {{"mu", mu}, {"lambda_", lbd}});

auto sigma = fem::Function<double>(S);
sigma.name = "cauchy_stress";
#if DOLFINX_VERSION_MINOR <= 8 
sigma.interpolate(sigma_expression, *(mesh_scaled[selected_mesh].mesh));
#else
sigma.interpolate(sigma_expression);
#endif
auto res = sigma.x()->mutable_array() ; 

result.resize(res.size()) ; 
std::copy(res.begin(), res.end(), result.begin()) ; 

// Save solution in VTK format
io::VTKFile file(MPI_COMM_WORLD, "u.pvd", "w");
file.write<T>({*u}, 0.0);
// Save Cauchy stress in XDMF format
io::XDMFFile file_sigma(mesh_scaled[selected_mesh].mesh->comm(), "sigma.xdmf", "w");
file_sigma.write_mesh(*(mesh_scaled[selected_mesh].mesh));
file_sigma.write_function(sigma, 0.0);
