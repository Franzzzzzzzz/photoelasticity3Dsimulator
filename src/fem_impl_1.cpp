auto element = basix::create_element<U>(
      basix::element::family::P, basix::cell::type::tetrahedron, 1,
      basix::element::lagrange_variant::unset,
      basix::element::dpc_variant::unset, false);
auto V = std::make_shared<fem::FunctionSpace<U>>(fem::create_functionspace(mesh_scaled[selected_mesh].mesh, element, {3}));

auto f = std::make_shared<fem::Constant<T>>(std::vector<T>{0, 0, 0});
auto g = std::make_shared<fem::Constant<T>>(std::vector<T>{0, 0, 0});
auto mu = std::make_shared<fem::Constant<T>>(1);
auto lbd = std::make_shared<fem::Constant<T>>(1.25);

// Define variational forms
auto a = std::make_shared<fem::Form<T>>(fem::create_form<T>(*form_sphere_a, {V, V}, {}, {{"mu", mu}, {"lambda_", lbd}}, {}, {}));
auto L = std::make_shared<fem::Form<T>>(fem::create_form<T>(*form_sphere_L, {V}, {}, {{"f", f}, {"T", g}}, {}, {}));
    
auto u = std::make_shared<fem::Function<T>>(V);
