auto IDENTIFIER( bdofs , THISID ) = fem::locate_dofs_geometrical(*V, [=](auto x) {
    std::vector<std::int8_t> marker(x.extent(1), false);
    double dst = 1e8 ; size_t prev =0 ; 
    
    for (std::size_t p = 0; p < x.extent(1); ++p)
    {
      double d = (x(0,p)-c0)*(x(0,p)-c0) + (x(1,p)-c1)*(x(1,p)-c1) + (x(2,p)-c2)*(x(2,p)-c2) ; 
      if (d<dst)
      {
        dst = d ; 
        marker[prev] = false ; 
        prev = p ;
        marker[prev] = true ; 
      }
    }
    return marker ;     
  }) ; 
auto IDENTIFIER( bc , THISID) = std::make_shared<const fem::DirichletBC<T>>(std::vector<T>{d0, d1, d2}, IDENTIFIER( bdofs, THISID ), V);
