#include "Geometry.h"


std::optional<entryexit> Geometry::intersection_ray_sphere (vec x0, vec xr, vec xc, double r) 
{
  double a = xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2] ; 
  double b = 2*((x0[0]-xc[0])*xr[0] + (x0[1]-xc[1])*xr[1] + (x0[2]-xc[2])*xr[2]) ;
  double c = (x0[0]-xc[0])*(x0[0]-xc[0]) + (x0[1]-xc[1])*(x0[1]-xc[1]) + (x0[2]-xc[2])*(x0[2]-xc[2]) - r*r;
  
  double delta = b*b - 4 * a * c ; 
  if (delta<=0) return {} ; 
  entryexit res ; 
  
  double alpha1 = (-b - sqrt(delta)) / (2*a) ; 
  double alpha2 = (-b + sqrt(delta)) / (2*a) ; 
  
  res.pt1 = x0 + alpha1*xr ; 
  res.pt2 = x0 + alpha2*xr ;

  return res ; 
}
//-----------------------------------------------------------------------
void Geometry::reverse_sort_from (std::vector<entryexit> & ee, const vec point) 
{
  std::vector<double> dst1 ; 
  std::vector<double> dst2 ;
  
  for (size_t i=0 ; i<ee.size() ; i++)
  {
    ee[i].dst1 = dstsqr(ee[i].pt1, point) ; 
    ee[i].dst2 = dstsqr(ee[i].pt2, point) ; 
  }
  std::sort(ee.begin(), ee.end(), [](auto i, auto j){return i.dst1>j.dst2;}) ; 
  
}
//-----------------------------------------------------------------------
std::optional<double> Geometry::intersection_ray_triangle_mollertrumbore (const Triangle & tri, const Eigen::Vector3d &x0, const Eigen::Vector3d &xr)
{
  const double EPSILON = 1e-8;

  //Vector3d edge1 = v1 - v0;
  //Vector3d edge2 = v2 - v0;

  Eigen::Vector3d h = xr.cross(tri.v);
  double a = tri.u.dot(h);

  if (fabs(a) < EPSILON) return std::nullopt; // ray coplanar with tri

  double f = 1.0 / a;
  Eigen::Vector3d s = x0 - tri.M;
  double u = f * s.dot(h);

  if (u < 0.0 || u > 1.0) return std::nullopt;

  Eigen::Vector3d q = s.cross(tri.u);
  double v = f * xr.dot(q);

  if (v < 0.0 || u + v > 1.0) return std::nullopt;

  double t = f * tri.v.dot(q);

  if (t > EPSILON) {return t ; }

  return std::nullopt; // No intersection
}
//-----------------------------------------------------------------------
std::optional<double> Geometry::intersection_ray_triangle_inversion (const Triangle & tri, const Eigen::Vector3d &x0, const Eigen::Vector3d &xr)
{
  Eigen::Matrix3d A ; 
  //Eigen::Vector3d b = x0-tri.M ; 
  A << -xr, tri.u, tri.v ; 
  Eigen::Vector3d sol = A.colPivHouseholderQr().solve(x0-tri.M);
  //double relative_error = (A*sol - b).norm() / b.norm();
  if (sol[1]>0 && sol[1]<1 && sol[2]>0 && sol[2]<1 && sol[1]+sol[2]<1)
    return sol[0] ; 
  return {} ; 
}

//-----------------------------------------------------------------------
char Geometry::intersection_ray_tetra (std::vector<std::tuple<double, double, int>> & inter, const Tetrahedron & tetra, int tetraid, const Eigen::Vector3d &x0, const Eigen::Vector3d &xr)
{
  std::vector<double> contacts ; 
  char tri_intersect = 0 ; 
   std::optional<double> ct ; 
  for (int i=0 ; i<4 ; i++)
  {
    if ( static_cast<int>(Parameters.strategy) & 1)
      ct=intersection_ray_triangle_inversion(tetra.triangles[i], x0, xr) ; 
    else 
      ct=intersection_ray_triangle_mollertrumbore(tetra.triangles[i], x0, xr) ; 
    //std::optional<double> ct = intersection_ray_triangle_B(tetra.triangles[i], x0, xr) ; 
    //if (ct.has_value() || ctB.has_value())
    //  printf("%g %g | ", ct.has_value()?ct.value():0, ctB.has_value()?ctB.value():0) ; 
    if (ct.has_value())
    {
      contacts.push_back(ct.value()) ; 
      tri_intersect |= 1<<i ; 
    }
  }
  if (contacts.size()!= 0  &&  contacts.size()!=2) 
  {
    printf("ERR: Expect 0 or 2 contacts with tetrahedron %ld \n", contacts.size()) ;
  }
  
  if (contacts.size()==2 && fabs(contacts[0]-contacts[1])>1e-8)
  {
    if (contacts[0]<contacts[1])
    {
      inter.push_back({contacts[0], contacts[1], tetraid}) ;
      return tri_intersect ; 
    }
    else      
    {
      inter.push_back({contacts[1], contacts[0], tetraid}) ;
      return tri_intersect ; 
    }
  }
  return 0;   
}

//-----------------------------------------------------------------------
void Geometry::intersection_ray_mesh (std::vector<std::tuple<double, double, int>> & intersect, const std::vector<Tetrahedron> & mesh, const Eigen::Vector3d &x0, const Eigen::Vector3d &xr)
{
  // TODO handle mesh offset (probably by offseting the ray, actually. 
  intersect.clear() ; 
  char contact=0 ; 
  for (size_t i=0 ; i<mesh.size() && !contact ; i++)
  {
    if (mesh[i].surface == 0) continue ; 
    contact = intersection_ray_tetra(intersect, mesh[i], i, x0, xr) ; 
    if (contact && ((contact & mesh[i].surface) == 0) )
    {
      intersect.pop_back() ;
      contact=0 ; 
    }
  }
  if (intersect.size()==0) return ; // No intersection with any of the surface tetrahedras
  
  int curtetra, newtetra ; 
  std::set<int> alreadyfound ; alreadyfound.insert(std::get<2>(intersect.back())) ; 
  do {
    curtetra=std::get<2>(intersect.back()) ; 
    for (auto i : mesh[curtetra].neigh)
    {
      if (alreadyfound.contains(i)) continue ;
      contact = intersection_ray_tetra(intersect, mesh[i], i, x0, xr) ;  
      if (contact)
      {
        if (fabs(std::get<0>(intersect[intersect.size()-2])-std::get<1>(intersect[intersect.size()-1])) > 1e-8
          && fabs(std::get<1>(intersect[intersect.size()-2])-std::get<0>(intersect[intersect.size()-1])) > 1e-8 )
        {
          intersect.pop_back() ; 
        }
        else
        {
          alreadyfound.insert(std::get<2>(intersect.back())) ; 
          break ; 
        }
      }
    }
    newtetra = std::get<2>(intersect.back()) ; 
    if (newtetra==curtetra)
    {
      for (auto i : mesh[curtetra].neigh_extended)
      {
        if (alreadyfound.contains(i)) continue ;
        contact = intersection_ray_tetra(intersect, mesh[i], i, x0, xr) ; 
        if (contact)
        {
          if (fabs(std::get<0>(intersect[intersect.size()-2])-std::get<1>(intersect[intersect.size()-1])) > 1e-8 &&
            fabs(std::get<1>(intersect[intersect.size()-2])-std::get<0>(intersect[intersect.size()-1])) > 1e-8)
          {
            intersect.pop_back() ; 
          }
          else
          {
            alreadyfound.insert(std::get<2>(intersect.back())) ; 
            break ; 
          }
        }
      }
      newtetra = std::get<2>(intersect.back()) ; 
      if (newtetra==curtetra)
      {
        printf("ERR: couldn't find an adjacent tetra traversed by the ray. May or may not be a problem\n") ; 
        break ; 
      }
    }
  } while( (contact & mesh[newtetra].surface) == 0) ;
  
  //for (int i=0 ; i<mesh.size() ; i++)
  //  intersection_ray_tetra(intersect, mesh[i], i, x0, xr) ; 
  
  std::sort(intersect.begin(), intersect.end(), [](auto i, auto j){return std::get<0>(i)<std::get<0>(j) ; }) ;
   
  //for (auto &v: intersect) printf("%g %g %d\n", std::get<0>(v), std::get<1>(v), std::get<2>(v)) ; 
  
  for (size_t i=1 ; i<intersect.size() ; i++)
    if (fabs(std::get<1>(intersect[i-1])-std::get<0>(intersect[i])) > 1e-8)
    {
      printf("\n") ; 
      for (auto &v: intersect) 
      {
        printf("%g %g %d | ", std::get<0>(v), std::get<1>(v), std::get<2>(v)) ; 
        for (auto &w: mesh[std::get<2>(v)].neigh)
          printf("%d ", w) ; 
        printf("\n") ; 
      }
      printf("ERR: too much distance between consecutive exit-entry points in tetrahedras %g %ld\n", fabs(std::get<1>(intersect[i-1])-std::get<0>(intersect[i])), i) ; 
      std::exit(0) ; 
    }
    
  return ;
}
