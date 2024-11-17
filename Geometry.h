#include <vector>
#include <optional>
#include <algorithm>
#include <complex>
#include <set>
#include "Eigen/Dense"

using namespace std::complex_literals ; 

#ifndef GEOMETRY_H
#define GEOMETRY_H

class vec {
public :
  vec() = default ; 
  vec(double a, double b, double c) {v[0]=a ; v[1]=b; v[2]=c;}
  double & operator[] (int i) { return v[i] ; }
  std::vector<double> v{0,0,0} ; 
  void normalise() {double l=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]) ; for (int i=0 ; i<3 ; i++) v[i]/=l ; }
} ; 

class vec_stokes {
public :
  vec_stokes() = default ; 
  vec_stokes(double a, double b, double c, double d) {v[0]=a ; v[1]=b; v[2]=c; v[3]=d ; }
  double & operator[] (int i) { return v[i] ; }
  std::vector<double> v{1,0,0,0} ;
  void disp() {for (auto &vv:v) printf("%g ", vv) ;  printf("\n") ; }
  void normalise_coherent()
  {
    double norm = sqrt(v[1]*v[1]+v[2]*v[2]+v[3]*v[3])/v[0] ; 
    v[1]/=norm ; 
    v[2]/=norm ; 
    v[3]/=norm ;     
  }
} ; 

class vec_jones {
public:
  vec_jones() = default ; 
  vec_jones(std::complex<double> a, std::complex<double> b) {v[0]=a ; v[1]=b;}
  std::complex<double> & operator[] (int i) { return v[i] ; }
  void disp() {for (auto &vv:v) std::cout << vv << " " ;  printf("\n") ; }
  std::complex<double> v[2] ;   
  void normalise_coherent()
  {
    double n = sqrt(norm(v[0])+norm(v[1])) ;
    v[0] /= n ; 
    v[1] /= n ; 
  }
} ; 



class tens {
public:
  tens() = default ; 
  tens(double a, double b, double c, double d, double e, double f, double g, double h, double i)
  {
    v[0] = a ; v[1] = b ; v[2] = c;
    v[3] = d ; v[4] = e ; v[5] = f;
    v[6] = g ; v[7] = h ; v[8] = i;
  }
  double & operator[] (int i) { return v[i] ; }
  std::vector<double> v{0,0,0,0,0,0,0,0,0} ; 
  double norm() {return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+
                                 v[3]*v[3]+v[4]*v[4]+v[5]*v[5]+
                                 v[6]*v[6]+v[7]*v[7]+v[8]*v[8]) ; }
  void disp() {for (auto &vv:v) printf("%g ", vv) ;  printf("\n") ; }
  void swap_row (int i, int j) 
  {
    double tmp ; 
    tmp = v[i*3+0] ; v[i*3+0]=v[j*3+0] ; v[j*3+0]=tmp ; 
    tmp = v[i*3+1] ; v[i*3+1]=v[j*3+1] ; v[j*3+1]=tmp ; 
    tmp = v[i*3+2] ; v[i*3+2]=v[j*3+2] ; v[j*3+2]=tmp ; 
  }
  /*void gauss_pivot()
  {
    std::array<int 3> p={0,1,2} ; 
    
    if (v[0]<v[3])
    {
      if (v[3]<v[6])
      {
        swap_row(2,0) ; 
        p[0] = 2 ; p[2] = 0 ; 
      }
      else 
      {
        swap_row(1,0) ;
        p[0] = 1 ; p[1] = 0 ;
      }
    }
    else //v[0]>v[3]
      if (v[0]<v[6])
      {
        swap_row(2,0) ;
        p[0] = 2 ; 
        p[2] = 0 ;
      }
    
    v[4] = v[4] - v[3]/v[0]*v[1] ; 
    v[5] = v[5] - v[3]/v[0]*v[2] ; 
    v[3] = 0 ; 
    v[7] = v[7] - v[6]/v[0]*v[1] ; 
    v[8] = v[8] - v[6]/v[0]*v[2] ; 
    v[6] = 0 ; 
    
    if (v[4]<v[5])
    {
       swap_row(2,1) ;
       int tmp ; tmp = p[1] ; p[1] = p[2] ; p[2] = tmp ; }
    }
    v[8] = v[8] - v[7]/v[4]*v[5] ; 
    v[7] = 0 ; 
    
    return p ;    
  }*/
  
} ; 


class mat_stokes{
public:
  mat_stokes() = default ; 
  mat_stokes(std::array<double,16> val)
  {
    for (int i=0 ; i<16 ; i++)
      v[i]=val[i] ;     
  }
  double & operator[] (int i) { return v[i] ; }
  std::vector<double> v{1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1} ; 
  void disp() {for(int i=0 ; i<4 ; i++) { for (int j=0 ; j<4 ; j++) printf("%g ", v[i*4+j]) ; printf("\n") ; }}
} ;

class mat_jones{
public:
  mat_jones() = default ; 
  mat_jones(std::array<std::complex<double>,4> val)
  {
    for (int i=0 ; i<4 ; i++)
      v[i]=val[i] ;     
  }
  std::complex<double> & operator[] (int i) { return v[i] ; }
  std::complex<double> v[4]={1,0, 0,1} ; 
  void disp() {for(int i=0 ; i<4 ; i++) { for (int j=0 ; j<4 ; j++) std::cout << v[i*4+j] << " " ; printf("\n") ; }}
} ;

namespace Polariser {
  static mat_stokes horizontal_lin({0.5, 0.5,0,0,  0.5,0.5,0,0, 0,0,0,0, 0,0,0,0}) ; 
  static mat_stokes vertical_lin({0.5,-0.5,0,0, -0.5,0.5,0,0, 0,0,0,0, 0,0,0,0}) ; 
  static mat_stokes deg45_lin({0.5,0,0.5,0, 0,0,0,0, 0.5,0,0.5,0, 0,0,0,0}) ; 
  
  /*static vec_jones horiz_jones({1,0}) ; 
  static vec_jones vert_jones({0,1}) ; 
  static vec_jones circright_jones({1/sqrt(2), -1i / sqrt(2)}) ; 
  static vec_jones circleft_jones({1/sqrt(2), -1i / sqrt(2)}) ; */
  
  static mat_jones horiz_jones({1,0,0,0}) ; 
  static mat_jones vert_jones({0,0,0,1}) ; 
  //static mat_jones circright_jones({1/sqrt(2), -1i / sqrt(2)}) ; 
  //static mat_jones circleft_jones({1/sqrt(2), -1i / sqrt(2)}) ; 
  
} ; 

class Triangle
{
public: 
  Eigen::Vector3d M, u, v ; // One of the triangle vertex, and the 2 side vectors from that vertex. 
  void disp() const {std::cout << M.transpose() << "|" << u.transpose() << "|" << v.transpose() << "\n";}
} ; 
class Tetrahedron
{ 
public:
  Triangle triangles[4] ; 
  Eigen::Vector3d com ; 
  std::array<int,4> vertex_id ; 
  std::set<int> neigh ;  
  std::set<int> neigh_extended ;  
  unsigned char surface = 0 ; 
  
  //----------------------------------
  char identify_tri (std::vector<int>& vids)
  {
    //printf("%d %d %d %d | %d %d %d\n", vertex_id[0], vertex_id[1], vertex_id[2], vertex_id[3], vids[0], vids[1], vids[2]) ; 
    if (vertex_id[0]==vids[0] || vertex_id[0]==vids[1] || vertex_id[0]==vids[2])
    {
      if (vertex_id[1]==vids[0] || vertex_id[1]==vids[1] || vertex_id[1]==vids[2])
      {
        if (vertex_id[2]==vids[0] || vertex_id[2]==vids[1] || vertex_id[2]==vids[2])
          return 0 ;
        else if (vertex_id[3]==vids[0] || vertex_id[3]==vids[1] || vertex_id[3]==vids[2])
          return 1 ;
        else
         printf("Inconsistent triangle finding 01?\n") ; 
      }
      else if (vertex_id[2]==vids[0] || vertex_id[2]==vids[1] || vertex_id[2]==vids[2])
      {
        if (vertex_id[3]==vids[0] || vertex_id[3]==vids[1] || vertex_id[3]==vids[2])
          return 2;
        else
          printf("Inconsistent triangle finding 02?\n") ; 
      }
      else 
        printf("Inconsistent triangle finding 0??\n") ; 
    }
    else if (vertex_id[1]==vids[0] || vertex_id[1]==vids[1] || vertex_id[1]==vids[2])
    {
      if (vertex_id[2]==vids[0] || vertex_id[2]==vids[1] || vertex_id[2]==vids[2])
      {
        if (vertex_id[3]==vids[0] || vertex_id[3]==vids[1] || vertex_id[3]==vids[2])
          return 3 ; 
        else
          printf("Inconsistent triangle finding 12?\n") ; 
      }
      else
        printf("Inconsistent triangle finding 1??\n") ; 
    }
    else
      printf("Inconsistent triangle finding ???\n") ; 
    return -1 ; 
  }
  
  //----------------------------------
  void write(FILE * out)
  {
    auto eigenwrite = [&](Eigen::Vector3d & v)
    {
      fwrite(&(v[0]), sizeof(double), 1, out) ; 
      fwrite(&(v[1]), sizeof(double), 1, out) ; 
      fwrite(&(v[2]), sizeof(double), 1, out) ;      
    } ; 
    
    for (int i=0 ; i<4 ; i++)
    {
      eigenwrite(triangles[i].M) ; 
      eigenwrite(triangles[i].u) ; 
      eigenwrite(triangles[i].v) ; 
    }
    eigenwrite(com) ;
    fwrite(&surface, sizeof(unsigned int), 1, out) ;    
    size_t tmp = neigh.size() ; 
    fwrite(&tmp, sizeof(size_t), 1, out) ; 
    for (auto v: neigh) 
      fwrite(&v, sizeof(int), 1, out) ;    
    tmp = neigh_extended.size() ; 
    fwrite(&tmp, sizeof(size_t), 1, out) ; 
    for (auto v: neigh_extended) 
      fwrite(&v, sizeof(int), 1, out) ; 
  }
  //----------------------------------
  void read(FILE * out)
  {
    size_t res ; 
    auto eigenread = [&](Eigen::Vector3d & v)
    {
      res=fread(&(v[0]), sizeof(double), 1, out) ; 
      res=fread(&(v[1]), sizeof(double), 1, out) ; 
      res=fread(&(v[2]), sizeof(double), 1, out) ;      
    } ; 
    
    for (int i=0 ; i<4 ; i++)
    {
      eigenread(triangles[i].M) ; 
      eigenread(triangles[i].u) ; 
      eigenread(triangles[i].v) ; 
    }
    eigenread(com) ; 
    fread(&surface, sizeof(unsigned int), 1, out) ;   
    size_t tmp ; 
    res=fread(&tmp , sizeof(size_t), 1, out) ; 
    for (size_t i=0 ; i<tmp ; i++)
    {
      int ii ; 
      res=fread(&ii, sizeof(int), 1, out) ;
      neigh.insert(ii) ; 
    }
    res=fread(&tmp , sizeof(size_t), 1, out) ;
    for (size_t i=0 ; i<tmp ; i++)
    {
      int ii ; 
      res=fread(&ii, sizeof(int), 1, out) ;  
      neigh_extended.insert(ii) ; 
    }
  }
} ;


inline double dstsqr (vec pt1, vec pt2) {return (pt1[0]-pt2[0])*(pt1[0]-pt2[0]) + (pt1[1]-pt2[1])*(pt1[1]-pt2[1]) + (pt1[2]-pt2[2])*(pt1[2]-pt2[2]) ; }
inline vec operator* (double s, vec v) {return {s*v[0],s*v[1],s*v[2]};}
inline vec operator+ (vec v1, vec v2) {return {v1[0]+v2[0],v1[1]+v2[1],v1[2]+v2[2]};}
inline vec operator- (vec v) {return {-v[0], -v[1], -v[2]} ; }
inline vec operator- (vec v1, vec v2) {return {v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]} ; }
inline tens operator* (tens & a, tens & b) 
{
  return tens(a[0]*b[0]+a[1]*b[3]+a[2]*b[6], a[0]*b[1]+a[1]*b[4]+a[2]*b[7], a[0]*b[2]+a[1]*b[5]+a[2]*b[8],
          a[3]*b[0]+a[4]*b[3]+a[5]*b[6], a[3]*b[1]+a[4]*b[4]+a[5]*b[7], a[0]*b[2]+a[4]*b[5]+a[5]*b[8],
          a[6]*b[0]+a[7]*b[3]+a[8]*b[6], a[6]*b[1]+a[7]*b[4]+a[8]*b[7], a[0]*b[2]+a[7]*b[5]+a[8]*b[8]) ; 
}
inline vec_stokes operator* (mat_stokes & a, vec_stokes & b) 
{
  return {a[0 ]*b[0]+a[1 ]*b[1]+a[2 ]*b[2]+a[3 ]*b[3], 
          a[4 ]*b[0]+a[5 ]*b[1]+a[6 ]*b[2]+a[7 ]*b[3], 
          a[8 ]*b[0]+a[9 ]*b[1]+a[10]*b[2]+a[11]*b[3], 
          a[12]*b[0]+a[13]*b[1]+a[14]*b[2]+a[15]*b[3]}; 
}
inline vec_jones operator* (mat_jones & a, vec_jones & b) 
{
  return {a[0 ]*b[0]+a[1 ]*b[1],
          a[2 ]*b[0]+a[3 ]*b[1]}; 
}

class entryexit {
public:
  vec pt1, pt2 ; 
  int objectid=-1; 
  double dst1, dst2 ;
  
  std::pair<std::vector<vec>, double> locations_inbetween (int N)
  {
    std::vector<vec> res ; 
    res.resize(N) ; 
    for (int i=0 ; i<N ; i++)
      res[i] = pt1 + (i/(double)N)*(pt2-pt1) ; 
    double ds = sqrt(dstsqr(pt2,pt1))/N ; 
    return {res,ds} ; 
  }
  
} ; 


namespace Geometry {
std::optional<entryexit> intersection_ray_sphere (vec x0, vec xr, vec xc, double r) 
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
  
void reverse_sort_from (std::vector<entryexit> & ee, const vec point) 
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

#include <iostream>
#include <vector>
#include <optional>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

// Function to compute intersection of a ray and a triangle
std::optional<double> intersection_ray_triangle_A (const Triangle & tri, const Eigen::Vector3d &x0, const Eigen::Vector3d &xr)
{
  const double EPSILON = 1e-8;

  //Vector3d edge1 = v1 - v0;
  //Vector3d edge2 = v2 - v0;

  Vector3d h = xr.cross(tri.v);
  double a = tri.u.dot(h);

  if (fabs(a) < EPSILON) return nullopt; // ray coplanar with tri

  double f = 1.0 / a;
  Vector3d s = x0 - tri.M;
  double u = f * s.dot(h);

  if (u < 0.0 || u > 1.0) return nullopt;

  Vector3d q = s.cross(tri.u);
  double v = f * xr.dot(q);

  if (v < 0.0 || u + v > 1.0) return nullopt;

  double t = f * tri.v.dot(q);

  if (t > EPSILON) {return t ; }

  return nullopt; // No intersection
}
//-----------------------------------------------------------------
std::optional<double> intersection_ray_triangle_B (const Triangle & tri, const Eigen::Vector3d &x0, const Eigen::Vector3d &xr)
{
  Eigen::Matrix3d A ; 
  Eigen::Vector3d b = x0-tri.M ; 
  A << -xr, tri.u, tri.v ; 
  Eigen::Vector3d sol = A.colPivHouseholderQr().solve(x0-tri.M);
  double relative_error = (A*sol - b).norm() / b.norm();
  if (sol[1]>0 && sol[1]<1 && sol[2]>0 && sol[2]<1 && sol[1]+sol[2]<1)
    return sol[0] ; 
  return {} ; 
}

char intersection_ray_tetra (std::vector<std::tuple<double, double, int>> & inter, const Tetrahedron & tetra, int tetraid, const Eigen::Vector3d &x0, const Eigen::Vector3d &xr)
{
  std::vector<double> contacts ; 
  char tri_intersect = 0 ; 
  for (int i=0 ; i<4 ; i++)
  {
    std::optional<double> ct = intersection_ray_triangle_A(tetra.triangles[i], x0, xr) ; 
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


void intersection_ray_mesh (std::vector<std::tuple<double, double, int>> & intersect, const std::vector<Tetrahedron> & mesh, const Eigen::Vector3d &x0, const Eigen::Vector3d &xr)
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


}
#endif









