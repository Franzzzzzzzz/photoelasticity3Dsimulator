#include <vector>
#include <optional>
#include <algorithm>

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

namespace Polariser {
  static mat_stokes horizontal_lin({0.5, 0.5,0,0,  0.5,0.5,0,0, 0,0,0,0, 0,0,0,0}) ; 
  static mat_stokes vertical_lin({0.5,-0.5,0,0, -0.5,0.5,0,0, 0,0,0,0, 0,0,0,0}) ; 
  static mat_stokes deg45_lin({0.5,0,0.5,0, 0,0,0,0, 0.5,0,0.5,0, 0,0,0,0}) ; 
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

/*std::optional intersection_ray_triangle (vec x0, vec xr, vec xt, vec u, vec v)
{
  
}*/

}
#endif









