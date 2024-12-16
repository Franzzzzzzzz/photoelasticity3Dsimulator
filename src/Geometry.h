#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include <optional>
#include <algorithm>
#include <complex>
#include <set>
#include <iostream>
#include "Eigen/Dense"
#include "Algebra.h"
#include "Parameters.h"
extern struct Parameters_ Parameters ; 

using namespace std::complex_literals ; 


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
    res=fread(&surface, sizeof(unsigned int), 1, out) ;   
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
std::optional<entryexit> intersection_ray_sphere (vec x0, vec xr, vec xc, double r) ;   
void reverse_sort_from (std::vector<entryexit> & ee, const vec point) ; 

// Function to compute intersection of a ray and a triangle
std::optional<double> intersection_ray_triangle_mollertrumbore (const Triangle & tri, const Eigen::Vector3d &x0, const Eigen::Vector3d &xr) ; 
std::optional<double> intersection_ray_triangle_inversion (const Triangle & tri, const Eigen::Vector3d &x0, const Eigen::Vector3d &xr) ; 
char intersection_ray_tetra (std::vector<std::tuple<double, double, int>> & inter, const Tetrahedron & tetra, int tetraid, const Eigen::Vector3d &x0, const Eigen::Vector3d &xr) ;
void intersection_ray_mesh (std::vector<std::tuple<double, double, int>> & intersect, const std::vector<Tetrahedron> & mesh, const Eigen::Vector3d &x0, const Eigen::Vector3d &xr) ;
}
#endif









