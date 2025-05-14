#ifndef RAY_H
#define RAY_H

#include "Geometry.h"
#include "Grains.h"
#include "unsupported/Eigen/MatrixFunctions"

class Ray {
public: 
  vec direction ; 
  vec destination ; 
  
  // Datas along ray
  double extra_value=0 ; 
  vec_jones B, Bout ; 
  std::complex<double> M[4] = {1,0,0,1} ; // Transformation matrix
  double path_length = 0 ; 
  
  tens get_rotation_matrix() 
  {
    tens res ; 
    double angle = atan2(direction[1], direction[0]) ; 
    res[0*3 + 0] = res[1*3 + 1] = cos(angle) ; 
    res[1*3 + 0] = sin(angle) ; 
    res[0*3 + 1] = -sin(angle) ; 
    res[2*3 + 2] =1 ; 
    return res ; 
  }
  tens get_inv_rotation_matrix() 
  {
    tens res ; 
    double angle = atan2(direction[1], direction[0]) ; 
    res[0*3 + 0] = res[1*3 + 1] = cos(angle) ; 
    res[1*3 + 0] = -sin(angle) ; 
    res[0*3 + 1] = sin(angle) ; 
    res[2*3 + 2] =1 ; 
    return res ; 
  }
  
  void set_direction(vec dir) {direction=dir ; }
  void set_destination(vec dest) {destination=dest; }
  void set_polarisation(vec_jones pol) {B = pol ; }
  double get_intensity() {return abs(Bout[0])*abs(Bout[0])+abs(Bout[1])*abs(Bout[1]) ; }
  
  std::vector<entryexit> find_grains_entryexit(const std::vector<Grains> & grains) 
  {
    std::vector<entryexit> res ; 
    for (size_t i = 0 ; i < grains.size() ; i++)
    {
      auto intersection = Geometry::intersection_ray_sphere (destination, direction, grains[i].center, grains[i].r) ;  
      if (intersection.has_value())
      {
        intersection.value().objectid = i ; 
        res.push_back(intersection.value()) ; 
      }
    }
    return res ;     
  }
  
  void propagate (std::vector<tens> &sigma, double C, double ds)
  {
    std::complex<double> Mnew[4] ; 
        
    for (size_t i=0 ; i<sigma.size() ; i++)
    {
      std::complex<double> Mtmp[4] = {1,0,0,1} ; 
      Mtmp[0]  = 1. - 1i * C * ds * (sigma[i][0]-sigma[i][8]) ; // Unclear if it should be the stress components (0,0) and (0,2) ; or (1,1) and (1,2) ...
      Mtmp[1]  =    - 1i * C * ds * (sigma[i][2]) ; 
      Mtmp[2]  =    - 1i * C * ds * (sigma[i][2]) ; 
      Mtmp[3]  = 1. + 1i * C * ds * (sigma[i][0]-sigma[i][8]); 
      
      Mnew[0] = Mtmp[0] * M[0] + Mtmp[1]*M[2] ; 
      Mnew[1] = Mtmp[0] * M[1] + Mtmp[1]*M[3] ; 
      Mnew[2] = Mtmp[2] * M[0] + Mtmp[3]*M[2] ; 
      Mnew[3] = Mtmp[2] * M[1] + Mtmp[3]*M[3] ; 
      
      M[0] = Mnew[0] ; 
      M[1] = Mnew[1] ; 
      M[2] = Mnew[2] ; 
      M[3] = Mnew[3] ; 
      
      //std::complex<double> dB0 = - 1i * C * ((sigma[i][4]-sigma[i][8]) * B[0] + sigma[i][5] * B[1]) ;
      //std::complex<double> dB1 = - 1i * C * ( sigma[i][5] * B[0]  - (sigma[i][4]-sigma[i][8]) * B[1] ) ; 
      
      //std::cout << dB0 << " " << dB1 << "\n" ;  
      
      //B[0] += dB0 ; 
      //B[1] += dB1 ;   
      //B.normalise_coherent() ; 
    }
  }
  
  void apply_propagation()
  { 
    Bout[0] = M[0]*B[0] + M[1]*B[1] ; 
    Bout[1] = M[2]*B[0] + M[3]*B[1] ; 
    Bout.normalise_coherent() ; 
  }
  void apply_absorption(double absorption)
  { 
    Bout[0] = B[0]*exp(-absorption*path_length) ; 
    Bout[1] = B[1]*exp(-absorption*path_length) ; 
  }
  
  
  void propagate_exp (std::vector<tens> &sigma, double C, std::vector<double> ds)
  {
    Eigen::Matrix2cd M({{1., 0.},{0., 1.}}) ; 
    Eigen::Matrix2cd Mtmp ; 
    
    for (size_t i=0 ; i<ds.size() ; i++)
    {
      Mtmp(0,0) = - 1i * C * ds[i] * (sigma[i][0]-sigma[i][8]) ; 
      Mtmp(0,1) = - 1i * C * ds[i] * (sigma[i][2]) ;  
      Mtmp(1,0) = - 1i * C * ds[i] * (sigma[i][2]) ;  
      Mtmp(1,1) = + 1i * C * ds[i] * (sigma[i][0]-sigma[i][8]);  
      M = M*(Mtmp.exp()) ; 
    }
    
    Bout[0] = M(0,0)*B[0]+M(0,1)*B[1] ; 
    Bout[1] = M(1,0)*B[0]+M(1,1)*B[1] ; 
    Bout.normalise_coherent() ;     
  }

  
  void apply_polariser (mat_jones pol)
  {
    Bout = pol * Bout; 
  }
  
} ; 


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
// Ray Stokes, kept for history, do not use.
class RayJones 
{
private: 
  vec direction ; 
  vec destination ; 
  vec_stokes polarisation ; 
  
  // Datas along ray
  double extra_value=0 ; 
  std::vector<double> ddeltas, dphis ; 
  std::vector<double> phis ; 
  
  tens get_rotation_matrix() 
  {
    tens res ; 
    double angle = atan2(direction[1], direction[0]) ; 
    res[0*3 + 0] = res[1*3 + 1] = cos(angle) ; 
    res[1*3 + 0] = sin(angle) ; 
    res[0*3 + 1] = -sin(angle) ; 
    res[2*3 + 2] =1 ; 
    return res ; 
  }
  tens get_inv_rotation_matrix() 
  {
    tens res ; 
    double angle = atan2(direction[1], direction[0]) ; 
    res[0*3 + 0] = res[1*3 + 1] = cos(angle) ; 
    res[1*3 + 0] = -sin(angle) ; 
    res[0*3 + 1] = sin(angle) ; 
    res[2*3 + 2] =1 ; 
    return res ; 
  }
  
  void set_direction(vec dir) {direction=dir ; }
  void set_destination(vec dest) {destination=dest; }
  void set_polarisation(vec_stokes pol) {polarisation = pol ; }
  double get_intensity () {return polarisation[0] ; }
  
  std::vector<entryexit> find_grains_entryexit(const std::vector<Grains> & grains) 
  {
    std::vector<entryexit> res ; 
    for (size_t i = 0 ; i < grains.size() ; i++)
    {
      auto intersection = Geometry::intersection_ray_sphere (destination, direction, grains[i].center, grains[i].r) ;  
      if (intersection.has_value())
      {
        intersection.value().objectid = i ; 
        res.push_back(intersection.value()) ; 
      }
    }
    return res ;     
  }
  
  void absorbe (double absorption, double ds)
  {
    polarisation[0] *= absorption*ds ; 
  }
  
  void propagate (double photoconstant) 
  {
    double s2tmp ; 
    //polarisation.disp() ; 
    for (size_t i=0 ; i<ddeltas.size() ; i++)
    {
      s2tmp = polarisation[2] ; 
      polarisation[2] += -2*dphis[i]*polarisation[1] - photoconstant*ddeltas[i]*polarisation[3] ; 
      polarisation[1] += 2*dphis[i]*s2tmp ; 
      polarisation[3] += photoconstant*ddeltas[i]*s2tmp ;   
      polarisation.normalise_coherent() ; 
    }    
    //polarisation.disp() ; 
  }
  
  void apply_polariser (mat_stokes pol)
  {
    polarisation = pol * polarisation; 
  }
  
  //---------------------------------------------
  void get_photoelastic_deltas(std::vector<tens> &sigma, double ds)
  {
    ddeltas.resize(sigma.size()) ; 
    for (size_t i = 0 ; i<ddeltas.size(); i++)
    {
      ddeltas[i] = sqrt( (sigma[i][4]-sigma[i][8])*(sigma[i][4]-sigma[i][8]) + 4 * sigma[i][5] * sigma[i][5] ) * 
               ds ;
    }
  }
  
  //-----------------------------
  std::vector<double> get_photoelastic_dphis(std::vector<tens> &sigma)
  {
    phis.resize(sigma.size()) ; 
    dphis.resize(sigma.size()) ; 
    for (size_t i=0 ; i<phis.size() ; i++)
    { 
      //phis[i] = 0.5*atan(2*sigma[i][5]/(sigma[i][4]-sigma[i][8])) ;
      phis[i] = 0.5*atan2(sigma[i][5], (sigma[i][4]-sigma[i][8])/2) ; 
    }
    
    for (size_t i=1 ; i<phis.size()-1 ; i++)
    {
      dphis[i] = (phis[i+1]-phis[i-1])/2 ; 
      if (dphis[i]>=M_PI/2) dphis[i]-=M_PI ; 
      if (dphis[i]<=-M_PI/2) dphis[i]+=M_PI ; 
      //dphis[i] /= (2.*ds) ;
    }
    dphis[0] = (dphis[1]-dphis[0]);///ds ; 
    dphis[dphis.size()-1] = (dphis[dphis.size()-1]-dphis[dphis.size()-2]);///ds ; 
    
    return dphis ;  
  }  
  
} ; 

#endif
