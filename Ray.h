#include "Geometry.h"

class Ray {
public: 
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
    /*for (int i=0 ; i<dphis.size() ; i++) printf("%g ", dphis[i]) ; */
    for (size_t i=0 ; i<ddeltas.size() ; i++)
    {
      s2tmp = polarisation[2] ; 
      polarisation[2] += -2*dphis[i]*polarisation[1] - photoconstant*ddeltas[i]*polarisation[3] ; 
      polarisation[1] += 2*dphis[i]*s2tmp ; 
      polarisation[3] += photoconstant*ddeltas[i]*s2tmp ;   
      //printf("%g %g ", dphis[i], ddeltas[i]) ; polarisation.disp() ; 
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

