#include "Geometry.h"


class Ray {
public: 
  vec direction ; 
  vec destination ; 
  vec_stokes polarisation ; 
  
  double extra_value=0 ; 
  
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
  void propagate (std::vector<double> ddeltas, std::vector<double> dphis, bool print=false) 
  {
    double s2tmp ; 
    //polarisation.disp() ; 
    if (print) for (int i=0 ; i<ddeltas.size() ; i++) printf("%g ", ddeltas[i]) ; 
    /*for (int i=0 ; i<dphis.size() ; i++) printf("%g ", dphis[i]) ; */
    for (size_t i=0 ; i<ddeltas.size() ; i++)
    {
      if (print) polarisation.disp() ; 
      s2tmp = polarisation[2] ; 
      polarisation[2] += -2*dphis[i]*polarisation[1] - ddeltas[i]*polarisation[3] ; 
      polarisation[1] += 2*dphis[i]*s2tmp ; 
      polarisation[3] += ddeltas[i]*s2tmp ;   
      //printf("%g %g ", dphis[i], ddeltas[i]) ; polarisation.disp() ; 
      polarisation.normalise_coherent() ; 
    }    
    //polarisation.disp() ; 
  }
  
  void apply_polariser (mat_stokes pol)
  {
    polarisation = pol * polarisation; 
  }
  
  
} ; 

