#ifndef IMAGE_H
#define IMAGE_H

#include <vector>
#include <SDL2/SDL.h>

#include "Parameters.h"
#include "Geometry.h"
#include "Ray.h"
#include "Fem.h"

extern struct Parameters_ Parameters ; 

class Image {
public: 
  Image (int W,int H) { I.resize(H, std::vector<double>(W,0)) ; extra = I ; w=W ; h=H ; }
  std::vector<std::vector<double>> I ; 
  std::vector<std::vector<double>> extra ; 
  std::vector<Ray> rays ; 
  double dx, dz ;
  int w, h ; 
  vec origin ; 
  vec normal ;
  
  void set_dxdz(double ddx, double ddz) {dx=ddx ; dz=ddz ; }
  void set_origin (vec v) {origin = v ; }
  void set_normal (vec v) {normal = v ; }
  void set_origin (double distance, double angle) {origin = {distance*cos(angle), distance*sin(angle),0.}; }
  void set_normal (double angle) { normal = {-cos(angle), -sin(angle), 0.} ; }
  int size () {return w*h ;}
  
  vec get_location(int i) ;
  void set_pixel(int i, int j, double intensity) {I[i][j]=intensity ; }
  void set_pixel(int i, double intensity) {set_pixel(i/w, i%w, intensity) ; }
  void set_extra_value(int i, int j, double v) {extra[i][j]=v ;  }
  void set_extra_value(int i, double v) {set_extra_value(i/w, i%w, v) ; }
  void set_rays () 
  {
    rays.resize(size()) ; 
    for (int i=0 ; i<size() ; i++)
    {
      rays[i].set_direction(normal) ; 
      rays[i].set_destination(get_location(i)) ; 
      rays[i].set_polarisation (Parameters.polarisation) ; 
    }
  }
  void reset_rays()
  {
    rays.clear() ; 
  }
  void apply_propagation()
  {
    for (int i=0 ; i<size() ; i++)
    {
      rays[i].apply_propagation() ;
      rays[i].apply_polariser(Parameters.post_polarisation) ; 
      set_pixel(i, rays[i].get_intensity()) ;
      set_extra_value(i, rays[i].extra_value) ; 
    }
  }
  
  void process_rays(FEsolver &FE, std::vector<Grains> &grains)
  {
    switch (Parameters.strategy){
      case Strategy::LINEAR_NEARESTNEIGHBOUR :
        process_rays_linnearest(FE, grains) ;         
        break ; 
      case Strategy::LINEAR_TETRAHEDRON_INVERSION:
      case Strategy::LINEAR_TETRAHEDRON_MOLLERTRUMBORE:
        process_rays_tetralin(FE, grains) ; 
        break ; 
      case Strategy::TETRAHEDRON_EXPONENTIAL_INVERSION:
      case Strategy::TETRAHEDRON_EXPONENTIAL_MOLLERTRUMBORE:
        process_rays_tetraexp(FE, grains) ; 
        break ; 
      default:
        printf("ERR: unknown ray processing strategy.\n") ; 
    }    
  }
  
  void process_rays_linnearest(FEsolver &FE, std::vector<Grains> &grains)
  {
    std::vector<int> ids ; ids.resize(Parameters.Ns) ; 
    for (auto & ray : rays)
    {
      auto R = ray.get_rotation_matrix() ; 
      auto Rinv = ray.get_inv_rotation_matrix() ; 
    
      auto entryexits = ray.find_grains_entryexit(grains) ; 
      Geometry::reverse_sort_from(entryexits, ray.destination) ; 
    
      for (auto entry : entryexits)
      {
        auto [locations, ds] = entry.locations_inbetween(Parameters.Ns) ; 
        ids = FE.interpolate(locations) ; 
        auto stress = grains[entry.objectid].stress_at(ids) ; 
        
        for (auto & s: stress)
        {
          s = Rinv*s;
          s = s*R ; 
        }     
        ray.propagate(stress, Parameters.photoelastic_constant, ds) ;
      }
    }
  }
  
  void process_rays_tetralin(FEsolver &FE, std::vector<Grains> &grains)
  {
    std::vector<std::tuple<double, double, int>> tetra_intersections ; 
    std::vector<int> ids ; ids.resize(Parameters.Ns) ; 
    for (auto & ray : rays)
    {
      auto R = ray.get_rotation_matrix() ; 
      auto Rinv = ray.get_inv_rotation_matrix() ; 
    
      auto entryexits = ray.find_grains_entryexit(grains) ; 
      Geometry::reverse_sort_from(entryexits, ray.destination) ; 
    
      for (auto entry : entryexits)
      {
        auto [locations, ds] = entry.locations_inbetween(Parameters.Ns) ;  
  
        Geometry::intersection_ray_mesh (tetra_intersections, FE.tetras, {ray.destination[0], ray.destination[1], ray.destination[2]} , {ray.direction[0], ray.direction[1], ray.direction[2]}) ;  
      
        if (tetra_intersections.size() == 0) {printf("There should be an intersection ...\n") ; continue ; }
        
        size_t curid = 0 ; 
        for (size_t j=0 ; j<ids.size() ; j++)
        {
          double alpha = (locations[j][0]-ray.destination[0])/(normal[0]) ; //WARNING INCORRECT
          
          for ( ; alpha > std::get<1>(tetra_intersections[curid]) && curid<tetra_intersections.size()-1 ; curid++) ;
          ids[j] = std::get<2>(tetra_intersections[curid]) ; 
        }
        
        auto stress = grains[entry.objectid].stress_at(ids) ; 
        
        for (auto & s: stress)
        {
          s = Rinv*s;
          s = s*R ; 
        }        
        ray.propagate(stress, Parameters.photoelastic_constant, ds) ;
      }
    }
  }
  
  void process_rays_tetraexp(FEsolver &FE, std::vector<Grains> &grains)
  {
    std::vector<std::tuple<double, double, int>> tetra_intersections ; 
    std::vector<double> lengths ;
    std::vector<int> ids ;
    for (auto & ray : rays)
    {
      auto R = ray.get_rotation_matrix() ; 
      auto Rinv = ray.get_inv_rotation_matrix() ; 
    
      auto entryexits = ray.find_grains_entryexit(grains) ; 
      Geometry::reverse_sort_from(entryexits, ray.destination) ; 
    
      for (auto entry : entryexits)
      {
        auto [locations, ds] = entry.locations_inbetween(Parameters.Ns) ;     

        Geometry::intersection_ray_mesh (tetra_intersections, FE.tetras, {ray.destination[0], ray.destination[1], ray.destination[2]} , {ray.direction[0], ray.direction[1], ray.direction[2]}) ;  
      
        if (tetra_intersections.size() == 0) {printf("There should be an intersection ...\n") ; continue ; }
        
        lengths.resize(tetra_intersections.size()) ;
        ids.resize(tetra_intersections.size()) ;
        for (size_t j=0 ; j<tetra_intersections.size() ; j++)
        {
          lengths[j] = fabs(std::get<1>(tetra_intersections[j]) - std::get<0>(tetra_intersections[j])) ; 
          ids[j]=std::get<2>(tetra_intersections[j]) ; 
        }
        
        auto stress = grains[entry.objectid].stress_at(ids) ; 
        
        for (auto & s: stress)
        {
          s = Rinv*s;
          s = s*R ; 
        }        
        ray.propagate_exp(stress, Parameters.photoelastic_constant, lengths) ; 
      }
    }
  }
  
  void display(SDL_Renderer** renderer, SDL_Texture** texture) ;  
  void display_extra(SDL_Renderer** renderer, SDL_Texture** texture) ;
} ; 

#endif
