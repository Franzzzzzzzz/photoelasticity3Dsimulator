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
  std::vector<std::vector<double>> I, I_r, I_g, I_b ; 
  std::vector<std::vector<double>> extra ; 
  std::vector<Ray> rays_white ; 
  std::vector<Ray> rays_red, rays_green, rays_blue ;    
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
  void set_pixel_rgb(int i, int j, double intensity, char color) 
  {
    if (color=='r') I_r[i][j]=intensity ; 
    if (color=='g') I_g[i][j]=intensity ; 
    if (color=='b') I_b[i][j]=intensity ; 
  }
  void set_pixel(int i, int j, double intensity) {I[i][j]=intensity ; }
  void set_pixel(int i, double intensity) {set_pixel(i/w, i%w, intensity) ; }
  void set_pixel_rgb(int i, double intensity, char color) {set_pixel_rgb(i/w, i%w, intensity, color) ; }
  void set_extra_value(int i, int j, double v) {extra[i][j]=v ;  }
  void set_extra_value(int i, double v) {set_extra_value(i/w, i%w, v) ; }
  
  void set_rays_rgb () 
  {
    set_rays(rays_red) ;
    set_rays(rays_green) ;
    set_rays(rays_blue) ;
  }
  void set_rays () {return set_rays(rays_white) ;}
  void set_rays (std::vector<Ray> & rays) 
  {
    rays.resize(size()) ; 
    for (int i=0 ; i<size() ; i++)
    {
      rays[i].set_direction(normal) ; 
      rays[i].set_destination(get_location(i)) ; 
      rays[i].set_polarisation (Parameters.polarisation) ; 
    }
  }
  void reset_rays_rgb() 
  {
    reset_rays(rays_red) ;
    reset_rays(rays_green) ;
    reset_rays(rays_blue) ;
  }
  void reset_rays() {return reset_rays(rays_white) ;}
  void reset_rays(std::vector<Ray> & rays)
  {
    rays.clear() ; 
  }
  
  void apply_propagation_rgb() 
  {
    I_r.resize(h, std::vector<double>(w,0)) ;
    I_g.resize(h, std::vector<double>(w,0)) ;
    I_b.resize(h, std::vector<double>(w,0)) ;
    for (int i=0 ; i<size() ; i++)
    {
      rays_red[i].apply_propagation() ;
      rays_red[i].apply_polariser(Parameters.post_polarisation) ; 
      rays_green[i].apply_propagation() ;
      rays_green[i].apply_polariser(Parameters.post_polarisation) ; 
      rays_blue[i].apply_propagation() ;
      rays_blue[i].apply_polariser(Parameters.post_polarisation) ; 
      set_pixel_rgb(i, rays_red[i].get_intensity(), 'r') ;
      set_pixel_rgb(i, rays_green[i].get_intensity(), 'g') ;
      set_pixel_rgb(i, rays_blue[i].get_intensity(), 'b') ;
    }
  }
  void apply_propagation()
  {
    for (int i=0 ; i<size() ; i++)
    {
      rays_white[i].apply_propagation() ;
      rays_white[i].apply_polariser(Parameters.post_polarisation) ; 
      set_pixel(i, rays_white[i].get_intensity()) ;
      set_extra_value(i, rays_white[i].extra_value) ; 
    }
  }
    
  void process_rays_rgb(FEsolver &FE, std::vector<Grains> &grains)
  {
    process_rays(FE, grains, rays_red, Parameters.photoelastic_constant_red) ; 
    process_rays(FE, grains, rays_green, Parameters.photoelastic_constant_green) ; 
    process_rays(FE, grains, rays_blue, Parameters.photoelastic_constant_blue) ; 
  }
  void process_rays(FEsolver &FE, std::vector<Grains> &grains)
  {
    return process_rays(FE, grains, rays_white, Parameters.photoelastic_constant) ; 
  }
  void process_rays(FEsolver &FE, std::vector<Grains> &grains,std::vector<Ray> & rays, double photoelastic_constant = Parameters.photoelastic_constant)
  {
    switch (Parameters.strategy){
      case Strategy::LINEAR_NEARESTNEIGHBOUR :
        process_rays_linnearest(FE, grains, rays, photoelastic_constant) ;         
        break ; 
      case Strategy::LINEAR_TETRAHEDRON_INVERSION:
      case Strategy::LINEAR_TETRAHEDRON_MOLLERTRUMBORE:
        process_rays_tetralin(FE, grains, rays, photoelastic_constant) ; 
        break ; 
      case Strategy::TETRAHEDRON_EXPONENTIAL_INVERSION:
      case Strategy::TETRAHEDRON_EXPONENTIAL_MOLLERTRUMBORE:
        process_rays_tetraexp(FE, grains, rays, photoelastic_constant) ; 
        break ; 
      default:
        printf("ERR: unknown ray processing strategy.\n") ; 
    }    
  }
  
  void process_rays_linnearest(FEsolver &FE, std::vector<Grains> &grains, std::vector<Ray> & rays, double photoelastic_constant)
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
        ray.propagate(stress, photoelastic_constant, ds) ;
      }
    }
  }
  
  void process_rays_tetralin(FEsolver &FE, std::vector<Grains> &grains, std::vector<Ray> & rays, double photoelastic_constant)
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
        ray.propagate(stress, photoelastic_constant, ds) ;
      }
    }
  }
  
  void process_rays_tetraexp(FEsolver &FE, std::vector<Grains> &grains, std::vector<Ray> & rays, double photoelastic_constant)
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
        ray.propagate_exp(stress, photoelastic_constant, lengths) ; 
      }
    }
  }
  
  void display(SDL_Renderer** renderer, SDL_Texture** texture) ;  
  void display_rgb(SDL_Renderer** renderer, SDL_Texture** texture) ; 
  void display_extra(SDL_Renderer** renderer, SDL_Texture** texture) ;
} ; 

#endif
