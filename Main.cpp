#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <thread>
#include <span>
#include <SDL2/SDL.h>

#include "Image.h"
#include "Grains.h"
#include "Ray.h"
#include "Fem.h"

struct Parameters_
{
  int imwidth=100, imheight=150 ;
  double im_dh=0.02, im_dv=0.02 ; 
  int screen_width=400, screen_height=600 ;
  int Ns = 100 ; 
  vec_jones polarisation = {1/sqrt(2),1/sqrt(2)} ; 
  double photoelastic_constant = 100; 
  double absorption = 1. ; 
  mat_jones post_polarisation = Polariser::vert_jones ; 
  
  std::string meshfile = "sphere_mesh.xdmf" ; 
  
  SDL_Window* window ;
  SDL_Renderer* renderer ;
  SDL_Texture* texture ;
} Parameters ; 

//----------------------------
int init_display()
{            
  int SCREEN_HEIGHT = Parameters.screen_height ; 
  int SCREEN_WIDTH = Parameters.screen_width ; 
  
  Parameters.window = SDL_CreateWindow("Display", 0, 0, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN|SDL_WINDOW_RESIZABLE);
  Parameters.renderer = SDL_CreateRenderer(Parameters.window, -1, SDL_RENDERER_ACCELERATED);
  Parameters.texture = SDL_CreateTexture(Parameters.renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, Parameters.imwidth, Parameters.imheight);
  
  if (!Parameters.window || !Parameters.renderer || !Parameters.texture) {
      std::cerr << "SDL window, renderer, or texture creation failed: " << SDL_GetError() << std::endl;
      return 0;
  }
  return 1 ; 
}
//==============================================================================
int main (int argc, char * argv[])
{
  init_display() ; 
  //ASSUMPTION: ray always perpendicular to vertical z axis. 
  
  Image image(Parameters.imwidth,Parameters.imheight) ; 
  std::vector<Grains> grains ; 
  
  image.set_normal({ 1,0,0}) ; 
  image.set_origin({-1,0,0}) ; 
  image.set_dxdz(Parameters.im_dh, Parameters.im_dv) ; 
  //set_grains_locations ();
  //set_grains_contacts () ;
  grains.push_back(Grains({0,0,0}, 0.5)) ; 
  
  FEsolver FE(argc, &argv) ; 
  
  auto start = std::chrono::high_resolution_clock::now();
  auto elapsed1 = std::chrono::high_resolution_clock::now()-start;
  auto elapsed2 = std::chrono::high_resolution_clock::now()-start;
  
  for (size_t i=0 ; i<grains.size() ; i++)
  {
    FE.prepare_mesh(Parameters.meshfile, grains[i].r) ;     
    FE.get_sigma(grains[i].stress) ;  
    printf("FEM finished\n") ; 
  }
  
  std::vector<std::tuple<double, double, int>> tetra_intersections ; 
  std::vector<int> ids ; ids.resize(Parameters.Ns) ; 
  std::vector<int> ids2 ; ids2.resize(Parameters.Ns) ; 
  
  
  for (int i=0 ; i<image.size() ; i++)
  {
    Ray ray ; 
    ray.set_direction(image.normal) ; 
    ray.set_destination(image.get_location(i)) ; 
    ray.set_polarisation (Parameters.polarisation) ; 
    
    auto R = ray.get_rotation_matrix() ; 
    auto Rinv = ray.get_inv_rotation_matrix() ; 
    /*for (int k=0 ; k<9 ; k++)
      printf("%g ", R[k]) ; 
    printf("\n") ; */
    
    auto entryexits = ray.find_grains_entryexit(grains) ; 
    Geometry::reverse_sort_from(entryexits, ray.destination) ; 
    
    for (auto entry : entryexits)
    {
      auto [locations, ds] = entry.locations_inbetween(Parameters.Ns) ; 
      auto ids = FE.interpolate(locations) ;       
      
      // Tetra intersection version
      Geometry::intersection_ray_mesh (tetra_intersections, FE.tetras, {ray.destination[0], ray.destination[1], ray.destination[2]} , {ray.direction[0], ray.direction[1], ray.direction[2]}) ;  
    
      /*if (tetra_intersections.size()<4)
      {
        for (auto & v: tetra_intersections)
          printf("%g %g %d\n", std::get<0>(v), std::get<1>(v), std::get<2>(v));
        printf("----------------\n") ; 
      }*/
      
      //std::vector<double> lengths ;  lengths.resize(tetra_intersections.size()) ;       
      size_t curid=0 ; 
      if (tetra_intersections.size() == 0) {printf("There should be an intersection ...\n") ; continue ; }
      
      for (size_t j=0 ; j<ids.size() ; j++)
      {
        double alpha = (locations[j][0]-image.get_location(i)[0])/(image.normal[0]) ; 
        
        for ( ; alpha > std::get<1>(tetra_intersections[curid]) && curid<tetra_intersections.size()-1 ; curid++) ;
        ids2[j] = std::get<2>(tetra_intersections[curid]) ; 
      }
      
      auto stress = grains[entry.objectid].stress_at(ids2) ; 
      
      for (auto & s: stress)
      {
        s = Rinv*s;
        s = s*R ; 
      }
      
      //ray.absorbe(Parameters.absorption, ds) ; 
      //start = std::chrono::high_resolution_clock::now();
      ray.propagate(stress, Parameters.photoelastic_constant, ds) ;
      //ray.propagate_exp(stress, Parameters.photoelastic_constant, lengths) ; 
      
      //elapsed2 += std::chrono::high_resolution_clock::now()-start;
           
      /*ray.extra_value = stress[0].norm() ; 
      for (size_t i=1 ; i<stress.size() ; i++)
        if (stress[i].norm()>ray.extra_value)
          ray.extra_value = stress[i].norm() ;*/
    }
    ray.apply_polariser (Parameters.post_polarisation) ; 
    image.set_pixel(i, ray.get_intensity()) ;
    image.set_extra_value(i, ray.extra_value) ; 
    //auto duration1= std::chrono::duration_cast<std::chrono::microseconds>(elapsed1).count() ; 
    //printf("%g\n",duration1/1000000.) ; 
  }  
  
  printf("ImageDisplayed") ; 
  //auto duration1= std::chrono::duration_cast<std::chrono::microseconds>(elapsed1).count() ; 
  //auto duration2= std::chrono::duration_cast<std::chrono::microseconds>(elapsed2).count() ; 
  //printf("%g %g\n", duration1/1000000., duration2/1000000.) ; 
  image.display(&Parameters.renderer, &Parameters.texture) ; 
  //--------------------------------------------
  SDL_Event event;
  while ( SDL_WaitEvent(&event) >= 0 ) {
        switch (event.type) {
          case SDL_WINDOWEVENT:
            image.display(&Parameters.renderer, &Parameters.texture) ;   
            break ;
          case SDL_MOUSEBUTTONDOWN:  
           {
            double xscale=Parameters.screen_width/(double)Parameters.imwidth ; 
            double yscale=Parameters.screen_height/(double)Parameters.imheight ; 
            printf("%g %g\n", event.button.y/yscale, event.button.x/xscale) ; 
           }
           break ; 
          case SDL_QUIT:
            exit(0);
          break;
        }
    }
  
}
