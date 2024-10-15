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
  int imwidth=200, imheight=300 ;
  int screen_width=400, screen_height=600 ;
  int Ns = 100 ; 
  vec_stokes polarisation{1,0,0,1} ; 
  double photoelastic_constant = 1; 
  double absorption = 1. ; 
  mat_stokes post_polarisation = Polariser::deg45_lin ; 
  
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
//-----------------------------
std::vector<double> get_photoelastic_deltas(std::vector<tens> &sigma, double ds)
{
  std::vector<double> res ;
  res.resize(sigma.size()) ; 
  for (size_t i = 0 ; i<res.size(); i++)
  {
    res[i] = Parameters.photoelastic_constant *
             sqrt( (sigma[i][4]-sigma[i][8])*(sigma[i][4]-sigma[i][8]) + 4 * sigma[i][5] * sigma[i][5] ) * 
             ds ;
  }
  return res ; 
}
//-----------------------------
std::vector<double> get_photoelastic_dphis(std::vector<tens> &sigma, double ds)
{
  std::vector<double> phis, dphis ; 
  phis.resize(sigma.size()) ; 
  dphis.resize(sigma.size()) ; 
  for (size_t i=0 ; i<phis.size() ; i++)
  { 
    phis[i] = 0.5*atan(2*sigma[i][5]/(sigma[i][4]-sigma[i][8])) ;
  }
  for (size_t i=1 ; i<phis.size()-1 ; i++)
  {
    dphis[i] = phis[i+1]+phis[i-1] ; 
    if (dphis[i]>=M_PI/2) dphis[i]-=M_PI ; 
    if (dphis[i]<=-M_PI/2) dphis[i]+=M_PI ; 
    dphis[i] /= (2.*ds) ;
  }
  dphis[0] = (dphis[1]-dphis[0])/ds ; 
  dphis[dphis.size()-1] = (dphis[dphis.size()-1]-dphis[dphis.size()-2])/ds ; 
  return dphis ;  
}

//==============================================================================
int main (int argc, char * argv[])
{
  init_display() ; 
  //ASSUMPTION: ray always perpendicular to vertical z axis. 
  
  Image image(Parameters.imwidth,Parameters.imheight) ; 
  std::vector<Grains> grains ; 
  
  image.set_normal({1,0,0}) ; 
  image.set_origin({0,0,0}) ; 
  image.set_dxdz(0.01, 0.01) ; 
  //set_grains_locations ();
  //set_grains_contacts () ;
  grains.push_back(Grains({0,0,0}, 0.5)) ; 
  
  FEsolver FE(argc, &argv) ; 
  
  for (size_t i=0 ; i<grains.size() ; i++)
  {
    FE.prepare_mesh(Parameters.meshfile, grains[i].r) ; 
    FE.get_sigma(grains[i].stress) ;  
    printf("FEM finished\n") ; 
  }
  
  std::vector<vec> a = {{-0.47905, -0.0938368, 0.0236803}} ; 
  printf("%d \n", FE.interpolate(a)[0]) ; 
  
  for (int i=0 ; i<image.size() ; i++)
  {
    Ray ray ; 
    ray.set_direction(-image.normal) ; 
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
      auto stress = grains[entry.objectid].stress_at(ids) ; 
      //auto stress = grains[entry.objectid].stress_at(locations) ; 
      
      for (auto & s: stress)
      {
        s = Rinv*s;
        s = s*R ; 
      }
      // from P. S. Theocaris et al., Matrix Theory of Photoelasticity
      auto ddelta = get_photoelastic_deltas(stress, ds) ; 
      auto dphi = get_photoelastic_dphis(stress, ds) ;
      
      //ray.absorbe(Parameters.absorption, ds) ; 
      ray.propagate(ddelta, dphi) ;
      
      /*for (int k=0 ; k<9 ; k++)
        printf("%g ", stress[0][k]) ; 
      printf("\n") ; */
      
      ray.extra_value = stress[0].norm() ; 
      for (size_t i=1 ; i<stress.size() ; i++)
        if (stress[i].norm()>ray.extra_value)
          ray.extra_value = stress[i].norm() ;
    }
    
    ray.apply_polariser (Parameters.post_polarisation) ; 
    image.set_pixel(i, ray.get_intensity()) ;  
    image.set_extra_value(i, ray.extra_value) ; 
  }  
  
  printf("ImageDisplayed") ; 
  image.display(&Parameters.renderer, &Parameters.texture) ;   
  //--------------------------------------------
  SDL_Event event;
  while ( SDL_WaitEvent(&event) >= 0 ) {
        switch (event.type) {
          case SDL_WINDOWEVENT:
            image.display(&Parameters.renderer, &Parameters.texture) ;   
            break ;
          case SDL_QUIT:
            exit(0);
          break;
        }
    }
  
}
