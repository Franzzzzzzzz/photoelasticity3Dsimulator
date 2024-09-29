#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <thread>
#include <SDL2/SDL.h>

#include "Image.h"
#include "Grains.h"
#include "Ray.h"

struct Parameters_
{
  int imwidth=200, imheight=300 ;
  int Ns = 10 ; 
  vec_stokes polarisation{1,0,0,0} ; 
  double photoelastic_constant = 0.1 ; 
  
  SDL_Window* window ;
  SDL_Renderer* renderer ;
  SDL_Texture* texture ;
} Parameters ; 

//----------------------------
int init_display()
{            
  int SCREEN_HEIGHT = Parameters.imheight ; 
  int SCREEN_WIDTH = Parameters.imwidth ; 
  
  Parameters.window = SDL_CreateWindow("Display", 0, 0, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN);
  Parameters.renderer = SDL_CreateRenderer(Parameters.window, -1, SDL_RENDERER_ACCELERATED);
  Parameters.texture = SDL_CreateTexture(Parameters.renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, SCREEN_WIDTH, SCREEN_HEIGHT);
  
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
             sqrt( (sigma[i][0]-sigma[i][2])*(sigma[i][0]-sigma[i][2]) + 4 * sigma[i][1] * sigma[i][1] ) * 
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
    phis[i] = 0.5*atan(2*sigma[i][1]/(sigma[i][0]-sigma[i][2])) ;
  for (size_t i=1 ; i<phis.size()-1 ; i++)
  {
    dphis[i] = phis[i+1]+phis[i-1] ; 
    if (dphis[i]>=M_PI/2) dphis[i]-=M_PI ; 
    if (dphis[i]<=M_PI/2) dphis[i]+=M_PI ; 
    dphis[i] /= (2.*ds) ;
  }
  dphis[0] = (dphis[1]-dphis[1])/ds ; 
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
  
  
  
  for (int i=0 ; i<image.size() ; i++)
  {
    Ray ray ; 
    ray.set_direction(-image.normal) ; 
    ray.set_destination(image.get_location(i)) ; 
    ray.set_polarisation (Parameters.polarisation) ; 
    
    auto R = ray.get_rotation_matrix() ; 
    auto Rinv = ray.get_inv_rotation_matrix() ; 
    
    auto entryexits = ray.find_grains_entryexit(grains) ; 
    Geometry::reverse_sort_from(entryexits, ray.destination) ; 
    
    for (auto entry : entryexits)
    {
      auto [locations, ds] = entry.locations_inbetween(Parameters.Ns) ; 
      auto stress = grains[entry.objectid].stress_at(locations) ; 
      
      for (auto & s: stress)
      {
        s = Rinv*s;
        s = s*R ; 
      }
      auto ddelta = get_photoelastic_deltas(stress, ds) ; 
      auto dphi = get_photoelastic_dphis(stress, ds) ;
      
      ray.propagate(ddelta, dphi) ;
    }
    
    //ray.apply_polarisor (Parameters.post_polarisation) ; 
    image.set_pixel(i, ray.get_intensity()) ;     
  }  
  
  image.display(&Parameters.renderer, &Parameters.texture) ;   
  //--------------------------------------------
  SDL_Event event;
  while ( SDL_WaitEvent(&event) >= 0 ) {
        switch (event.type) {            
            case SDL_QUIT: {
                exit(0);
            }
            break;
        }
    }
  
}
