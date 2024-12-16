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

struct Parameters_ Parameters ; 

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
  if (argc>1)
  {
    auto j = Parameters.parameters_from_file(argv[1]) ; 
  }
  
  init_display() ; 
  //ASSUMPTION: ray always perpendicular to vertical z axis. 
  
  Image image(Parameters.imwidth,Parameters.imheight) ; 
   
  image.set_normal_and_origin(Parameters.distance, Parameters.azimuth) ; 
  image.set_dxdz(Parameters.im_dh, Parameters.im_dv) ; 
  
  FEsolver FE(argc, &argv) ; 
  
  auto start = std::chrono::high_resolution_clock::now();
  auto elapsed1 = std::chrono::high_resolution_clock::now()-start;
  auto elapsed2 = std::chrono::high_resolution_clock::now()-start;
  
  for (size_t i=0 ; i<Parameters.grains.size() ; i++)
  {
    FE.prepare_mesh(Parameters.meshfile, Parameters.grains[i].r) ;     
    switch (Parameters.grains[i].contactpoints.size()) {
      case 1: FE.get_sigma<1>(Parameters.grains[i].stress, Parameters.grains[i].contactpoints, Parameters.grains[i].displacements) ; break ;
      case 2: FE.get_sigma<2>(Parameters.grains[i].stress, Parameters.grains[i].contactpoints, Parameters.grains[i].displacements) ; break ;
      case 3: FE.get_sigma<3>(Parameters.grains[i].stress, Parameters.grains[i].contactpoints, Parameters.grains[i].displacements) ; break ;
      case 4: FE.get_sigma<4>(Parameters.grains[i].stress, Parameters.grains[i].contactpoints, Parameters.grains[i].displacements) ; break ;
      case 5: FE.get_sigma<5>(Parameters.grains[i].stress, Parameters.grains[i].contactpoints, Parameters.grains[i].displacements) ; break ;
      case 6: FE.get_sigma<6>(Parameters.grains[i].stress, Parameters.grains[i].contactpoints, Parameters.grains[i].displacements) ; break ;
      default :
        printf("ERR: this number contacts is not implemented.\n") ; 
        FE.get_sigma<0>(Parameters.grains[i].stress, Parameters.grains[i].contactpoints, Parameters.grains[i].displacements) ; break ;
    }
    printf("FEM finished\n") ; 
  }
  
  
  image.set_rays() ; 
  
  //for (int k=0 ; k<24; k++)
  //{
  image.reset_rays() ; 
  image.set_normal_and_origin(Parameters.distance, Parameters.azimuth) ; 
  image.set_rays() ; 
  
  image.process_rays(FE, Parameters.grains) ; 
  image.apply_propagation() ; 
  image.display(&Parameters.renderer, &Parameters.texture) ; 
  
  printf("ImageDisplayed") ; 
  
  /*SDL_Surface* pScreenShot = SDL_CreateRGBSurface(0, Parameters.screen_width, Parameters.screen_height, 32, 0x00ff0000, 0x0000ff00, 0x000000ff, 0xff000000);
  if (pScreenShot)
  {
    SDL_RenderReadPixels(Parameters.renderer, NULL, SDL_GetWindowPixelFormat(Parameters.window), pScreenShot->pixels, pScreenShot->pitch);
    // Create the bmp screenshot file
    std::string filename = "Screenshot-" + std::to_string(int(k)) + ".bmp" ; 
    SDL_SaveBMP(pScreenShot, filename.c_str());
    SDL_FreeSurface(pScreenShot);
  }
  
  Parameters.azimuth -= 15./2. * M_PI/180. ;
  printf("Current azimuth: %g \n", Parameters.azimuth/M_PI*180.) ;
  }*/
  
  //auto duration1= std::chrono::duration_cast<std::chrono::microseconds>(elapsed1).count() ; 
  //auto duration2= std::chrono::duration_cast<std::chrono::microseconds>(elapsed2).count() ; 
  //printf("%g %g\n", duration1/1000000., duration2/1000000.) ; 
  //image.display(&Parameters.renderer, &Parameters.texture) ; 
  //--------------------------------------------
  SDL_Event event;
  while ( SDL_WaitEvent(&event) >= 0 ) {
        switch (event.type) {
          case SDL_WINDOWEVENT:
            //image.display(&Parameters.renderer, &Parameters.texture) ;   
            break ;
          case SDL_KEYDOWN:
            switch(event.key.keysym.sym) {
              case SDLK_r: 
                if (event.key.keysym.mod & KMOD_SHIFT)
                {
                  Parameters.polarisation = Polarisation::circright_jones ; 
                  image.set_rays() ;                  
                } 
                else
                  Parameters.post_polarisation = Polariser::circright_jones ; 
                image.apply_propagation() ; 
                Parameters.disp_pol() ; 
                break ; 
              case SDLK_l:
                if (event.key.keysym.mod & KMOD_SHIFT)
                {
                  Parameters.polarisation = Polarisation::circleft_jones ; 
                  image.set_rays() ;                  
                } 
                else
                  Parameters.post_polarisation = Polariser::circleft_jones ; 
                image.apply_propagation() ; 
                Parameters.disp_pol() ; 
                break ; 
              case SDLK_SLASH:
              case SDLK_KP_DIVIDE:
                if (event.key.keysym.mod & KMOD_SHIFT)
                {
                  Parameters.polarisation = Polarisation::Lp45_jones ; 
                  image.set_rays() ;                  
                } 
                else
                  Parameters.post_polarisation = Polariser::Lp45_jones ; 
                image.apply_propagation() ; 
                Parameters.disp_pol() ; 
                break ; 
              case SDLK_BACKSLASH:
                if (event.key.keysym.mod & KMOD_SHIFT)
                {
                  Parameters.polarisation = Polarisation::Lm45_jones ; 
                  image.set_rays() ;                  
                } 
                else
                  Parameters.post_polarisation = Polariser::Lm45_jones ; 
                image.apply_propagation() ; 
                Parameters.disp_pol() ; 
                break ; 
              case SDLK_MINUS:
              case SDLK_KP_MINUS:
                if (event.key.keysym.mod & KMOD_SHIFT)
                {
                  Parameters.polarisation = Polarisation::horiz_jones ; 
                  image.set_rays() ;                  
                } 
                else
                  Parameters.post_polarisation = Polariser::horiz_jones ; 
                image.apply_propagation() ; 
                Parameters.disp_pol() ; 
                break ; 
              case SDLK_EXCLAIM:
              case SDLK_KP_VERTICALBAR:
                if (event.key.keysym.mod & KMOD_SHIFT)
                {
                  Parameters.polarisation = Polarisation::vert_jones ; 
                  image.set_rays() ;                  
                } 
                else
                  Parameters.post_polarisation = Polariser::vert_jones ; 
                image.apply_propagation() ; 
                Parameters.disp_pol() ; 
                break ; 
              case SDLK_LEFT:
                Parameters.azimuth += 15. * M_PI/180. ;
                printf("%g ", Parameters.azimuth/M_PI*180.) ; 
                image.reset_rays() ; 
                image.set_normal_and_origin(Parameters.distance, Parameters.azimuth) ; 
                image.set_rays() ; 
                image.process_rays(FE, Parameters.grains) ; 
                image.apply_propagation() ; 
                image.display(&Parameters.renderer, &Parameters.texture) ; 
                break ; 
              case SDLK_RIGHT:
                Parameters.azimuth -= 15. * M_PI/180. ;
                printf("%g ", Parameters.azimuth/M_PI*180.) ;
                image.reset_rays() ; 
                image.set_normal_and_origin(Parameters.distance, Parameters.azimuth) ; 
                image.set_rays() ; 
                image.process_rays(FE, Parameters.grains) ; 
                image.apply_propagation() ; 
                image.display(&Parameters.renderer, &Parameters.texture) ; 
                break ;
              case SDLK_c:
                image.reset_rays_rgb() ; 
                image.set_normal_and_origin(Parameters.distance, Parameters.azimuth) ; 
                image.set_rays_rgb() ; 
                image.process_rays_rgb(FE, Parameters.grains) ; 
                image.apply_propagation_rgb() ; 
                image.display_rgb(&Parameters.renderer, &Parameters.texture) ; 
                break ;
              case SDLK_s:
                {
                  // Read the pixels from the current render target and save them onto the surface
                  SDL_Surface* pScreenShot = SDL_CreateRGBSurface(0, Parameters.screen_width, Parameters.screen_height, 32, 0x00ff0000, 0x0000ff00, 0x000000ff, 0xff000000);
                  if (pScreenShot)
                  {
                    SDL_RenderReadPixels(Parameters.renderer, NULL, SDL_GetWindowPixelFormat(Parameters.window), pScreenShot->pixels, pScreenShot->pitch);
                    // Create the bmp screenshot file
                    SDL_SaveBMP(pScreenShot, "Screenshot.bmp");
                    SDL_FreeSurface(pScreenShot);
                  }
                }
                break ; 
            }
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
