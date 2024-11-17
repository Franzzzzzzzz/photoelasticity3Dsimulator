#include "Image.h"

vec Image::get_location(int i) 
{
  vec pt ; 
  pt[0] = origin[0] - (i%w-w/2)*dz * normal[1] ;
  pt[1] = origin[1] + (i%w-w/2)*dz * normal[0] ;
  pt[2] = origin[2] + (i/w-h/2)*dz ;
  return pt ; 
}

void Image::display(SDL_Renderer** renderer, SDL_Texture** texture) 
{
  
  
  SDL_Color* pixels = nullptr;
  int pitch;
  SDL_LockTexture(*texture, nullptr, (void**)&pixels, &pitch);
  SDL_Color color ; 
  
  // Scale 
  double min=I[0][0], max=I[0][0] ; 
  for (int i =0 ; i<h ; i++)
    for (int j=0 ; j<w ; j++)
    {
      if (I[i][j]<min) min=I[i][j] ; 
      if (I[i][j]>max) max=I[i][j] ; 
    }
  min=0 ; max=1 ; 

  for (int i =0 ; i<h ; i++)
    for (int j=0 ; j<w ; j++)
    {
      uint8_t c= (I[i][j]-min)/(max-min)*255 ; 
      color = {0,c,c,c} ;         
      pixels[ i * (pitch / sizeof(Uint32)) + j] = color;
      
    }

  SDL_UnlockTexture(*texture);

  // Clear the renderer and render the texture
  SDL_RenderClear(*renderer);
  SDL_RenderCopy(*renderer, *texture, nullptr, nullptr);
  SDL_RenderPresent(*renderer);
} 

void Image::display_extra(SDL_Renderer** renderer, SDL_Texture** texture) 
{
  
  
  SDL_Color* pixels = nullptr;
  int pitch;
  SDL_LockTexture(*texture, nullptr, (void**)&pixels, &pitch);
  SDL_Color color ; 
  
  // Scale 
  double min=extra[0][0], max=extra[0][0] ; 
  for (int i =0 ; i<h ; i++)
    for (int j=0 ; j<w ; j++)
    {
      if (extra[i][j]<min) min=extra[i][j] ; 
      if (extra[i][j]>max) max=extra[i][j] ; 
    }
    
  printf("%g %g\n", min, max) ;
  min=0 ; max=0.4 ; 

  for (int i =0 ; i<h ; i++)
    for (int j=0 ; j<w ; j++)
    {
      uint8_t c= (extra[i][j]-min)/(max-min)*255 ; 
      color = {0,c,c,c} ;         
      pixels[ i * (pitch / sizeof(Uint32)) + j] = color;
      
    }

  SDL_UnlockTexture(*texture);

  // Clear the renderer and render the texture
  SDL_RenderClear(*renderer);
  SDL_RenderCopy(*renderer, *texture, nullptr, nullptr);
  SDL_RenderPresent(*renderer);
}  
