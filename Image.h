#include <vector>
#include <SDL2/SDL.h>

#include "Geometry.h"







class Image {
public: 
  Image (int W,int H) { I.resize(H, std::vector<double>(W,0)) ; w=W ; h=H ; }
  std::vector<std::vector<double>> I ; 
  double dx, dz ;
  int w, h ; 
  vec origin ; 
  vec normal ;
  
  void set_dxdz(double ddx, double ddz) {dx=ddx ; dz=ddz ; }
  void set_origin (vec v) {origin = v ; }
  void set_normal (vec v) {normal = v ; }
  int size () {return w*h ;}
  
  vec get_location(int i) 
  { 
    vec pt ; 
    pt[0] = origin[2] - (i%w-w/2)*dz * normal[1] ;
    pt[1] = origin[1] + (i%w-w/2)*dz * normal[0] ;
    pt[2] = origin[2] + (i/w-h/2)*dz ;
    return pt ; 
  }
  void set_pixel(int i, int j, double intensity) {I[i][j]=intensity ; }
  void set_pixel(int i, double intensity) {set_pixel(i/w, i%w, intensity) ; }
  
  void display(SDL_Renderer** renderer, SDL_Texture** texture) 
  {
    
    
    SDL_Color* pixels = nullptr;
    int pitch;
    SDL_LockTexture(*texture, nullptr, (void**)&pixels, &pitch);
    SDL_Color color ; 
    
    // Scale 
    double min=1000000, max=0 ; 
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
} ; 

