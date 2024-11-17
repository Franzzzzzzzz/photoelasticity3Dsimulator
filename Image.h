#ifndef IMAGE_H
#define IMAGE_H

#include <vector>
#include <SDL2/SDL.h>

#include "Geometry.h"

class Image {
public: 
  Image (int W,int H) { I.resize(H, std::vector<double>(W,0)) ; extra = I ; w=W ; h=H ; }
  std::vector<std::vector<double>> I ; 
  std::vector<std::vector<double>> extra ; 
  double dx, dz ;
  int w, h ; 
  vec origin ; 
  vec normal ;
  
  void set_dxdz(double ddx, double ddz) {dx=ddx ; dz=ddz ; }
  void set_origin (vec v) {origin = v ; }
  void set_normal (vec v) {normal = v ; }
  int size () {return w*h ;}
  
  vec get_location(int i) ;
  void set_pixel(int i, int j, double intensity) {I[i][j]=intensity ; }
  void set_pixel(int i, double intensity) {set_pixel(i/w, i%w, intensity) ; }
  void set_extra_value(int i, int j, double v) {extra[i][j]=v ;  }
  void set_extra_value(int i, double v) {set_extra_value(i/w, i%w, v) ; }
  
  void display(SDL_Renderer** renderer, SDL_Texture** texture) ;  
  void display_extra(SDL_Renderer** renderer, SDL_Texture** texture) ;
} ; 

#endif
