#ifndef PARAMETERS_
#define PARAMETERS_

#include "Algebra.h"
#include <SDL2/SDL.h>

/* Strategies:
 * LINEAR_NEARESTNEIGHBOUR: equidistant sampling along the ray within the grain, each point is allocated the stress in the "nearest" tetrahedron, calculated from the center of mass of the tetrahedron.
 * LINEAR_TETRAHEDRON_INVERSION: equidistant sampling along the ray within the grain, however the intersection with each tetrahedra within the ray are actually calculated, and the stress used at each point is the one from the tetrahedron it belongs to.
 * LINEAR_TETRAHEDRON_MOLLER–TRUMBORE: same as previous, but the ray-tetrahedron intersection is calculated using the Möller Trumbore algorithm, rather than the matrix inversion (my own algo). 
 * TETRAHEDRON_EXPONENTIAL_INVERSION: the intersection with tetrahedra is calculated, and the polarisation transformation matrix is calculated from the actual matrix exponentiation, rather than repeated exponentiation as with the previous methods.
 * TETRAHEDRON_EXPONENTIAL_MOLLER–TRUMBORE: same as previous, but the ray-tetrahedron intersection is calculated using the Möller Trumbore algorithm, rather than the matrix inversion (my own algo). 
 */

enum Strategy {LINEAR_NEARESTNEIGHBOUR=0, LINEAR_TETRAHEDRON_INVERSION=1, LINEAR_TETRAHEDRON_MOLLERTRUMBORE=2, TETRAHEDRON_EXPONENTIAL_INVERSION=3, TETRAHEDRON_EXPONENTIAL_MOLLERTRUMBORE=4} ; 

struct Parameters_
{
  int imwidth=100, imheight=150 ;
  double im_dh=0.02, im_dv=0.02 ; 
  int screen_width=400, screen_height=600 ;
  int Ns = 100 ; 
  vec_jones polarisation = Polarisation::Lp45_jones ; 
  double photoelastic_constant = 100; 
  double absorption = 1. ; 
  mat_jones post_polarisation = Polariser::vert_jones ; 
  Strategy strategy = LINEAR_TETRAHEDRON_MOLLERTRUMBORE; 
  
  std::string meshfile = "sphere_mesh.xdmf" ; 
  
  SDL_Window* window ;
  SDL_Renderer* renderer ;
  SDL_Texture* texture ;
  
  void disp_pol() {
    printf(" -> ") ; 
    if (polarisation == Polarisation::Lp45_jones)
      printf("/") ; 
    else if (polarisation == Polarisation::Lm45_jones)
      printf("\\") ; 
    else if (polarisation == Polarisation::horiz_jones)
      printf("-") ; 
    else if (polarisation == Polarisation::vert_jones)
      printf("|") ; 
    else if (polarisation == Polarisation::circleft_jones)
      printf("\u21bb") ; 
    else if (polarisation == Polarisation::circright_jones)
      printf("\u21ba") ;
    printf(" ooooo ") ; 
    
    if (post_polarisation == Polariser::Lp45_jones)
      printf("/") ; 
    else if (post_polarisation == Polariser::Lm45_jones)
      printf("\\") ; 
    else if (post_polarisation == Polariser::horiz_jones)
      printf("-") ; 
    else if (post_polarisation == Polariser::vert_jones)
      printf("|") ; 
    else if (post_polarisation == Polariser::circleft_jones)
      printf("\u21bb") ; 
    else if (post_polarisation == Polariser::circright_jones)
      printf("\u21ba") ;
    printf("\n") ; 
  } ; 
  
} ; 

#endif
