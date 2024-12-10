#ifndef PARAMETERS_
#define PARAMETERS_

#include "Algebra.h"
#include <SDL2/SDL.h>
#include "json.hpp"
#include <fstream>
#include <regex>
#include "Grains.h"

/* Strategies:
 * LINEAR_NEARESTNEIGHBOUR: equidistant sampling along the ray within the grain, each point is allocated the stress in the "nearest" tetrahedron, calculated from the center of mass of the tetrahedron.
 * LINEAR_TETRAHEDRON_INVERSION: equidistant sampling along the ray within the grain, however the intersection with each tetrahedra within the ray are actually calculated, and the stress used at each point is the one from the tetrahedron it belongs to.
 * LINEAR_TETRAHEDRON_MOLLERTRUMBORE: same as previous, but the ray-tetrahedron intersection is calculated using the Möller Trumbore algorithm, rather than the matrix inversion (my own algo). 
 * TETRAHEDRON_EXPONENTIAL_INVERSION: the intersection with tetrahedra is calculated, and the polarisation transformation matrix is calculated from the actual matrix exponentiation, rather than repeated exponentiation as with the previous methods.
 * TETRAHEDRON_EXPONENTIAL_MOLLERTRUMBORE: same as previous, but the ray-tetrahedron intersection is calculated using the Möller Trumbore algorithm, rather than the matrix inversion (my own algo). 
 */

enum Strategy {LINEAR_NEARESTNEIGHBOUR=0, LINEAR_TETRAHEDRON_INVERSION=1, LINEAR_TETRAHEDRON_MOLLERTRUMBORE=2, TETRAHEDRON_EXPONENTIAL_INVERSION=3, TETRAHEDRON_EXPONENTIAL_MOLLERTRUMBORE=4} ; 

struct Parameters_
{  
  SDL_Window* window ;
  SDL_Renderer* renderer ;
  SDL_Texture* texture ;
   
  // Parameters for the user
  int imwidth=100, imheight=150 ;
  double im_dh=0.01, im_dv=0.01 ; 
  int screen_width=400, screen_height=600 ;
  int Ns = 100 ; 
  vec_jones polarisation = Polarisation::Lp45_jones ; 
  double photoelastic_constant = 100 ; 
  double photoelastic_constant_red   =  100; 
  double photoelastic_constant_green =  80; 
  double photoelastic_constant_blue  =  60; 
  bool color = false ; 
  double absorption = 0. ; 
  mat_jones post_polarisation = Polariser::vert_jones ; 
  Strategy strategy = LINEAR_TETRAHEDRON_MOLLERTRUMBORE; 
  
  std::vector<Grains> grains ; 
  std::string meshfile = "sphere_mesh_coarse.xdmf" ; 
  
  void parse_json_grains (const nlohmann::json & j)
  {
    for (auto & g: j)
    {
      grains.push_back( {vec(g["center"][0],g["center"][1],g["center"][2]), g["radius"].get<double>()} ) ;
      
      if (g.contains("contacts"))
      {
        for (auto &c: g["contacts"])
        {
          grains.back().contactpoints.push_back({c["location"][0], c["location"][1], c["location"][2]}) ; 
          grains.back().forces.push_back({c["force"][0], c["force"][1], c["force"][2]}) ; 
        }
      }
    }
  }
  
  void parse_json (const nlohmann::json &j)
  {
    for (auto& [key, value] : j.items())
    {
      if (key == "color") 
        color = value ; 
      else if (key == "sampling")
        Ns = value ; 
      else if (key == "absorption")
        absorption = value ; 
      else if (key == "mesh")
        meshfile = value ; 
      else if (key == "grains")
        parse_json_grains(j["grains"]) ; 
      else if (key=="image size") 
      {
        imwidth = value[0] ; 
        imheight = value[1] ;          
      }
      else if (key == "pixel size")
      {
        if (value.is_array())
        {
          im_dh = value[0] ; 
          im_dv = value[1] ; 
        }
        else 
          im_dh = im_dv = value ; 
      }
      else if (key== "display size") 
      {
        screen_width = value[0] ; 
        screen_height = value[1] ;          
      }
      else if (key == "photoelastic constant")
      {
        if (value.is_array())
        {
          int idx=0 ; 
          if (value.size()==4)
            photoelastic_constant = value[idx++] ; 
          photoelastic_constant_red = value[idx++] ; 
          photoelastic_constant_green = value[idx++] ; 
          photoelastic_constant_blue = value[idx++] ;
          color = true ; 
        }
        else 
        {
          photoelastic_constant = value ; 
          color = false ; 
        }
      }
      else if (key == "pre-polarisation")
      {
        if (value.is_number())
          polarisation = Polarisation::linear_jones(value) ; 
        else if (value == "left")
          polarisation = Polarisation::circleft_jones ; 
        else if (value == "right")
          polarisation = Polarisation::circright_jones ; 
        else if (value == "horizontal")
          polarisation = Polarisation::horiz_jones ; 
        else if (value == "vertical")
          polarisation = Polarisation::vert_jones ; 
        else 
          printf("Unknown polarisation %s\n", value.get<std::string>().c_str()) ; 
      }
      else if (key == "post_polarisation")
      {
        if (value.is_number())
          post_polarisation = Polariser::linear_jones(value) ; 
        else if (value == "left")
          post_polarisation = Polariser::circleft_jones ; 
        else if (value == "right")
          post_polarisation = Polariser::circright_jones ; 
        else if (value == "horizontal")
          post_polarisation = Polariser::horiz_jones ; 
        else if (value == "vertical")
          post_polarisation = Polariser::vert_jones ; 
        else 
          printf("Unknown polariser %s\n", value.get<std::string>().c_str()) ; 
      }
      else if (key == "strategy")
      {
        if (value == "LINEAR_NEARESTNEIGHBOUR")
          strategy = Strategy::LINEAR_NEARESTNEIGHBOUR ; 
        else if (value == "LINEAR_TETRAHEDRON_INVERSION")          
          strategy = Strategy::LINEAR_TETRAHEDRON_INVERSION ; 
        else if (value == "LINEAR_TETRAHEDRON_MOLLERTRUMBORE")
          strategy = Strategy::LINEAR_TETRAHEDRON_MOLLERTRUMBORE ; 
        else if (value == "TETRAHEDRON_EXPONENTIAL_INVERSION")
          strategy = Strategy::TETRAHEDRON_EXPONENTIAL_INVERSION ; 
        else if (value == "TETRAHEDRON_EXPONENTIAL_MOLLERTRUMBORE")
          strategy = Strategy::TETRAHEDRON_EXPONENTIAL_MOLLERTRUMBORE ; 
        else 
          printf("Unknown strategy %s.\n", value.get<std::string>().c_str()) ; 
      }      
      else
        printf("Unknown key: %s \n", key.c_str()) ; 
    }
  }
    
  nlohmann::json parameters_from_file(const std::string filename) {
      std::ifstream file(filename); // Open the file
      if (!file.is_open()) { throw std::runtime_error("Unable to open file: " + filename); }
      std::ostringstream ss;
      ss << file.rdbuf(); 
      std::string jsontxt=ss.str();    
    
      // Converting from json5 to json
      jsontxt = std::regex_replace(jsontxt, std::regex(R"(//.*?$)"), "" );
      jsontxt = std::regex_replace(jsontxt, std::regex(R"(/\*[\s\S]*?\*/)"), "");
      jsontxt = std::regex_replace(jsontxt, std::regex(R"('([^'\\]*(?:\\.[^'\\]*)*)')"), R"("$1")");
      jsontxt = std::regex_replace(jsontxt, std::regex(R"(,\s*([\]}]))"), "$1");
      jsontxt = std::regex_replace(jsontxt, std::regex(R"(^\s+|\s+$)"), "");
  
      
      nlohmann::json j ;
      try 
      {
        j = nlohmann::json::parse(jsontxt); 
        parse_json(j) ; 
      }
      catch (...) {
        std::cout << jsontxt ;  std::cout << "ERR: json parsing failed. Here is what we got: " << j ; fflush(stdout) ;}
      
      return j ; 
  }
  
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
