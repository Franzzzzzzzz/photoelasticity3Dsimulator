#ifndef FEM_H
#define FEM_H

#include "sphere.h"
#include <basix/finite-element.h>
#include <cmath>
#include <dolfinx.h>
#include <dolfinx/fem/Constant.h>
#include <dolfinx/fem/petsc.h>
#include <dolfinx/la/petsc.h>
#include <dolfinx/io/XDMFFile.h>
#include <petscmat.h>
#include <petscsys.h>
#include <petscsystypes.h>
#include <utility>
#include <vector>
#include "Cells.h"

using namespace dolfinx;
using T = PetscScalar;
using U = typename dolfinx::scalar_value_type_t<T>;

class FEsolver{
public:
  FEsolver(int argc, char **argv[]) 
  {
    dolfinx::init_logging(argc, *argv);
    PetscInitialize(&argc, argv, nullptr, nullptr);
  }
  //------------------------------------------------  
  void prepare_mesh (std::string filename, double scaling=1) ; 
  void get_sigma(std::vector<double> &result, const std::vector<std::array<double,3>> &contact_points, const std::vector<std::array<double,3>> & forces) ; 
  std::vector<int> interpolate(std::vector<vec> & pos) ;
  int interpolate_nearestneighbour_cell (vec & pos) ;
  int interpolate_nearestneighbour(vec & pos) ; 
  std::vector<int> get_closest_points (const std::vector<std::array<double,3>> & contact_points) ;
  
  std::shared_ptr<mesh::Mesh<double>> mesh ;
  int ncell ; 
  std::vector<double> com ;
  std::vector<Tetrahedron> tetras ; 
  Cells cells ; 
} ; 

#endif



