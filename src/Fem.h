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
  void prepare_mesh (std::string filename, std::vector<double> scalings) ; 
  //void get_sigma(std::vector<double> &result, const std::vector<std::array<double,3>> &contact_points, const std::vector<std::array<double,3>> & forces) ; 
  std::vector<int> interpolate(std::vector<vec> & pos, vec graincenter, double radius) ;
  //int interpolate_nearestneighbour(vec & pos) ; 
  //std::vector<int> get_closest_points (const std::vector<std::array<double,3>> & contact_points) ;
  
  struct mesh_wrapper {
    double scaling ; 
    std::shared_ptr<mesh::Mesh<double>> mesh ;
    std::vector<double> com ;
    std::vector<Tetrahedron> tetras ;
    int ncell ;  
    Cells cells ; 
    int interpolate_nearestneighbour_cell (vec & pos) ;
  } ; 
  std::vector<mesh_wrapper> mesh_scaled; 
  int selected_mesh = -1 ; 
  
  int select_mesh(double radius)
  {
    auto it = find_if(mesh_scaled.begin(), mesh_scaled.end(), [=](mesh_wrapper &e){return e.scaling==radius; }) ; 
    if ( it == mesh_scaled.end())
      selected_mesh = -1 ; 
    else
      selected_mesh = it - mesh_scaled.begin() ; 
    return selected_mesh ; 
  }
  
  std::vector<Tetrahedron> & get_tetras(double radius) 
  {
    int s = select_mesh(radius) ; 
    return mesh_scaled[s].tetras ; 
  }
  
  
  template <int N>
  void get_sigma(std::vector<double> &result, const std::vector<std::array<double,3>> &contact_points, const std::vector<std::array<double,3>> & displacements) {}
  
} ; 
template <> void FEsolver::get_sigma<0>([[maybe_unused]] std::vector<double> &result, [[maybe_unused]] const std::vector<std::array<double,3>> &contact_points, [[maybe_unused]] 		const std::vector<std::array<double,3>> & displacements) ; 
template <> void FEsolver::get_sigma<1>(std::vector<double> &result, const std::vector<std::array<double,3>> &contact_points, const std::vector<std::array<double,3>> & displacements) ;
template <> void FEsolver::get_sigma<2>(std::vector<double> &result, const std::vector<std::array<double,3>> &contact_points, const std::vector<std::array<double,3>> & displacements) ;
template <> void FEsolver::get_sigma<3>(std::vector<double> &result, const std::vector<std::array<double,3>> &contact_points, const std::vector<std::array<double,3>> & displacements) ;
template <> void FEsolver::get_sigma<4>(std::vector<double> &result, const std::vector<std::array<double,3>> &contact_points, const std::vector<std::array<double,3>> & displacements) ;
template <> void FEsolver::get_sigma<5>(std::vector<double> &result, const std::vector<std::array<double,3>> &contact_points, const std::vector<std::array<double,3>> & displacements) ; 
template <> void FEsolver::get_sigma<6>(std::vector<double> &result, const std::vector<std::array<double,3>> &contact_points, const std::vector<std::array<double,3>> & displacements) ; 
template <> void FEsolver::get_sigma<7>(std::vector<double> &result, const std::vector<std::array<double,3>> &contact_points, const std::vector<std::array<double,3>> & displacements) ;
#endif



