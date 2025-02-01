#ifndef CELLCONTACTDETECTION
#define CELLCONTACTDETECTION

#include <boost/random.hpp>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <chrono>
#include <algorithm>
#include <numeric>
#include "ternary.h"
#include "Geometry.h"

/** \brief Individual cell for contact detection
 */
class Cell {
public: 
  std::vector<int> incell ; ///< list of particles in the cell
  std::vector<int> neighbours ;  ///< list of cell neighbours to the current cell
} ;

/** \brief All the cells making the space, with related function for creating the cell array, neighbour arrays etc. 
 * Currently non implemented TODO: no pbc handling, no handling of moving boundaries. 
 */
class Cells {
public: 
  int init_cells (std::vector<std::pair<double,double>> &bounds, double ds) ; 
  //int init_cells (std::vector<Boundary<d>> &boundaries, std::vector<double> & r) {double ds = 2.1*(*std::max_element(r.begin(), r.end())) ; return (init_cells(boundaries, ds)) ; }
  int allocate_to_cells (std::vector<double> & X) ; 
  int closest(std::vector<double> & X, vec pos) ; 
  int belonging_to(vec pos) ;
  //-------------------------------------------
  int x2id (const std::vector<int>&v) 
  {
    int res=0 ; 
    for (int dd=0; dd<3 ; dd++)
    {
      res += cum_n_cell[dd] * ( v[dd]<0? 0 : (v[dd]>=n_cell[dd] ? n_cell[dd]-1 : v[dd])) ; 
    }
    return res ; 
  }
  //-------------------------------------------
  std::vector<int> id2x (int id)  
  {
    std::vector<int> ids (3,0) ;
    for (int dd=0 ; dd<3 ; dd++)
    {
      ids[dd] = id % n_cell[dd] ; 
      id /= n_cell[dd] ; 
    }    
    return ids ;
  }
  //-------------------------------------------
  bool test_idempotent ()
  {
    for (int i=0 ; i<cum_n_cell[3] ; i++)
    {
      if (x2id(id2x(i))!=i)
      {
        printf("ERROR %d\n", i) ; 
        return false ; 
      }
    }
    return true ;         
  }
  //---------------------------------------------
  double dsqr (std::vector<double> &X1, std::vector<double> &X2)
  {
    double sum = 0 ; 
    for (int k=0 ; k<3 ; k++)
      sum += (X1[k]-X2[k])*(X1[k]-X2[k]) ; 
    return sum ; 
  }
  //---------------------------------------------
  void rescale(double scale)
  {
    for (auto &v: origin) v *= scale ; 
    for (auto &v: delta) v *= scale ; 
  }
  
  std::vector<int> n_cell , cum_n_cell ; 
  std::vector<double> origin ; 
  std::vector<Cell> cells ;
  double delta[3] ;
  int planesize=0 ; 
} ; 

#endif
