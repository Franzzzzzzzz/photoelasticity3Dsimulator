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
      if (v[dd] < 0 || v[dd]>=n_cell[dd]) {/*printf("WARN: particle outside the cell volume\n") ;*/ fflush(stdout) ; return -1 ; }
      res += cum_n_cell[dd] * v[dd] ; 
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
  std::vector<int> n_cell , cum_n_cell ; 
  std::vector<double> origin ; 
  std::vector<Cell> cells ;
  double delta[3] ;
  int planesize=0 ; 
} ; 

/** @}*/

/*****************************************************************************************************
 *                                                                                                   *
 *                                                                                                   *
 *                                                                                                   *
 * IMPLEMENTATIONS                                                                                   *
 *                                                                                                   *
 *                                                                                                   *
 *                                                                                                   *
 * ***************************************************************************************************/

//------------------------------------------------------------------------------------
int Cells::init_cells(std::vector<std::pair<double,double>> &bounds, double ds)
{
  printf("CELL SIZE: %g\n", ds) ; 

  n_cell.resize(3,0);
  cum_n_cell.resize(3+1,0) ; 
  origin.resize(3,0) ;
  
  cum_n_cell[0] = 1 ;  
  for (int i=0 ; i<3 ; i++)
  {
    auto [xmin, xmax] = bounds[i] ; 
    n_cell[i] = floor((xmax-xmin)/ds) ; 
    delta[i] = (xmax-xmin)/n_cell[i] ; 
    if (i>0) 
      cum_n_cell[i] = cum_n_cell[i-1]*n_cell[i-1] ; 
    origin[i] = std::get<0>(bounds[i]) ; 
  }
  cum_n_cell[3] = cum_n_cell[3-1]*n_cell[3-1] ; 
  planesize = cum_n_cell[2] ; 
  
  cells.resize(cum_n_cell[3]) ; 
  for (auto &v: n_cell) printf("%d ", v) ;
  printf("\n") ; 
  for (auto &v: cum_n_cell) printf("%d ", v) ; 
  printf("\n") ;   
  for (int i=0 ; i<3 ; i++)
    printf("%g ", delta[i]) ; 
  printf("\n") ; 
  
  for (size_t i=0 ; i<cells.size() ; i++)
  {
    auto x_base = id2x(i) ; 
    auto x = x_base ;
    
    cells[i].neighbours.reserve(27) ; 
    ternary t ;
    t++ ;
    for ( ; t<27 ; t++)
    {
      bool noadd = false ;
      x = x_base ; 
      
      for (int dd=0 ; dd<3 ; dd++)
      {
        if ( t[dd]==2 )
            noadd = true ; 
        else
          x[dd] += t[dd] ; 
        
        if ( x[dd]<0 ) 
        {
          noadd=true ; 
          break ; 
        }
        else if (x[dd]>=n_cell[dd])
        {
          noadd = true ; 
          break ; 
        }
      }
      if (!noadd)
        cells[i].neighbours.push_back(x2id(x)) ;
    }
    
    sort( cells[i].neighbours.begin(), cells[i].neighbours.end() );
    int bef = cells[i].neighbours.size() ; 
    cells[i].neighbours.erase( unique( cells[i].neighbours.begin(), cells[i].neighbours.end() ), cells[i].neighbours.end() );
    int aft = cells[i].neighbours.size() ; 
    if (bef-aft != 0) printf("ERR: no duplication of cells should happen with the algorithm.\n") ;  
    //cells[i].neighbours.erase(std::remove_if(cells[i].neighbours.begin(), cells[i].neighbours.end(), 
    //                                        [=](size_t x) { return x<=i ; }), cells[i].neighbours.end());
  }
  /*for (int i=0 ; i<cells.size() ; i++)
  {
    printf("%d | ", i) ; 
    for (int j=0 ; j<cells[i].neighbours.size() ; j++)
    {
      printf("%d ", cells[i].neighbours[j]) ; 
    }
    printf("\n") ; 
  }*/

  return 0 ; 
}
//=================================================================
int Cells::allocate_to_cells (std::vector<double> & X)
{
  std::vector<int> v (3,0) ; 
  
  for (size_t i=0 ; i<cells.size() ; i++)
    cells[i].incell.clear() ; 
  
  for (size_t i=0 ; i<X.size()/3 ; i++)
  {
    for (int dd=0 ; dd<3 ; dd++)
      v[dd] = floor((X[i*3+dd]-origin[dd])/delta[dd]) ; 
    int id = x2id(v) ; 
    if (id != -1)
      cells[id].incell.push_back(i) ;
    else      
      printf("Cannot find a cell to allocate the tetrahedron to ...") ; 
  }
  return 0 ; 
}
//=================================================================
int Cells::belonging_to(vec pos)
{
  std::vector<int> v ; v.resize(3) ;  
  for (int dd=0 ; dd<3 ; dd++)
      v[dd] = floor((pos[dd]-origin[dd])/delta[dd]) ; 
  return x2id(v) ; 
}
//----------------------------------
int Cells::closest(std::vector<double> & X, vec pos)
{
  int c = belonging_to(pos) ; 
  double mindst = 100000 ; int minid=-1 ; 
  
  // Inner cell contact
  for (size_t ii=0 ; ii<cells[c].incell.size() ; ii++)
  {
      double sum = 0 ; 
      int i = cells[c].incell[ii] ;
      sum = (X[i*3+0]-pos[0])*(X[i*3+0]-pos[0]) +
            (X[i*3+1]-pos[1])*(X[i*3+1]-pos[1]) +
            (X[i*3+2]-pos[2])*(X[i*3+2]-pos[2]) ;             
      if (sum<mindst) 
      {
        mindst = sum ; 
        minid = i ; 
      }
  }
      
  //neighbours contacts
  for (size_t cc=0 ; cc<cells[c].neighbours.size() ; cc++)
  {
    int c2 = cells[c].neighbours[cc] ; 
    for (size_t jj=0 ; jj<cells[c2].incell.size() ; jj++)
    {
      double sum = 0 ; 
      int i = cells[c2].incell[jj] ;
      sum = (X[i*3+0]-pos[0])*(X[i*3+0]-pos[0]) +
            (X[i*3+1]-pos[1])*(X[i*3+1]-pos[1]) +
            (X[i*3+2]-pos[2])*(X[i*3+2]-pos[2]) ;             
      if (sum<mindst) 
      {
        mindst = sum ; 
        minid = i ; 
      }
    }
  }
  return minid ; 
}

#endif
