#ifndef GRAINS_H
#define GRAINS_H





class Grains {
public: 
  Grains(vec location, double radius) {center=location ; r = radius ; }
  vec center ; 
  double r ;
  std::vector<std::array<double, 3>>  contactpoints ;
  //std::vector<std::array<double, 3>>  forces ; 
  std::vector<std::array<double, 3>>  displacements ; 
  
  std::vector<double> stress ; 
  
  std::vector<tens> stress_at(std::vector<int> & idx)
  {
    std::vector<tens> res ; 
    res.resize(idx.size()) ; 
    for (size_t i=0 ; i<idx.size() ; i++)
      for (int j=0 ; j<9 ; j++)  
        res[i][j] = stress[idx[i]*9+j]; 
    return res ; 
  }
  
  
  std::vector<tens> stress_at(std::vector<vec> loc)
  {
    // PLACEHOLDER
    std::vector<tens> res ; 
    res.resize(loc.size()) ; 
    for (size_t i=0 ; i<res.size() ; i++)
    {
      res[i]=tens({1,0,0,0,2,0,0,0,3}) ; 
    }
    return res ; 
  }
  
} ;

#endif
