





class Grains {
public: 
  Grains(vec location, double radius) {center=location ; r = radius ; }
  vec center ; 
  double r ; 
  
  std::vector<tens> stress_at(std::vector<vec> loc)
  {
    // PLACEHOLDER
    std::vector<tens> res ; 
    res.resize(loc.size()) ; 
    return res ; 
  }
  
} ;
