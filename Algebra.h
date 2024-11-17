#ifndef ALGEBRA_H
#define ALGEBRA_H


class vec {
public :
  vec() = default ; 
  vec(double a, double b, double c) {v[0]=a ; v[1]=b; v[2]=c;}
  double & operator[] (int i) { return v[i] ; }
  std::vector<double> v{0,0,0} ; 
  void normalise() {double l=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]) ; for (int i=0 ; i<3 ; i++) v[i]/=l ; }
} ; 

class vec_stokes {
public :
  vec_stokes() = default ; 
  vec_stokes(double a, double b, double c, double d) {v[0]=a ; v[1]=b; v[2]=c; v[3]=d ; }
  double & operator[] (int i) { return v[i] ; }
  std::vector<double> v{1,0,0,0} ;
  void disp() {for (auto &vv:v) printf("%g ", vv) ;  printf("\n") ; }
  void normalise_coherent()
  {
    double norm = sqrt(v[1]*v[1]+v[2]*v[2]+v[3]*v[3])/v[0] ; 
    v[1]/=norm ; 
    v[2]/=norm ; 
    v[3]/=norm ;     
  }
} ; 

class vec_jones {
public:
  vec_jones() = default ; 
  vec_jones(std::complex<double> a, std::complex<double> b) {v[0]=a ; v[1]=b;}
  std::complex<double> & operator[] (int i) { return v[i] ; }
  void disp() {for (auto &vv:v) std::cout << vv << " " ;  printf("\n") ; }
  std::complex<double> v[2] ;   
  void normalise_coherent()
  {
    double n = sqrt(norm(v[0])+norm(v[1])) ;
    v[0] /= n ; 
    v[1] /= n ; 
  }
} ; 

class tens {
public:
  tens() = default ; 
  tens(double a, double b, double c, double d, double e, double f, double g, double h, double i)
  {
    v[0] = a ; v[1] = b ; v[2] = c;
    v[3] = d ; v[4] = e ; v[5] = f;
    v[6] = g ; v[7] = h ; v[8] = i;
  }
  double & operator[] (int i) { return v[i] ; }
  std::vector<double> v{0,0,0,0,0,0,0,0,0} ; 
  double norm() {return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+
                                 v[3]*v[3]+v[4]*v[4]+v[5]*v[5]+
                                 v[6]*v[6]+v[7]*v[7]+v[8]*v[8]) ; }
  void disp() {for (auto &vv:v) printf("%g ", vv) ;  printf("\n") ; }
  void swap_row (int i, int j) 
  {
    double tmp ; 
    tmp = v[i*3+0] ; v[i*3+0]=v[j*3+0] ; v[j*3+0]=tmp ; 
    tmp = v[i*3+1] ; v[i*3+1]=v[j*3+1] ; v[j*3+1]=tmp ; 
    tmp = v[i*3+2] ; v[i*3+2]=v[j*3+2] ; v[j*3+2]=tmp ; 
  }
} ; 

class mat_stokes{
public:
  mat_stokes() = default ; 
  mat_stokes(std::array<double,16> val)
  {
    for (int i=0 ; i<16 ; i++)
      v[i]=val[i] ;     
  }
  double & operator[] (int i) { return v[i] ; }
  std::vector<double> v{1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1} ; 
  void disp() {for(int i=0 ; i<4 ; i++) { for (int j=0 ; j<4 ; j++) printf("%g ", v[i*4+j]) ; printf("\n") ; }}
} ;

class mat_jones{
public:
  mat_jones() = default ; 
  mat_jones(std::array<std::complex<double>,4> val)
  {
    for (int i=0 ; i<4 ; i++)
      v[i]=val[i] ;     
  }
  std::complex<double> & operator[] (int i) { return v[i] ; }
  std::complex<double> v[4]={1,0, 0,1} ; 
  void disp() {for(int i=0 ; i<4 ; i++) { for (int j=0 ; j<4 ; j++) std::cout << v[i*4+j] << " " ; printf("\n") ; }}
} ;

namespace Polariser {
  static mat_stokes horizontal_lin({0.5, 0.5,0,0,  0.5,0.5,0,0, 0,0,0,0, 0,0,0,0}) ; 
  static mat_stokes vertical_lin({0.5,-0.5,0,0, -0.5,0.5,0,0, 0,0,0,0, 0,0,0,0}) ; 
  static mat_stokes deg45_lin({0.5,0,0.5,0, 0,0,0,0, 0.5,0,0.5,0, 0,0,0,0}) ; 
  
  /*static vec_jones horiz_jones({1,0}) ; 
  static vec_jones vert_jones({0,1}) ; 
  static vec_jones circright_jones({1/sqrt(2), -1i / sqrt(2)}) ; 
  static vec_jones circleft_jones({1/sqrt(2), -1i / sqrt(2)}) ; */
  
  static mat_jones horiz_jones({1,0,0,0}) ; 
  static mat_jones vert_jones({0,0,0,1}) ; 
  //static mat_jones circright_jones({1/sqrt(2), -1i / sqrt(2)}) ; 
  //static mat_jones circleft_jones({1/sqrt(2), -1i / sqrt(2)}) ; 
  
} ;

#endif
