#include "Fem.h"

void FEsolver::prepare_mesh (std::string filename, double scaling) 
{
  io::XDMFFile meshfile(MPI_COMM_WORLD, filename.c_str(), "r");
  auto cel = fem::CoordinateElement<double>(mesh::CellType::tetrahedron, 1, basix::element::lagrange_variant::unset);
  mesh = std::make_shared<mesh::Mesh<double>>(meshfile.read_mesh(cel, mesh::GhostMode::shared_facet, "sphere"));
  
  
  // Calculate cells center of mass
  auto meshpoints=mesh->geometry().x() ;
  if (scaling != 1)
    for (size_t i=0 ; i<meshpoints.size() ; i++)
      meshpoints[i] *= scaling ;
  auto topology=mesh->topology()->connectivity(3,0)->array() ; 
      
  const int CONN = 4 ; // TODO do not assume tetrahedron
  const int DIM = 3 ; 
  ncell=topology.size()/CONN ; 
  com.resize(ncell*DIM,0) ; 
  
  std::vector<double> dst ; dst.resize(6) ; 
  double dstmax=0 ; 
  for (int i=0 ; i<ncell ; i++)
  {
    auto c = [&](int pt, int dim){return meshpoints[topology[i*CONN+pt]*DIM+dim] ; };
    for (int j=0 ; j<CONN ; j++)
      for (int dd=0 ; dd<DIM ; dd++)
        com[i*DIM+dd] += meshpoints[topology[i*CONN+j]*DIM+dd]/4. ; 
    
    dst[0] = (c(0,0)-c(1,0))*(c(0,0)-c(1,0))+(c(0,1)-c(1,1))*(c(0,1)-c(1,1))+(c(0,2)-c(1,2))*(c(0,2)-c(1,2)) ; 
    dst[1] = (c(0,0)-c(2,0))*(c(0,0)-c(2,0))+(c(0,1)-c(2,1))*(c(0,1)-c(2,1))+(c(0,2)-c(2,2))*(c(0,2)-c(2,2)) ; 
    dst[2] = (c(0,0)-c(3,0))*(c(0,0)-c(3,0))+(c(0,1)-c(3,1))*(c(0,1)-c(3,1))+(c(0,2)-c(3,2))*(c(0,2)-c(3,2)) ; 
    dst[3] = (c(1,0)-c(2,0))*(c(1,0)-c(2,0))+(c(1,1)-c(2,1))*(c(1,1)-c(2,1))+(c(1,2)-c(2,2))*(c(1,2)-c(2,2)) ; 
    dst[4] = (c(1,0)-c(3,0))*(c(1,0)-c(3,0))+(c(1,1)-c(3,1))*(c(1,1)-c(3,1))+(c(1,2)-c(3,2))*(c(1,2)-c(3,2)) ; 
    dst[5] = (c(2,0)-c(3,0))*(c(2,0)-c(3,0))+(c(2,1)-c(3,1))*(c(2,1)-c(3,1))+(c(2,2)-c(3,2))*(c(2,2)-c(3,2)) ; 
    if (*std::max_element(dst.begin(), dst.end())>dstmax)
      dstmax = *std::max_element(dst.begin(), dst.end()); 
  }
  std::vector<std::pair<double,double>> boundaries = {{-1,1},{-1,1},{-1,1}} ;
  cells.init_cells(boundaries, 1.1*sqrt(dstmax)) ; 
  cells.allocate_to_cells(com) ; 
  
  // getting connectivity and tetrahedra
  filename = filename+".con" ; 
  FILE *in=fopen(filename.c_str(), "rb") ;
  if (in)
  {
    printf("%s exists, using it for the mesh connectivity.\n", filename.c_str()) ; 
    // Reading the connectivity file
    size_t tmp; size_t res ; 
    res=fread(&tmp, sizeof(size_t), 1, in) ; 
    if (res!=1) printf("Error reading size") ; 
    tetras.resize(tmp) ; 
    for (size_t i=0 ; i<tetras.size() ; i++)
      tetras[i].read(in) ;
    fclose(in); 
  }
  else
  {
    printf("Creating the connectivity file. This may take time, but if you do not change the mesh the .con file created may be reused.\n") ; 
    tetras.resize(ncell) ; 
    auto coord = [&](int cell,int pt,int dim){return meshpoints[(topology[cell*CONN+pt]*DIM+dim)];} ; 
    for (int i=0 ; i<ncell ; i++)
    {
      tetras[i].vertex_id={topology[i*CONN+0], topology[i*CONN+1], topology[i*CONN+2], topology[i*CONN+3]} ; 
      
      tetras[i].triangles[0].M = {coord(i,0,0), coord(i,0,1), coord(i,0,2)} ; 
      tetras[i].triangles[1].M = tetras[i].triangles[0].M ;
      tetras[i].triangles[2].M = tetras[i].triangles[0].M ;
      tetras[i].com = tetras[i].triangles[0].M ; 
      Eigen::Vector3d pt = {coord(i,1,0), coord(i,1,1), coord(i,1,2)} ;
      tetras[i].triangles[0].u = pt - tetras[i].triangles[0].M ; 
      tetras[i].triangles[1].u = pt - tetras[i].triangles[1].M ; 
      tetras[i].triangles[3].M = pt ;
      tetras[i].com += pt ; 
      pt = Eigen::Vector3d(coord(i,2,0), coord(i,2,1), coord(i,2,2)) ;
      tetras[i].triangles[0].v = pt - tetras[i].triangles[0].M ; 
      tetras[i].triangles[2].u = pt - tetras[i].triangles[2].M ; 
      tetras[i].triangles[3].u = pt - tetras[i].triangles[3].M ;
      tetras[i].com += pt ; 
      pt = Eigen::Vector3d(coord(i,3,0), coord(i,3,1), coord(i,3,2)) ;
      tetras[i].triangles[1].v = pt - tetras[i].triangles[1].M ; 
      tetras[i].triangles[2].v = pt - tetras[i].triangles[2].M ; 
      tetras[i].triangles[3].v = pt - tetras[i].triangles[3].M ;
      tetras[i].com += pt ; 
      tetras[i].com /= 4 ;
    }
    
    // Finding connectivity
    std::vector<int> intersection ;
    for (int i=0 ; i<ncell ; i++)
    {
      std::set<int> tetra1{topology[i*CONN+0], topology[i*CONN+1], topology[i*CONN+2], topology[i*CONN+3]} ; 
      for (size_t j=i+1 ; j<topology.size()/CONN ; j++)
      {
        std::set<int> tetra2{topology[j*CONN+0], topology[j*CONN+1], topology[j*CONN+2], topology[j*CONN+3]} ; 
        intersection.clear() ; 
        std::set_intersection(tetra1.begin(), tetra1.end(), tetra2.begin(), tetra2.end(), std::back_inserter(intersection)) ;
        if (intersection.size()==3)
        {
          char tri_i= tetras[i].identify_tri(intersection) ; 
          char tri_j= tetras[j].identify_tri(intersection) ; 
          
          tetras[i].neigh.insert(j) ; 
          tetras[i].surface |= 1<<tri_i ; 
          tetras[j].neigh.insert(i) ; 
          tetras[j].surface |= 1<<tri_j ; 
        }
        if (intersection.size()>=1)
        {
          tetras[i].neigh_extended.insert(j) ; 
          tetras[j].neigh_extended.insert(i) ; 
        }
      }
    }
    for (int i=0 ; i<ncell ; i++)
      tetras[i].surface = (~tetras[i].surface)&(0x0F) ; 
    
    // Writing the connectivity file
    FILE *out=fopen(filename.c_str(), "wb") ;
    size_t tmp = tetras.size() ; 
    fwrite(&tmp, sizeof(size_t), 1, out) ; 
    for (size_t i=0 ; i<tetras.size() ; i++)
      tetras[i].write(out) ;
    fclose(out);      
  }

  printf("Finished preparing the mesh\n") ; fflush(stdout) ; 
}
//------------------------------------------------
std::vector<int> FEsolver::get_closest_points (const std::vector<std::array<double,3>> & contact_points)
{
  std::vector<int> res(contact_points.size(),0) ; 
  std::vector<double> dist(contact_points.size(),0) ; 
  
  auto meshpoints=mesh->geometry().x() ;
  int npt = meshpoints.size()/3 ; 
  
  auto dst = [&](int a,int b){ return ((contact_points[a][0]-meshpoints[b*3+0]) * (contact_points[a][0]-meshpoints[b*3+0]) + 
                               (contact_points[a][1]-meshpoints[b*3+1]) * (contact_points[a][1]-meshpoints[b*3+1]) + 
                               (contact_points[a][2]-meshpoints[b*3+2]) * (contact_points[a][2]-meshpoints[b*3+2])) ; } ; 
                               
  FILE * out ; 
  out = fopen("log_meshgeom", "w") ; 
  for (size_t i=0 ; i<contact_points.size() ; i++) 
  {
    dist[i] = dst(i,0) ;    
    fprintf(out, "%g, %g, %g\n", meshpoints[0], meshpoints[1], meshpoints[2]) ; 
  }
  
  for (int i=1 ; i<npt ; i++)
  {
    fprintf(out, "%g, %g, %g\n", meshpoints[i*3+0], meshpoints[i*3+1], meshpoints[i*3+2]) ; 
    for (size_t j=0 ; j<contact_points.size() ; j++)
      if (dst(j,i) < dist[j])
      {
        dist[j]=dst(j,i) ;
        res[j]=i ; 
      }
  }
  fclose(out) ; 
  return res ;   
}
// ====================================================
// I apologise for this horrible, horrible code ...
//------------------------------------------------
template <>
void FEsolver::get_sigma<1>(std::vector<double> &result, const std::vector<std::array<double,3>> &contact_points, const std::vector<std::array<double,3>> & displacements)
{
#define IDENTIFIER(name,id) name ## id
#include "fem_impl_1.cpp"
double c0=contact_points[0][0], c1=contact_points[0][1] , c2=contact_points[0][2] ;
double d0=displacements[0][0], d1=displacements[0][1] , d2=displacements[0][2] ;
#include "fem_impl_2_0.cpp"
#undef THISBCS
#define THISBCS() {bc0}
#include "fem_impl_3.cpp"  
}
//------------------------------------------------
template <>
void FEsolver::get_sigma<2>(std::vector<double> &result, const std::vector<std::array<double,3>> &contact_points, const std::vector<std::array<double,3>> & displacements)
{
#define IDENTIFIER(name,id) name ## id
#include "fem_impl_1.cpp"
int fid=0 ; 
double c0=contact_points[fid][0], c1=contact_points[fid][1] , c2=contact_points[fid][2] ;
double d0=displacements[fid][0], d1=displacements[fid][1] , d2=displacements[fid][2] ;
#include "fem_impl_2_0.cpp"
fid=1 ;
c0=contact_points[fid][0], c1=contact_points[fid][1] , c2=contact_points[fid][2] ;
d0=displacements[fid][0], d1=displacements[fid][1] , d2=displacements[fid][2] ;
#include "fem_impl_2_1.cpp"
#undef THISBCS
#define THISBCS() {bc0, bc1}
#include "fem_impl_3.cpp"  
}
//------------------------------------------------
template <>
void FEsolver::get_sigma<3>(std::vector<double> &result, const std::vector<std::array<double,3>> &contact_points, const std::vector<std::array<double,3>> & displacements)
{
#define IDENTIFIER(name,id) name ## id
#include "fem_impl_1.cpp"
int fid=0 ; 
double c0=contact_points[fid][0], c1=contact_points[fid][1] , c2=contact_points[fid][2] ;
double d0=displacements[fid][0], d1=displacements[fid][1] , d2=displacements[fid][2] ;
#include "fem_impl_2_0.cpp"
fid=1 ;
c0=contact_points[fid][0], c1=contact_points[fid][1] , c2=contact_points[fid][2] ;
d0=displacements[fid][0], d1=displacements[fid][1] , d2=displacements[fid][2] ;
#include "fem_impl_2_1.cpp"
fid=2 ;
c0=contact_points[fid][0], c1=contact_points[fid][1] , c2=contact_points[fid][2] ;
d0=displacements[fid][0], d1=displacements[fid][1] , d2=displacements[fid][2] ;
#include "fem_impl_2_2.cpp"
#undef THISBCS
#define THISBCS() {bc0, bc1, bc2}
#include "fem_impl_3.cpp"  
}
//------------------------------------------------
template <>
void FEsolver::get_sigma<4>(std::vector<double> &result, const std::vector<std::array<double,3>> &contact_points, const std::vector<std::array<double,3>> & displacements)
{
#define IDENTIFIER(name,id) name ## id
#include "fem_impl_1.cpp"
int fid=0 ; 
double c0=contact_points[fid][0], c1=contact_points[fid][1] , c2=contact_points[fid][2] ;
double d0=displacements[fid][0], d1=displacements[fid][1] , d2=displacements[fid][2] ;
#include "fem_impl_2_0.cpp"
fid=1 ;
c0=contact_points[fid][0], c1=contact_points[fid][1] , c2=contact_points[fid][2] ;
d0=displacements[fid][0], d1=displacements[fid][1] , d2=displacements[fid][2] ;
#include "fem_impl_2_1.cpp"
fid=2 ;
c0=contact_points[fid][0], c1=contact_points[fid][1] , c2=contact_points[fid][2] ;
d0=displacements[fid][0], d1=displacements[fid][1] , d2=displacements[fid][2] ;
#include "fem_impl_2_2.cpp"
fid=3 ;
c0=contact_points[fid][0], c1=contact_points[fid][1] , c2=contact_points[fid][2] ;
d0=displacements[fid][0], d1=displacements[fid][1] , d2=displacements[fid][2] ;
#include "fem_impl_2_3.cpp"
#undef THISBCS
#define THISBCS() {bc0, bc1, bc2, bc3}
#include "fem_impl_3.cpp"  
}
//------------------------------------------------
template <>
void FEsolver::get_sigma<5>(std::vector<double> &result, const std::vector<std::array<double,3>> &contact_points, const std::vector<std::array<double,3>> & displacements)
{
#define IDENTIFIER(name,id) name ## id
#include "fem_impl_1.cpp"
int fid=0 ; 
double c0=contact_points[fid][0], c1=contact_points[fid][1] , c2=contact_points[fid][2] ;
double d0=displacements[fid][0], d1=displacements[fid][1] , d2=displacements[fid][2] ;
#include "fem_impl_2_0.cpp"
fid=1 ;
c0=contact_points[fid][0], c1=contact_points[fid][1] , c2=contact_points[fid][2] ;
d0=displacements[fid][0], d1=displacements[fid][1] , d2=displacements[fid][2] ;
#include "fem_impl_2_1.cpp"
fid=2 ;
c0=contact_points[fid][0], c1=contact_points[fid][1] , c2=contact_points[fid][2] ;
d0=displacements[fid][0], d1=displacements[fid][1] , d2=displacements[fid][2] ;
#include "fem_impl_2_2.cpp"
fid=3 ;
c0=contact_points[fid][0], c1=contact_points[fid][1] , c2=contact_points[fid][2] ;
d0=displacements[fid][0], d1=displacements[fid][1] , d2=displacements[fid][2] ;
#include "fem_impl_2_3.cpp"
fid=4 ;
c0=contact_points[fid][0], c1=contact_points[fid][1] , c2=contact_points[fid][2] ;
d0=displacements[fid][0], d1=displacements[fid][1] , d2=displacements[fid][2] ;
#include "fem_impl_2_4.cpp"
#undef THISBCS
#define THISBCS() {bc0, bc1, bc2, bc3, bc4}
#include "fem_impl_3.cpp"  
}//------------------------------------------------
template <>
void FEsolver::get_sigma<6>(std::vector<double> &result, const std::vector<std::array<double,3>> &contact_points, const std::vector<std::array<double,3>> & displacements)
{
#define IDENTIFIER(name,id) name ## id
#include "fem_impl_1.cpp"
int fid=0 ; 
double c0=contact_points[fid][0], c1=contact_points[fid][1] , c2=contact_points[fid][2] ;
double d0=displacements[fid][0], d1=displacements[fid][1] , d2=displacements[fid][2] ;
#include "fem_impl_2_0.cpp"
fid=1 ;
c0=contact_points[fid][0], c1=contact_points[fid][1] , c2=contact_points[fid][2] ;
d0=displacements[fid][0], d1=displacements[fid][1] , d2=displacements[fid][2] ;
#include "fem_impl_2_1.cpp"
fid=2 ;
c0=contact_points[fid][0], c1=contact_points[fid][1] , c2=contact_points[fid][2] ;
d0=displacements[fid][0], d1=displacements[fid][1] , d2=displacements[fid][2] ;
#include "fem_impl_2_2.cpp"
fid=3 ;
c0=contact_points[fid][0], c1=contact_points[fid][1] , c2=contact_points[fid][2] ;
d0=displacements[fid][0], d1=displacements[fid][1] , d2=displacements[fid][2] ;
#include "fem_impl_2_3.cpp"
fid=4 ;
c0=contact_points[fid][0], c1=contact_points[fid][1] , c2=contact_points[fid][2] ;
d0=displacements[fid][0], d1=displacements[fid][1] , d2=displacements[fid][2] ;
#include "fem_impl_2_4.cpp"
fid=5 ;
c0=contact_points[fid][0], c1=contact_points[fid][1] , c2=contact_points[fid][2] ;
d0=displacements[fid][0], d1=displacements[fid][1] , d2=displacements[fid][2] ;
#include "fem_impl_2_5.cpp"
#undef THISBCS
#define THISBCS() {bc0, bc1, bc2, bc3, bc4, bc5}
#include "fem_impl_3.cpp"  
}
//------------------------------------------------
std::vector<int> FEsolver::interpolate(std::vector<vec> & pos) 
{
  std::vector<int> res ; 
  res.resize(pos.size()) ; 
  for (size_t i=0 ; i<pos.size() ; i++)
    res[i] = interpolate_nearestneighbour_cell(pos[i]) ; 
  return res ; 
}
//------------------------------------------------
int FEsolver::interpolate_nearestneighbour_cell (vec & pos)
{
  return cells.closest(com,pos) ; 
}
//------------------------------------------------
int FEsolver::interpolate_nearestneighbour(vec & pos)
{
  double mindst ; 
  int idmin = 0 ; 
  double dst = 0 ; 
  for (int i=0 ; i<ncell ; i++)
  {
    dst = (pos[0]-com[i*3+0])*(pos[0]-com[i*3+0])+
          (pos[1]-com[i*3+1])*(pos[1]-com[i*3+1])+
          (pos[2]-com[i*3+2])*(pos[2]-com[i*3+2]) ; 
    if (i==0) mindst = dst ; 
    if (dst<mindst)
    {
      mindst=dst ; 
      idmin = i ; 
    }
  }
  return idmin ;     
}
