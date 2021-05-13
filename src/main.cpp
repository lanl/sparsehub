#include "config.hpp"

#include "crs.hpp"

#include <cmath>
#include <cstring>
#include <cctype>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <locale>
#include <map>
#include <set>
#include <sstream>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <vector>

#ifdef ENABLE_FP_EXCEPTIONS
#include <fenv.h>
#endif

#include <exodusII.h>

using ex_index_t = int64_t;

constexpr auto epsilon = std::numeric_limits<double>::epsilon();
constexpr auto digits = std::numeric_limits<double>::digits10;

// linear advection params
constexpr double alpha = 1;
constexpr double coef[] = {0.5, 0.5, 0};
constexpr double mult[] = {1.0, 1.0, 0.0};
constexpr double center[] = {2.0, 2.0, 0.0};

constexpr double final_time = 10.0;
constexpr double initial_deltat = 1.e-6;
constexpr double cfl = 1;

constexpr int max_iters = 100;
constexpr int output_freq = 1;

const std::string prefix = "out";

constexpr int sigfigs = digits;
constexpr int width = digits + 7;

////////////////////////////////////////////////////////////////////////////////
/// Soution struct
////////////////////////////////////////////////////////////////////////////////
struct solution_t {
  double coef[3] = {0, 0, 0};
  double val = 0;
};

////////////////////////////////////////////////////////////////////////////////
/// Print usage
////////////////////////////////////////////////////////////////////////////////
void print_usage() {
  std::cout << std::endl;
  std::cout << "help: sparsehub [-h] [-l LAYOUT] FILENAME" << std::endl;
  std::cout << std::endl;
  std::cout << "Sparsehub solves the advection equation using a sparse";
  std::cout << " matrix representation of the data." << std::endl;
  std::cout << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << " -h \t\t Display this help message." << std::endl; 
  std::cout << " -l LAYOUT \t Specify the sparse layout type "
    << "[0=dense (default), 1=fixed bandwith, 2=compressed row]."
    << std::endl; 
  std::cout << std::endl;
  std::cout << "Arguments:" << std::endl;
  std::cout << " FILENAME \t The mesh file name.";
  std::cout << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Return the wall time counter.
//! \return the wall time
////////////////////////////////////////////////////////////////////////////////
double wall_time(void)
{
  struct timeval time;
  if (gettimeofday(&time,NULL)){
      //  Handle error
      return 0;
  }
  return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if file exists
////////////////////////////////////////////////////////////////////////////////
bool file_exists(const char * filename) {
  std::ifstream file(filename);
  return file.is_open();
}

////////////////////////////////////////////////////////////////////////////////
/// Square of a number
////////////////////////////////////////////////////////////////////////////////
template<typename T>
T sqr(T x) { return x*x; }

///////////////////////////////////////////////////////////////////////////////
// Tack on an iteration number to a string
///////////////////////////////////////////////////////////////////////////////
std::string zero_padded(unsigned int n, unsigned int padding=6)
{
  std::stringstream ss;
  ss << std::setw( padding ) << std::setfill( '0' ) << n;
  return ss.str();
}
  
////////////////////////////////////////////////////////////////////////////////
/// Flip column to row major
////////////////////////////////////////////////////////////////////////////////
template<typename T>
void column2row_major(T * data, size_t nrows, size_t ncols)
{
  std::vector<T> temp(nrows*ncols);
  for (size_t r=0; r<nrows; ++r)
    for (size_t c=0; c<ncols; ++c)
      temp[r*ncols+c] = data[c*nrows+r];

  for (size_t i=0; i<nrows*ncols; ++i)
    data[i] = temp[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the gometry info for a triangular face
////////////////////////////////////////////////////////////////////////////////
void triangle_normal(
    const double * ax,
    const double * bx,
    const double * cx,
    double * nx)
{
  double u[] = { 
    bx[0] - ax[0],
    bx[1] - ax[1],
    bx[2] - ax[2] };
  double v[] = {
    cx[0] - ax[0],
    cx[1] - ax[1],
    cx[2] - ax[2] };

  nx[0] = 0.5 * (u[1]*v[2] - u[2]*v[1]);
  nx[1] = 0.5 * (u[2]*v[0] - u[0]*v[2]);
  nx[2] = 0.5 * (u[0]*v[1] - u[1]*v[0]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the gometry info for a triangular face
////////////////////////////////////////////////////////////////////////////////
void triangle_geometry(
    const double * ax,
    const double * bx,
    const double * cx,
    double * fx,
    double * nx,
    double & a)
{
  constexpr auto third = 1. / 3.;
  fx[0] = third*(ax[0]+bx[0]+cx[0]);
  fx[1] = third*(ax[1]+bx[1]+cx[1]);
  fx[2] = third*(ax[2]+bx[2]+cx[2]);;

  double u[] = { 
    bx[0] - ax[0],
    bx[1] - ax[1],
    bx[2] - ax[2] };
  double v[] = {
    cx[0] - ax[0],
    cx[1] - ax[1],
    cx[2] - ax[2] };

  nx[0] = 0.5 * (u[1]*v[2] - u[2]*v[1]);
  nx[1] = 0.5 * (u[2]*v[0] - u[0]*v[2]);
  nx[2] = 0.5 * (u[0]*v[1] - u[1]*v[0]);

  a = std::sqrt( nx[0]*nx[0] + nx[1]*nx[1] + nx[2]*nx[2] );
}


////////////////////////////////////////////////////////////////////////////////
/// Compute the gometry info for a quad face
////////////////////////////////////////////////////////////////////////////////
void quad_geometry(
    const double * ax,
    const double * bx,
    const double * cx,
    const double * dx,
    double * fx,
    double * nx,
    double & a)
{
  double fm[] = {
    0.25*(ax[0] + bx[0] + cx[0] + dx[0]),
    0.25*(ax[1] + bx[1] + cx[1] + dx[1]),
    0.25*(ax[2] + bx[2] + cx[2] + dx[2])};

  fx[0] = 0;
  fx[1] = 0;
  fx[2] = 0;
  nx[0] = 0;
  nx[1] = 0;
  nx[2] = 0;
  a = 0;

  double n[3], f[3], atmp;

  triangle_geometry(ax, bx, fm, f, n, atmp);
  nx[0] += n[0];
  nx[1] += n[1];
  nx[2] += n[2];
  fx[0] += atmp*f[0];
  fx[1] += atmp*f[1];
  fx[2] += atmp*f[2];
  a += atmp;
  
  triangle_geometry(bx, cx, fm, f, n, atmp);
  nx[0] += n[0];
  nx[1] += n[1];
  nx[2] += n[2];
  fx[0] += atmp*f[0];
  fx[1] += atmp*f[1];
  fx[2] += atmp*f[2];
  a += atmp;
  
  triangle_geometry(cx, dx, fm, f, n, atmp);
  nx[0] += n[0];
  nx[1] += n[1];
  nx[2] += n[2];
  fx[0] += atmp*f[0];
  fx[1] += atmp*f[1];
  fx[2] += atmp*f[2];
  a += atmp;
  
  triangle_geometry(dx, ax, fm, f, n, atmp);
  nx[0] += n[0];
  nx[1] += n[1];
  nx[2] += n[2];
  fx[0] += atmp*f[0];
  fx[1] += atmp*f[1];
  fx[2] += atmp*f[2];
  a += atmp;

  if ( a > epsilon ) {
    double fact = 1/a;
    fx[0] *= fact;
    fx[1] *= fact;
    fx[2] *= fact;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Helper to compute the gometry info for a hex[0] cell
////////////////////////////////////////////////////////////////////////////////
void tet_geometry(
    const double * ax,
    const double * bx,
    const double * cx,
    double * xc,
    double & v)
{
  double n[3];
  triangle_normal(ax, bx, cx, n);
  xc[0] += n[0] * (sqr(ax[0] + bx[0]) + sqr(bx[0] + cx[0]) + sqr(ax[0] + cx[0]));
  xc[1] += n[1] * (sqr(ax[1] + bx[1]) + sqr(bx[1] + cx[1]) + sqr(ax[1] + cx[1]));
  xc[2] += n[2] * (sqr(ax[2] + bx[2]) + sqr(bx[2] + cx[2]) + sqr(ax[2] + cx[2]));
  v += n[0]*cx[0] + n[1]*cx[1] + n[2]*cx[2];
}

////////////////////////////////////////////////////////////////////////////////
/// Helper to compute the gometry info for a hex[0] cell
////////////////////////////////////////////////////////////////////////////////
void hex_geometry(
    const double * ax,
    const double * bx,
    const double * cx,
    const double * dx,
    double * xc,
    double & v)
{
  // face midpoint
  double xm[] = {
    0.25*(ax[0]+bx[0]+cx[0]+dx[0]),
    0.25*(ax[1]+bx[1]+cx[1]+dx[1]),
    0.25*(ax[2]+bx[2]+cx[2]+dx[2])
  };
 
  tet_geometry(ax, bx, xm, xc, v);
  tet_geometry(bx, cx, xm, xc, v);
  tet_geometry(cx, dx, xm, xc, v);
  tet_geometry(dx, ax, xm, xc, v);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the cell geometry
////////////////////////////////////////////////////////////////////////////////
void cell_geometry(
    const ex_index_t * vs,
    const double * vx,
    double * xc,
    double & v)
{
  auto ax = vx + vs[0]*3;
  auto bx = vx + vs[3]*3;
  auto cx = vx + vs[2]*3;
  auto dx = vx + vs[1]*3;
  auto ex = vx + vs[4]*3;
  auto fx = vx + vs[7]*3;
  auto gx = vx + vs[6]*3;
  auto hx = vx + vs[5]*3;

  // x=i
  hex_geometry(
    ax, bx, cx,  dx,
    xc, v);
  // x=i+1
  hex_geometry(
    hx, gx, fx, ex,
    xc, v);
  
  // y=i
  hex_geometry(
    ax, ex, fx, bx,
    xc, v);
  // y=i+1
  hex_geometry(
    dx, cx, gx, hx,
    xc, v);
  
  // z=i
  hex_geometry(
    ax, dx, hx, ex,
    xc, v);
  // z=i+1
  hex_geometry(
    bx, fx, gx, cx,
    xc, v);
  
  auto vabs = std::abs(v);
  
  if ( vabs > epsilon ) {
    auto fact = 1 / (8*v);
    xc[0] *= fact;
    xc[1] *= fact;
    xc[2] *= fact;
  }

  v = vabs / 3;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the face geometry
////////////////////////////////////////////////////////////////////////////////
void face_geometry(
    const ex_index_t * vs,
    const double * vx,
    double * x,
    double * n,
    double & a)
{
  quad_geometry(
      vx + vs[0]*3,
      vx + vs[1]*3,
      vx + vs[2]*3,
      vx + vs[3]*3,
      x, n, a);
  if (a > epsilon) { 
    n[0] /= a;
    n[1] /= a;
    n[2] /= a;
  }
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Transpose a connectivity array.
////////////////////////////////////////////////////////////////////////////////
template<typename T, typename U>
void
transpose(const crs_t<T, U> & in, crs_t<T, U> & out) {
  
  // determine offsets
  size_t num_from = in.size();

  const auto & in_offsets = in.offsets;
  const auto & in_indices = in.indices;
  size_t num_to = *std::max_element(in_indices.begin(), in_indices.end()) + 1;
  
  auto & out_offsets = out.offsets;
  out_offsets.clear();
  out_offsets.resize(num_to+1, 0);
  
  for (size_t from = 0; from < num_from; ++from) {
    for (auto i=in_offsets[from]; i<in_offsets[from+1]; ++i) {
      auto to = in_indices[i];
      out_offsets[to+1]++;
    }
  }

  out_offsets[0] = 0;
  for (size_t to = 0; to < num_to; ++to)
    out_offsets[to+1] += out_offsets[to];


  // determine final connectivity
  auto & out_indices = out.indices;
  out_indices.clear();
  out_indices.resize( out_offsets.back() );

  std::vector<int> pos(num_to, 0);

  for (size_t from = 0; from < num_from; ++from) {
    for (auto i=in_offsets[from]; i<in_offsets[from+1]; ++i) {
      auto to = in_indices[i];
      out_indices[ out_offsets[to] + pos[to] ] = from;
      pos[to]++;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Output vtk solution
////////////////////////////////////////////////////////////////////////////////
void output_vtk(
    const char * prefix,
    const std::vector<double> & node_coords,
    const crs_t<ex_index_t, ex_index_t> & cell2vertices,
    const std::vector<solution_t> & cell_solution)
{
  std::string filename(prefix);
  filename += ".vtk";
  std::ofstream file(filename);
  
  file.setf(std::ios::scientific);
  file.precision(sigfigs);

  auto nverts = node_coords.size() / 3;
  auto ncells = cell_solution.size();

  file << "# vtk DataFile Version 2.0" << std::endl;
  file << "Advection results" << std::endl;
  file << "ASCII" << std::endl;
  file << "DATASET UNSTRUCTURED_GRID" << std::endl;

  file << "POINTS " << nverts << " double" << std::endl;
  for (size_t v=0; v<nverts; ++v) {
    for (size_t d=0; d<3; ++d)
      file << node_coords[v*3 + d] << " ";
    file << std::endl;
  }
  
  file << "CELLS " << ncells << " " << ncells*(1+8) << std::endl;
  for (size_t c=0; c<ncells; ++c) {
    auto vs = cell2vertices.at(c);
    file << vs.size() << " ";
    for (auto v : vs) file << v << " ";
    file << std::endl;
  }
  
  file << "CELL_TYPES " << ncells << std::endl;
  for (size_t c=0; c<ncells; ++c)
    file << 12 << " ";
  file << std::endl;
  
  file << "CELL_DATA " << ncells << std::endl;
  
  file << "SCALARS solution double 1" << std::endl;
  file << "LOOKUP_TABLE default" << std::endl;
  for (size_t c=0; c<ncells; ++c) 
    file << cell_solution[c].val << " ";
  file << std::endl;
  
  file << "VECTORS coefficient double" << std::endl;
  for (size_t c=0; c<ncells; ++c) {
    const auto & coef = cell_solution[c].coef;
    file << coef[0] << " " << coef[1] << " " << coef[2] << std::endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Output csv solution
////////////////////////////////////////////////////////////////////////////////
void output_csv(
    const char * prefix,
    const std::vector<double> & cell_centroids,
    const std::vector<solution_t> & cell_solution)
{
  std::string filename(prefix);
  filename += ".csv";
  std::ofstream file(filename);
  
  file.setf(std::ios::scientific);
  file.precision(sigfigs);

  auto ncells = cell_solution.size();
  auto ndims = cell_centroids.size() / ncells;

  for (size_t d=0; d<ndims; ++d) {
    std::stringstream ss;
    ss << "x" << d;
    file << std::setw(width) << ss.str() << ",";
  }
  for (size_t d=0; d<ndims; ++d) {
    std::stringstream ss;
    ss << "coefficient" << d;
    file << std::setw(width) << ss.str() << ",";
  }
  file << std::setw(width) << "solution" << std::endl; 

  for (size_t c=0; c<ncells; ++c) {
    
    for (size_t d=0; d<ndims; ++d) {
      auto & x = cell_centroids[c*ndims + d];
      file << std::setw(width) << x << ",";
    }

    const auto & u = cell_solution[c];
    for (size_t d=0; d<ndims; ++d)
      file << std::setw(width) << u.coef[d] << ",";
    file << std::setw(width) << u.val << std::endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Compute time step size
////////////////////////////////////////////////////////////////////////////////
double time_step(
    size_t ndims,
    size_t ncells,
    const ex_index_t * c2f_offsets,
    const ex_index_t * c2f_indices,
    const double * face_areas,
    const double * face_normals,
    const double * cell_volumes,
    const solution_t * cell_solution)
{
  double dt = 0;
  
  for (size_t c=0; c<ncells; ++c) {
    auto start = c2f_offsets[c];
    auto end = c2f_offsets[c+1];

    for (auto i=start; i<end; ++i) {
      auto f = c2f_indices[i];

      const auto n = &face_normals[f*ndims];
      const auto a = &cell_solution[c].coef[0];

      double dot = 0;
      for (size_t d=0; d<ndims; ++d)
        dot += n[d] * a[d];

      auto dx = cell_volumes[c] / face_areas[f];
      auto dtc = dot / dx;

      dt = std::max(dt, dtc);
    }
  }

  return 1 / dt;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute flux
////////////////////////////////////////////////////////////////////////////////
double flux(
    const solution_t & ul,
    const solution_t & ur,
    const double * n,
    size_t ndims)
{
  const auto al = &ul.coef[0];
  const auto ar = &ur.coef[0];

  // dot product
  double dot = 0;
  for (size_t d=0; d<ndims; ++d)
    dot += 0.5 * (al[d] + ar[d]) * n[d];

  // upwind
  double u = dot > 0 ? ul.val : ur.val;

  return u*dot;
}

    
////////////////////////////////////////////////////////////////////////////////
/// Compute residual
////////////////////////////////////////////////////////////////////////////////
void residual(
    size_t ndims,
    size_t nfaces,
    size_t nbnd,
    const ex_index_t * f2c_offsets,
    const ex_index_t * f2c_indices,
    const ex_index_t * bnd_faces,
    const double * face_normals,
    const double * face_areas,
    const solution_t * cell_solution,
    const solution_t * bnd_solution,
    double * cell_residual)
{

  for (size_t f=0; f<nfaces; ++f) {
    auto start = f2c_offsets[f];
    auto end = f2c_offsets[f+1];
    auto nc = end - start;

    if (nc == 2) {
      auto ileft = f2c_indices[start];
      auto iright = f2c_indices[start+1];

      const auto & ul = cell_solution[ileft];
      const auto & ur = cell_solution[iright];

      const auto n = &face_normals[f*ndims];
      const auto a = face_areas[f];

      auto flx = flux(ul, ur, n, ndims);
      
      cell_residual[ileft]  -= flx*a;
      cell_residual[iright] += flx*a;
    }

  }
  
  for (size_t bnd=0; bnd<nbnd; ++bnd) {
    auto f = bnd_faces[bnd];
    auto start = f2c_offsets[f];
    auto ileft = f2c_indices[start];

    const auto & ul = cell_solution[ileft];
    const auto & ur = bnd_solution[bnd];

    const auto n = &face_normals[f*ndims];
    const auto a = face_areas[f];

    auto flx = flux(ul, ur, n, ndims);
      
    cell_residual[ileft]  -= flx*a;
  }
}
   

////////////////////////////////////////////////////////////////////////////////
/// Update solution
////////////////////////////////////////////////////////////////////////////////
void update(
    size_t ncells,
    double dt,
    const double * cell_volumes,
    const double * cell_residual,
    solution_t * cell_solution)
{
  for (size_t c=0; c<ncells; ++c) {
    auto fact = dt / cell_volumes[c];
    cell_solution[c].val += fact * cell_residual[c];
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Main
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char ** argv) {

#ifdef ENABLE_FP_EXCEPTIONS
  feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
#endif

  //============================================================================
  // Read command line arguments
  //============================================================================
  
  std::string meshfile;
  int layout;

  opterr = 0;

  std::set<char> options = {'h'};
  std::set<char> arguments = {'l'};

  std::string option_string;
  for (auto opt : options) option_string.push_back(opt);
  for (auto opt : arguments) {
    option_string.push_back(opt);
    option_string.push_back(':');
  }

  int c;
  while ((c = getopt(argc, argv, "hl:")) != -1) {
    switch (c) {
    case 'h':
      print_usage();
      return 0;
    case 'l':
      layout = atoi(optarg);
      break;
    case '?':
      if (arguments.find(optopt) != arguments.end()) {
        std::cout << "Option -" << std::string(1,optopt);
        std::cout << " requires an argument." << std::endl;
      }
      else if (isprint(optopt))
        std::cout << "Unknown option `-" << std::string(1,optopt) << "'." << std::endl;
      else {
        std::cout << "Unknown option character '" << static_cast<int>(optopt);
        std::cout << "'." << std::endl;
      }
      print_usage();
      return 1;
    default:
      std::cout << "Unknown character detected" << std::endl;
      print_usage();
      abort();
    }
  }

  if (optind == argc) { 
    std::cout << "No positional argument specified." << std::endl;
    print_usage();
    return 1;
  }
  else if (argc - optind > 1) { 
    std::cout << "Only one positional argument allowed." << std::endl;
    print_usage();
    return 1;
  }

  meshfile = argv[optind];
  
  //============================================================================
  // Read mesh file
  //============================================================================

  std::cout << "Reading mesh file '" << meshfile << "'" << std::endl;

  if(!file_exists(meshfile.c_str())) {
    std::cout << "File does not exist." << std::endl;
    return 1;
  }

  // size of floating point variables used in app.
  int app_word_size = sizeof(double);

  // size of floating point stored in name.
  int exo_word_size = 0;
  // the version number
  float version;

  // open the file
  auto exoid = ex_open(meshfile.c_str(), EX_READ, &app_word_size, &exo_word_size, &version);
  
  //This sets the file to read IDs as 64 bit.  If the file does not have 
  //64 bit IDs, it should have no effect. 
  ex_set_int64_status(exoid, EX_ALL_INT64_API);

  //----------------------------------------------------------------------------
  // Read vertices
  ex_init_params params;
  ex_get_init_ext(exoid, &params);

  auto num_dims = params.num_dim;
  auto num_nodes = params.num_nodes;
  auto num_blocks = params.num_elem_blk;

  if (num_dims != 3) {
    std::cout << "Only 3d meshes are supported." << std::endl;
    return 1;
  }

  std::vector<double> node_coords(num_nodes*num_dims);
  ex_get_coord(
      exoid, 
      &node_coords[0],
      &node_coords[num_nodes],
      &node_coords[2*num_nodes]);

  // transpose the vertices
  column2row_major(&node_coords[0], num_nodes, num_dims);

  //----------------------------------------------------------------------------
  // Read blocks
      
  std::vector<ex_index_t> block_ids(num_blocks);
  ex_get_ids(exoid, EX_ELEM_BLOCK, block_ids.data());
    
  crs_t<ex_index_t, ex_index_t> cell2verts;
  cell2verts.offsets.emplace_back(0);
  
  crs_t<ex_index_t, ex_index_t> cell2faces;

  std::map<std::vector<ex_index_t>, ex_index_t> sorted_verts2faces;
  crs_t<ex_index_t, ex_index_t> face2verts;
  std::vector<ex_index_t> face_owner;

  auto add_face = [&](const auto & vs, auto cid)
  {
    std::vector<ex_index_t> sorted_vs(vs.begin(), vs.end());
    std::sort(sorted_vs.begin(), sorted_vs.end());
    auto res = sorted_verts2faces.emplace(sorted_vs, sorted_verts2faces.size());
    auto fid = res.first->second;
    if (res.second) {
      face2verts.push_back(vs);
      face_owner.push_back(cid);
    }
    return fid;
  };

  size_t num_cells = 0;
    
  for (int b=0; b<num_blocks; ++b) {
    auto block_id = block_ids[b];
  
    //----------------------------------
    // read block

    // params
    ex_block block_params;
    block_params.id = block_id;
    block_params.type = EX_ELEM_BLOCK;
    ex_get_block_param(exoid, &block_params);

    auto num_nodes_per_elem = block_params.num_nodes_per_entry;
    auto num_elem = block_params.num_entry;
          
    if (!(strncmp("hex", block_params.topology, 3) == 0 ||
          strncmp("HEX", block_params.topology, 3) == 0))
    {
      std::cout << "Only support hexes!." << std::endl;
      return 1;
    }
  
    // get block data
    auto & indices = cell2verts.indices;
    auto & offsets = cell2verts.offsets;

    auto start = indices.size();
    auto end = start + num_nodes_per_elem*num_elem;
    indices.resize(end);
    ex_get_conn(exoid, EX_ELEM_BLOCK, block_id, &indices[start], nullptr, nullptr);
  
    // convert indices to 0-based
    for (auto i=start; i<end; ++i) indices[i]--;

    // compute offsets
    start = offsets.size();
    end = start + num_elem;
    offsets.resize(end);
    for (auto i=start; i<end; ++i)
      offsets[i] = offsets[i-1] + num_nodes_per_elem;

    //----------------------------------
    // Loop over cells

    auto cells_end = num_cells + num_elem;

    for (auto c=num_cells; c<cells_end; ++c)
    {
      auto vs = cell2verts.at(c); 
      ex_index_t fs[6];
      fs[0] = add_face( std::vector<ex_index_t>{vs[0], vs[1], vs[5], vs[4]}, c );
      fs[1] = add_face( std::vector<ex_index_t>{vs[1], vs[2], vs[6], vs[5]}, c );
      fs[2] = add_face( std::vector<ex_index_t>{vs[2], vs[3], vs[7], vs[6]}, c );
      fs[3] = add_face( std::vector<ex_index_t>{vs[3], vs[0], vs[4], vs[7]}, c );
      fs[4] = add_face( std::vector<ex_index_t>{vs[0], vs[3], vs[2], vs[1]}, c );
      fs[5] = add_face( std::vector<ex_index_t>{vs[4], vs[5], vs[6], vs[7]}, c );
      cell2faces.push_back(&fs[0], &fs[6]);
    }

    // increment cell counter
    num_cells += num_elem;

  } // blocks
  
  //============================================================================
  // Compute cell quantities
  //============================================================================
  
  std::cout << "Compute cell geometry." << std::endl;
  
  std::vector<double> cell_centroids(num_cells*num_dims);
  std::vector<double> cell_volumes(num_cells);
  
  for (size_t c=0; c<num_cells; ++c) {
    auto vs = cell2verts.at(c);
    cell_geometry(
        &vs[0],
        &node_coords[0],
        &cell_centroids[c*num_dims],
        cell_volumes[c]);
  }
  
  //============================================================================
  // Compute face quantities
  //============================================================================
  
  std::cout << "Compute face geometry." << std::endl;

  size_t num_faces = face2verts.size();

  crs_t<ex_index_t, ex_index_t> face2cells;
  transpose(cell2faces, face2cells);
  
  std::vector<double> face_centroids(num_faces*num_dims);
  std::vector<double> face_normals(num_faces*num_dims);
  std::vector<double> face_areas(num_faces);

  std::vector<ex_index_t> bnd_faces;

  for (size_t f=0; f<num_faces; ++f) {
    auto vs = face2verts.at(f);
    face_geometry(
        &vs[0],
        &node_coords[0],
        &face_centroids[f*num_dims],
        &face_normals[f*num_dims],
        face_areas[f]);
    auto cs = face2cells.at(f);
    if (cs.size()==1) bnd_faces.emplace_back(f);
  }
  
  auto num_bnd_faces = bnd_faces.size();

  //============================================================================
  // Initialize solution
  //============================================================================
  
  std::cout << "Initialize solution." << std::endl;

  std::vector<solution_t> cell_solution(num_cells);
  std::vector<solution_t> bnd_face_solution(num_bnd_faces);
  std::vector<double> cell_residual(num_cells);

  auto ics = [](const auto x, auto & u) {
    u.coef[0] = coef[0];
    u.coef[1] = coef[1];
    u.coef[2] = coef[2];
    auto fact =
      mult[0]*sqr(x[0]-center[0]) +
      mult[1]*sqr(x[1]-center[1]) +
      mult[2]*sqr(x[2]-center[2]);
    u.val = exp( - fact * alpha );
  };

  for (size_t c=0; c<num_cells; ++c)
    ics(&cell_centroids[c*num_dims], cell_solution[c]);

  for (size_t bnd=0; bnd<num_bnd_faces; ++bnd) {
    auto f = bnd_faces[bnd];
    ics(&face_centroids[f*num_dims], bnd_face_solution[bnd]);
  }
  std::cout << "done" << std::endl;
  
  size_t iter=0;
  double time=0;

  size_t output_counter = 0;
  auto do_output = output_freq > 0;
  auto has_cfl = cfl > 0;
  
  if (do_output) {
    std::stringstream ss;
    ss << prefix << zero_padded(output_counter++);
    std::cout << "Outputing: " << ss.str() << std::endl;
    output_vtk(ss.str().c_str(), node_coords, cell2verts, cell_solution);
    output_csv(ss.str().c_str(), cell_centroids, cell_solution);
  }
  
  //============================================================================
  // Iterate
  //============================================================================
    
  std::cout << std::string(70, '=') << std::endl;
  std::cout << std::setw(8) << "Step"
    << std::setw(17) << "Step Size (s)"
    << std::setw(20) << "Physical Time (s)"
    << std::setw(20) << "Elapsed Time (s)"
    << std::endl;
  std::cout << std::string(70, '=') << std::endl;

  auto clock_start = wall_time();

  double step_size = initial_deltat;
  
  for (; iter < max_iters && time < final_time;) {

    if (has_cfl) {
      step_size = time_step(
          num_dims,
          num_cells,
          cell2faces.offsets.data(),
          cell2faces.indices.data(),
          face_areas.data(),
          face_normals.data(),
          cell_volumes.data(),
          cell_solution.data());
      step_size *= cfl;
    }

    step_size = std::min( step_size, final_time - time);

    std::fill(cell_residual.begin(), cell_residual.end(), 0);

    residual(
        num_dims,
        num_faces,
        num_bnd_faces,
        face2cells.offsets.data(),
        face2cells.indices.data(),
        bnd_faces.data(),
        face_normals.data(),
        face_areas.data(),
        cell_solution.data(),
        bnd_face_solution.data(),
        cell_residual.data());
    
    update(
        num_cells,
        step_size,
        cell_volumes.data(),
        cell_residual.data(),
        cell_solution.data());

      iter++; 
      time += step_size;
      
      auto tdelta = wall_time() - clock_start;
      auto ss = std::cout.precision();
      std::cout.setf( std::ios::scientific );
      std::cout.precision(6);
      std::cout << std::setw(8) << iter
        << std::setw(17) << step_size
        << std::setw(20) << time
        << std::setw(20) << tdelta
        << std::endl;
      std::cout.unsetf( std::ios::scientific );
      std::cout.precision(ss);
  
    if (do_output && iter % output_freq == 0) {
      std::stringstream ss;
      ss << prefix << zero_padded(output_counter++);
      std::cout << "Outputing: " << ss.str() << std::endl;
      output_vtk(ss.str().c_str(), node_coords, cell2verts, cell_solution);
      output_csv(ss.str().c_str(), cell_centroids, cell_solution);
    }

  }

  auto tdelta = wall_time() - clock_start;
  std::cout << "Final solution time is " << std::scientific
    << std::setprecision(2) << time << " after " << iter
    << " steps." << std::endl;

  std::cout << "Elapsed wall time for the last " << iter << " steps is "
    << std::setprecision(4) << std::scientific << tdelta << " s." << std::endl;

  return 0;

}

