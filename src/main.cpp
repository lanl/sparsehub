#include "config.hpp"

#include "crs.hpp"
#include "sparse_layout.hpp"

#include <cmath>
#include <cstring>
#include <cctype>
#include <fstream>
#include <functional>
#include <iostream>
#include <iomanip>
#include <iterator>
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

using index_t = uint64_t;
using layout_t = sparse_layout_t<index_t, int>;

constexpr auto epsilon = std::numeric_limits<double>::epsilon();
constexpr auto digits = std::numeric_limits<double>::digits10;

constexpr int sigfigs = digits;
constexpr int width = digits + 7;

////////////////////////////////////////////////////////////////////////////////
/// Boundary struct
////////////////////////////////////////////////////////////////////////////////
struct boundary_t {
  index_t id;
  index_t ghost;
  boundary_t(index_t fid, index_t g) : id(fid), ghost(g) {}
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
    const index_t * vs,
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
    const index_t * vs,
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
    const crs_t<index_t, index_t> & cell2vertices,
    const std::vector<double> & cell_solution)
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
    file << cell_solution[c] << " ";
  file << std::endl;
  
}

////////////////////////////////////////////////////////////////////////////////
/// Output csv solution
////////////////////////////////////////////////////////////////////////////////
void output_csv(
    const char * prefix,
    const std::vector<double> & cell_centroids,
    const std::vector<double> & cell_solution)
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
  file << std::setw(width) << "solution" << std::endl; 

  for (size_t c=0; c<ncells; ++c) {
    
    for (size_t d=0; d<ndims; ++d) {
      auto & x = cell_centroids[c*ndims + d];
      file << std::setw(width) << x << ",";
    }

    const auto & u = cell_solution[c];
    file << std::setw(width) << u << std::endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Compute time step size
////////////////////////////////////////////////////////////////////////////////
double time_step(
    size_t ndims,
    size_t ncells,
    const layout_t & layout,
    const index_t * c2f_offsets,
    const index_t * c2f_indices,
    const double * face_areas,
    const double * face_normals,
    const double * cell_volumes,
    const double * cell_solution,
    const double * coef)
{
  double dt = 0;
  
  for (size_t c=0; c<ncells; ++c) {
    auto start = c2f_offsets[c];
    auto end = c2f_offsets[c+1];

    auto my_mats = layout.row_size(c);

    for (auto i=start; i<end; ++i) {
      auto f = c2f_indices[i];

      const auto n = &face_normals[f*ndims];

      for (int m=0; m<my_mats; ++m) {
        auto mat = layout.column(c, m);
        const auto a = &coef[mat*ndims];

        double dot = 0;
        for (size_t d=0; d<ndims; ++d)
          dot += n[d] * a[d];

        auto dx = cell_volumes[c] / face_areas[f];
        auto dtc = dot / dx;

        dt = std::max(dt, dtc);
      }
    }
  }

  return 1 / dt;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute flux
////////////////////////////////////////////////////////////////////////////////
double flux(
    const double ul,
    const double ur,
    const double * coef,
    const double * n,
    size_t ndims)
{
  // dot product
  double dot = 0;
  for (size_t d=0; d<ndims; ++d)
    dot += coef[d] * n[d];

  // upwind
  double u = dot > 0 ? ul : ur;

  return u*dot;
}

////////////////////////////////////////////////////////////////////////////////
/// Materiall flux
////////////////////////////////////////////////////////////////////////////////
double material_flux(
    int m,
    size_t ndims,
    const layout_t & layout,
    const index_t ileft,
    const index_t iright,
    const double * coef,
    const double * cell_solution,
    const double * n,
    bool & exists)
{
  double ul;
  double ur;

  auto ml = layout.find_column(ileft,  m);
  auto mr = layout.find_column(iright, m);

  bool has_left = ml!=-1;
  bool has_right = mr!=-1;
  exists = (has_left || has_right);
  if (!exists) return 0;

  if (has_right) {
    auto posr = layout(iright, mr);
    ur = cell_solution[posr];
    if (has_left) ul = 0;
  }

  if (has_left) {
    auto posl = layout(ileft, ml);
    ul = cell_solution[posl];
    if (has_right) ur = 0;
  }

  return flux(ul, ur, &coef[m*ndims], n, ndims);
}
    
////////////////////////////////////////////////////////////////////////////////
/// Compute residual
////////////////////////////////////////////////////////////////////////////////
void residual(
    size_t ndims,
    size_t num_mats,
    size_t nfaces,
    size_t nbnd,
    const layout_t & layout,
    const index_t * f2c_offsets,
    const index_t * f2c_indices,
    const boundary_t * bnd_faces,
    const double * face_normals,
    const double * face_areas,
    const double * coef,
    const double * cell_solution,
    double * cell_residual)
{

  bool exists;

  for (size_t f=0; f<nfaces; ++f) {
    auto start = f2c_offsets[f];
    auto end = f2c_offsets[f+1];
    auto nc = end - start;

    if (nc == 2) {
      auto ileft = f2c_indices[start];
      auto iright = f2c_indices[start+1];
        
      const auto n = &face_normals[f*ndims];
      const auto a = face_areas[f];

      for (size_t m=0; m<num_mats; ++m) {

        auto flx = material_flux(
            m, ndims, layout,
            ileft, iright,
            coef, cell_solution, n,
            exists);
        
        if (exists) {
          cell_residual[ileft *num_mats + m] -= flx*a;
          cell_residual[iright*num_mats + m] += flx*a;
        }

      } // mats

    }

  } // faces
  
  for (size_t bnd=0; bnd<nbnd; ++bnd) {
    const auto & f = bnd_faces[bnd];
    auto start = f2c_offsets[f.id];
    auto ileft = f2c_indices[start];

    const auto n = &face_normals[f.id*ndims];
    const auto a = face_areas[f.id];

    for (size_t m=0; m<num_mats; ++m) {
      auto flx = material_flux(
          m, ndims, layout,
          ileft, f.ghost,
          coef, cell_solution, n,
          exists);
      if (exists)
        cell_residual[ileft *num_mats + m] -= flx*a;
    } // mat
      
  } // bnd
}
   

////////////////////////////////////////////////////////////////////////////////
/// Update solution
////////////////////////////////////////////////////////////////////////////////
void update(
    size_t ncells,
    size_t nmats,
    double dt,
    const layout_t & layout,
    const double * cell_volumes,
    const double * cell_solution,
    double * cell_residual)
{
  for (size_t c=0; c<ncells; ++c) {
    auto fact = dt / cell_volumes[c];
    
    for (size_t m=0; m<nmats; ++m)
      cell_residual[c*nmats+m] *= fact;

    auto my_mats = layout.row_size(c);

    for (int m=0; m<my_mats; ++m) {
      auto mat = layout.column(c, m);
      auto matpos = layout(c,m);
      cell_residual[c*nmats+mat] += cell_solution[matpos];
    }

  }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the solution
////////////////////////////////////////////////////////////////////////////////
void to_sparse(
    size_t ncells,
    size_t nmats,
    const layout_t & layout,
    const double * dense,
    double * sparse)
{
  for (size_t c=0; c<ncells; ++c) {

    auto my_mats = layout.row_size(c);
    auto start = layout(c,0);
    for (int m=0; m<my_mats; ++m) {
      auto mat = layout.column(c, m);
      auto matpos = c*nmats + mat;
      sparse[start+m] = dense[matpos];
    }

  } // cells
}

////////////////////////////////////////////////////////////////////////////////
/// Set the solution
////////////////////////////////////////////////////////////////////////////////
void to_dense(
    size_t ncells,
    size_t nmats,
    const layout_t & layout,
    const double * sparse,
    double * dense)
{
  for (size_t c=0; c<ncells; ++c) {

    auto my_mats = layout.row_size(c);
    auto start = layout(c,0);
    for (int m=0; m<my_mats; ++m) {
      auto mat = layout.column(c, m);
      auto matpos = c*nmats + mat;
      dense[matpos] = sparse[start+m];
    }

  } // cells
}

////////////////////////////////////////////////////////////////////////////////
/// Count materials
////////////////////////////////////////////////////////////////////////////////
void count_materials(
    size_t ncells,
    size_t nmats,
    double tol,
    const double * cell_residual,
    char * mat_active,
    int * cell_num_mats,
    int & max_bandwidth)
{
  max_bandwidth = 0;

  for (size_t c=0; c<ncells; ++c) {

    cell_num_mats[c] = 0;
    for (size_t m=0; m<nmats; ++m) {
      auto matpos = c*nmats+m;
      mat_active[matpos] = cell_residual[matpos] > tol;
      cell_num_mats[c] += mat_active[matpos];
    }

    max_bandwidth = std::max( cell_num_mats[c], max_bandwidth );

  }
}


////////////////////////////////////////////////////////////////////////////////
/// Is this a number
////////////////////////////////////////////////////////////////////////////////
bool is_number(const std::string& str)
{
  for (char const &c : str)
    if (std::isdigit(c) == 0) return false;
  return true;
}

void check_number(const std::string & str) {
  if (!is_number(str)) {
    std::cout << "Value '" << str << "' is not representable ";
    std::cout << " as a number." << std::endl;
    abort();
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Parse input file
////////////////////////////////////////////////////////////////////////////////
int parse_input_file(
    const char * filename,
    std::map<std::string, std::string> & input)
{
  std::ifstream file(filename);
 
  while (!file.eof()) {
    std::string key, value, sep;
    file >> key;
    file >> sep;
    file >> value;
    input[key] = value;
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief split a string using a list of delimeters
//! \param [in] str  the input string
//! \param [in] delim  the list of delimeters
//! \return the list of split strings
////////////////////////////////////////////////////////////////////////////////
std::vector<std::string> split(
  const std::string & str, 
  std::vector<char> delim)
{

  if (str.empty()) return {};

  struct tokens_t : std::ctype<char>
  {
    using ctype_base = std::ctype_base;
    using cctype = std::ctype<char>;
    using ccmask = cctype::mask;
    
    tokens_t(const std::vector<char> & delims) 
      : cctype(get_table(delims)) {}

    static ctype_base::mask const * get_table(
      const std::vector<char> & delims
    ) {
      static const ccmask * const_rc = cctype::classic_table();
      static ccmask rc[cctype::table_size];
      std::memcpy(rc, const_rc, cctype::table_size*sizeof(ccmask));
      for (const auto & d : delims) 
        rc[static_cast<int>(d)] = ctype_base::space;
      return &rc[0];
    }
  };

  std::stringstream ss(str);
  ss.imbue(std::locale(std::locale(), new tokens_t(delim)));
  std::istream_iterator<std::string> begin(ss);
  std::istream_iterator<std::string> end;
  std::vector<std::string> vstrings(begin, end);
  return vstrings;
}
  
////////////////////////////////////////////////////////////////////////////////
/// Get input
////////////////////////////////////////////////////////////////////////////////
template<typename T>
T to_val(const std::string & str);

template<>
double to_val<double>(const std::string & str)
{
  check_number(str);
  return atof(str.c_str());
}

template<>
int to_val<int>(const std::string & str)
{
  check_number(str);
  return atoi(str.c_str());
}

template<>
size_t to_val<size_t>(const std::string & str)
{
  check_number(str);
  return atoll(str.c_str());
}

template<>
std::string to_val<std::string>(const std::string & str)
{ return str; }

template<typename T>
T as_scalar(
    std::map<std::string, std::string> & map,
    const std::string & key,
    const T & defval)
{
  bool exists = map.count(key);
  auto res = exists ? to_val<T>(map.at(key)) : defval;
  std::cout << "'" << key << "' = " << res << (exists?"":" (default)") << std::endl;
  return res;
}

template<typename T>
std::vector<T> as_vector(
    std::map<std::string, std::string> & map,
    const std::string & key,
    const std::vector<T> & defval)
{
  std::vector<T> tmp = defval;
  bool exists = map.count(key);

  if (exists) {
    auto toks = split(map.at(key), {','});
    for (size_t t=0; t<std::min(toks.size(), defval.size()); ++t) 
      tmp[t] = to_val<T>(toks[t]);
  }

  std::cout << "'" << key << "' = [ ";
  for (const auto & i : tmp) std::cout << i << " ";
  std::cout << "]" << (exists?"":" (default)") << std::endl;

  return tmp;
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
  
  int layout = 0;

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

  if (argc - optind > 1) { 
    std::cout << "Only one positional argument allowed." << std::endl;
    print_usage();
    return 1;
  }

  std::string inputfile;
  if (argc - optind==1) inputfile = argv[optind];
  
  //============================================================================
  // Read input file
  //============================================================================
  
  std::map<std::string, std::string> input_map;

  if (!inputfile.empty()) {

    std::cout << "Reading input file '" << inputfile << "'" << std::endl;
  
    if(!file_exists(inputfile.c_str())) {
      std::cout << "File '" << inputfile << "'does not exist." << std::endl;
      return 1;
    }
  
    parse_input_file(inputfile.c_str(), input_map);

    std::cout << std::endl;

  }

  // mesh params
  std::cout << "Mesh params:" << std::endl;
  auto dims = as_vector<size_t>(input_map, "dims", {100, 100, 1});
  auto lower = as_vector<double>(input_map, "lower", {0., 0., 0.});
  auto upper = as_vector<double>(input_map, "upper", {10., 10., 1.});
  
  // linear advection params
  std::cout << std::endl << "Problem params:" << std::endl;
  auto alpha = as_scalar<double>(input_map, "alpha", 1.0);
  auto coef = as_vector<double>(input_map, "coef", {0.5, 0.5, 0});
  auto mult = as_vector<double>(input_map, "mult", {1.0, 1.0, 0.0});
  auto center = as_vector<double>(input_map, "center", {2.0, 2.0, 0.0});
  
  std::cout << std::endl << "Solver params:" << std::endl;
  auto final_time = as_scalar<double>(input_map, "time", 10.0);
  auto initial_deltat = as_scalar<double>(input_map, "dt", 1.e-6);
  auto cfl = as_scalar<double>(input_map, "cfl", 0.5);
  
  auto max_iters = as_scalar<size_t>(input_map, "iterations", 100);
  auto output_freq = as_scalar<int>(input_map, "output_freq", 1);

  auto num_stages = as_scalar<int>(input_map, "stages", 2);
  
  auto prefix = as_scalar<std::string>(input_map, "prefix", "out");

  auto num_mats = as_scalar<int>(input_map, "materials", 1);
  auto tol = as_scalar<double>(input_map, "tolerance", 1e-6);

  auto sparse_layout = layout_t(layout);
  

  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // Create mesh
  
  //------------------------------------
  // Mesh stats

  auto num_dims = dims.size();
  if (num_dims != 3) {
    std::cout << "Only 3d meshes are supported." << std::endl;
    return 1;
  }
  
  std::vector<double> delta(num_dims);
  for ( size_t i=0; i<num_dims; ++i ) {
    if ( upper[i] < lower[i] ) {
      std::cout
        << "Bounding box is invalid: for i=" << i << ", " << lower[i]
        << " should be less than " << upper[i]
        << std::endl;
      return 1;
    }
    delta[i] =  (upper[i] - lower[i]) / dims[i];
  }
  
  size_t num_nodes = 1;
  size_t num_cells = 1;
  size_t num_nodes_per_cell = 1;

  for ( size_t i=0; i<num_dims; ++i ) {
    num_nodes *= dims[i]+1;
    num_cells *= dims[i];
    num_nodes_per_cell *= 2;
  }
  
  std::cout << "Creating mesh:" << std::endl;
  std::cout << "  " << num_nodes << " vertices." << std::endl;
  std::cout << "  " << num_cells << " cells." << std::endl;
  
  
  //------------------------------------
  // Build cells
  
  std::vector<double> node_coords(num_nodes*num_dims);
  
  crs_t<index_t, index_t> cell2verts;
  cell2verts.offsets.emplace_back(0);
  
  crs_t<index_t, index_t> cell2faces;

  std::map<std::vector<index_t>, index_t> sorted_verts2faces;
  crs_t<index_t, index_t> face2verts;
  std::vector<index_t> face_owner;

  auto add_face = [&](const auto & vs, auto cid)
  {
    std::vector<index_t> sorted_vs(vs.begin(), vs.end());
    std::sort(sorted_vs.begin(), sorted_vs.end());
    auto res = sorted_verts2faces.emplace(sorted_vs, sorted_verts2faces.size());
    auto fid = res.first->second;
    if (res.second) {
      face2verts.push_back(vs);
      face_owner.push_back(cid);
    }
    return fid;
  };
  
  auto vertex_id = [&](auto i, auto j, auto k) {
    return i + (dims[0]+1)*(j + (dims[1]+1)*k);
  };
  auto cell_id = [&](auto i, auto j, auto k) {
    return i + dims[0]*(j + dims[1]*k);
  };
  
  for ( size_t k=0; k<dims[2]+1; ++k ) {
    for ( size_t j=0; j<dims[1]+1; ++j ) {
      for ( size_t i=0; i<dims[0]+1; ++i ) {
        auto id = vertex_id(i, j, k);
        node_coords[id*num_dims+0] = lower[0] + i*delta[0];
        node_coords[id*num_dims+1] = lower[1] + j*delta[1];
        node_coords[id*num_dims+2] = lower[2] + k*delta[2];
      }
    }
  }

  std::vector<boundary_t> bnd_faces;

  for ( size_t k=0; k<dims[2]; ++k ) {
    for ( size_t j=0; j<dims[1]; ++j ) {
      for ( size_t i=0; i<dims[0]; ++i ) {

        index_t vs[8];
        vs[0] = vertex_id(i  ,   j, k  );
        vs[1] = vertex_id(i+1,   j, k  );
        vs[2] = vertex_id(i+1, j+1, k  );
        vs[3] = vertex_id(i  , j+1, k  );
        vs[4] = vertex_id(i  ,   j, k+1);
        vs[5] = vertex_id(i+1,   j, k+1);
        vs[6] = vertex_id(i+1, j+1, k+1);
        vs[7] = vertex_id(i  , j+1, k+1);
        cell2verts.push_back(&vs[0], &vs[8]);
      
        index_t fs[6];
        fs[0] = add_face( std::vector<index_t>{vs[0], vs[1], vs[5], vs[4]}, c );
        fs[1] = add_face( std::vector<index_t>{vs[1], vs[2], vs[6], vs[5]}, c );
        fs[2] = add_face( std::vector<index_t>{vs[2], vs[3], vs[7], vs[6]}, c );
        fs[3] = add_face( std::vector<index_t>{vs[3], vs[0], vs[4], vs[7]}, c );
        fs[4] = add_face( std::vector<index_t>{vs[0], vs[3], vs[2], vs[1]}, c );
        fs[5] = add_face( std::vector<index_t>{vs[4], vs[5], vs[6], vs[7]}, c );
        cell2faces.push_back(&fs[0], &fs[6]);

        if (i==0)
          bnd_faces.emplace_back( fs[3], cell_id(dims[0]-1,j,k) );
        if (i==dims[0]-1)
          bnd_faces.emplace_back( fs[1], cell_id(0,j,k) );
        if (j==0)
          bnd_faces.emplace_back( fs[0], cell_id(i,dims[1]-1,k) );
        if (j==dims[1]-1)
          bnd_faces.emplace_back( fs[2], cell_id(i,0,k) );
        if (k==0)
          bnd_faces.emplace_back( fs[4], cell_id(i,j,dims[2]-1) );
        if (k==dims[2]-1)
          bnd_faces.emplace_back( fs[5], cell_id(i,j,0) );
      }
    }
  }

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

  crs_t<index_t, index_t> face2cells;
  transpose(cell2faces, face2cells);
  
  std::vector<double> face_centroids(num_faces*num_dims);
  std::vector<double> face_normals(num_faces*num_dims);
  std::vector<double> face_areas(num_faces);

  for (size_t f=0; f<num_faces; ++f) {
    auto vs = face2verts.at(f);
    face_geometry(
        &vs[0],
        &node_coords[0],
        &face_centroids[f*num_dims],
        &face_normals[f*num_dims],
        face_areas[f]);
  }


  auto num_bnd_faces = bnd_faces.size();

  //============================================================================
  // Initialize solution
  //============================================================================
  
  std::cout << "Initialize solution." << std::endl;

  std::vector<std::vector<int>> matids(num_cells);
  std::vector<std::vector<double>> initial_solution(num_cells);

  std::vector<int> cell_num_mats(num_cells, 0);
  std::vector<char> cell_mat_active(num_cells*num_mats, 0);

  auto ics = [&](const auto m, const auto x, auto & u) {
    double fact = 0;
    for (size_t d=0; d<num_dims; ++d)
      fact += mult[d]*sqr(x[d]-center[m*num_dims + d]);
    u = exp( - fact * alpha );
  };

  for (int m=0; m<num_mats; ++m) {

    for (size_t c=0; c<num_cells; ++c) {
      double sol;
      ics(m, &cell_centroids[c*num_dims], sol);
      if (sol > tol) {
        matids[c].emplace_back(m);
        initial_solution[c].emplace_back(sol);
        cell_num_mats[c]++;
        cell_mat_active[c*num_mats+m] = 1;
      }
    }

  } // mats

  std::vector<index_t> mat_offsets;
  std::vector<int> mat_counts;
  std::vector<int> mat_indices;

  sparse_layout.setup(
      num_mats,
      matids,
      mat_offsets,
      mat_counts,
      mat_indices);

  std::vector<double> cell_solution;
  sparse_layout.compress(matids, initial_solution,  cell_solution);
  

  std::vector<double> cell_residual(num_mats*num_cells);

  size_t iter=0;
  double time=0;

  size_t output_counter = 0;
  auto has_cfl = cfl > 0;
  
  std::function<void(void)> do_output;

  if (output_freq > 0) {
    do_output = [&]() {
      std::stringstream ss;
      ss << prefix << zero_padded(output_counter++);
      auto eff = static_cast<double>(cell_solution.size()) / (num_mats*num_cells);
      std::cout << "Outputing: " << ss.str() << ", Efficiency=" << eff*100.;
      std::cout << " %" <<  std::endl;
      std::vector<double> dense_solution(num_mats*num_cells, 0);
      to_dense(
          num_cells,
          num_mats,
          sparse_layout,
          cell_solution.data(),
          dense_solution.data());
      output_vtk(ss.str().c_str(), node_coords, cell2verts, dense_solution);
      output_csv(ss.str().c_str(), cell_centroids, dense_solution);
    };

    do_output();
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
          sparse_layout,
          cell2faces.offsets.data(),
          cell2faces.indices.data(),
          face_areas.data(),
          face_normals.data(),
          cell_volumes.data(),
          cell_solution.data(),
          coef.data());
      step_size *= cfl;
    }

    step_size = std::min( step_size, final_time - time);

    auto stage_fact = static_cast<double>(1) / num_stages;

    for (int s=0; s<num_stages; ++s) {

      std::fill(cell_residual.begin(), cell_residual.end(), 0);

      residual(
          num_dims,
          num_mats,
          num_faces,
          num_bnd_faces,
          sparse_layout,
          face2cells.offsets.data(),
          face2cells.indices.data(),
          bnd_faces.data(),
          face_normals.data(),
          face_areas.data(),
          coef.data(),
          cell_solution.data(),
          cell_residual.data());
      
      update(
          num_cells,
          num_mats,
          stage_fact*step_size,
          sparse_layout,
          cell_volumes.data(),
          cell_solution.data(),
          cell_residual.data());
  
      int max_bandwidth = 0;
      count_materials(
          num_cells,
          num_mats,
          tol,
          cell_residual.data(),
          cell_mat_active.data(),
          cell_num_mats.data(),
          max_bandwidth);

      if (sparse_layout.needs_resize(num_mats, max_bandwidth, cell_num_mats))
        sparse_layout.resize(num_mats, max_bandwidth, cell_num_mats, cell_solution);
  
      sparse_layout.expand(
          num_mats,
          max_bandwidth,
          cell_num_mats,
          mat_offsets,
          mat_counts,
          mat_indices);
  
      sparse_layout.shuffle(cell_mat_active, cell_solution);
      sparse_layout.reconfigure(cell_mat_active);

      to_sparse(
          num_cells,
          num_mats,
          sparse_layout,
          cell_residual.data(),
          cell_solution.data());

    } // stages

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
  
    if (do_output && iter % output_freq == 0) 
      do_output();

  }

  auto tdelta = wall_time() - clock_start;
  std::cout << "Final solution time is " << std::scientific
    << std::setprecision(2) << time << " after " << iter
    << " steps." << std::endl;

  std::cout << "Elapsed wall time for the last " << iter << " steps is "
    << std::setprecision(4) << std::scientific << tdelta << " s." << std::endl;

  return 0;

}

