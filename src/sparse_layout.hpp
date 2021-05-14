#ifndef SPARSE_LAYOUT_HPP
#define SPARSE_LAYOUT_HPP

#include <cassert>
#include <iostream>
#include <numeric>
#include <set>
#include <vector>

enum class sparse_layout_type_t {
  dense,
  itpack,
  crs
};

////////////////////////////////////////////////////////////////////////////////
/// Spars layout data structure
////////////////////////////////////////////////////////////////////////////////
template<typename ROW_INDEX, typename COL_INDEX>
struct sparse_layout_t {

  using row_t = ROW_INDEX;
  using col_t = COL_INDEX;

  size_t num_rows = 0;
  size_t num_columns = 0;
  size_t bandwidth = 0;
  sparse_layout_type_t type = sparse_layout_type_t::dense;
  row_t * row_offsets = nullptr;
  col_t * row_counts = nullptr;
  col_t * col_ids = nullptr;

  sparse_layout_t() = default;
  sparse_layout_t(int ty) : type(static_cast<sparse_layout_type_t>(ty)) {}
  sparse_layout_t(sparse_layout_type_t ty) : type(ty) {}
  
  //============================================================================
  /// Setup a sparse layout
  //============================================================================
  void setup(
      size_t ncols,
      const std::vector<std::vector<col_t>> & ids,
      std::vector<row_t> & off,
      std::vector<col_t> & cnts,
      std::vector<col_t> & ind,
      col_t padding = 0)
  {
    num_columns = ncols;
    bandwidth = 0;
    row_offsets = nullptr;
    row_counts = nullptr;
    col_ids = nullptr;
      
    num_rows = ids.size();
  
    off.clear();
    cnts.clear();
    ind.clear();
      
    num_columns += padding;
      
    // max bandwidth
    for (const auto & row : ids)
      bandwidth = std::max<col_t>(row.size(), bandwidth);
  
    bandwidth += padding;
  
  
    //----------------------------------
    // Itpack
    if (type == sparse_layout_type_t::itpack) {
  
      cnts.reserve(num_rows);
      ind.reserve(bandwidth*num_rows);
  
      for (const auto & row : ids) {
        auto n = row.size();
        for (size_t j=0; j<n; ++j)
          ind.emplace_back(row[j]);
        for (size_t j=n; j<bandwidth; ++j)
          ind.emplace_back(-1);
        cnts.emplace_back(n);
      }
  
      col_ids = ind.data();
      row_counts = cnts.data();
  
    }
    //----------------------------------
    // crs
    else if (type == sparse_layout_type_t::crs) {
  
      // count entries
      size_t len = num_rows * padding;
      for (const auto & row : ids)
        len += row.size();
  
      // resize
      ind.reserve(len);
      off.reserve(num_rows + 1);
      cnts.reserve(num_rows);
  
      // fill
      off.emplace_back(0);
      for (const auto & row : ids) {
        auto n = row.size();
        for (size_t j=0; j<n; ++j)
          ind.emplace_back(row[j]);
        for (col_t j=0; j<padding; ++j)
          ind.emplace_back(-1);
        cnts.emplace_back(n);
        off.emplace_back(ind.size());
      }
  
      col_ids = ind.data();
      row_offsets = off.data();
      row_counts = cnts.data();
  
    }
  
  } // make
  
  
  //============================================================================
  /// Re-Setup sparse representation
  //============================================================================
  void expand(
      size_t tot_cols,
      size_t max_bandwidth,
      const std::vector<col_t> & ncols,
      std::vector<row_t> & off,
      std::vector<col_t> & cnts,
      std::vector<col_t> & ind)
  {
    //----------------------------------
    // Itpack
    if (type == sparse_layout_type_t::itpack) {
      resize(tot_cols, max_bandwidth, ncols, ind, -1);
  
      for (size_t i=0; i<num_rows; ++i)
        cnts[i] = ncols[i];
  
      row_counts = cnts.data();
      col_ids = ind.data();
    }
    //----------------------------------
    // crs
    else if (type == sparse_layout_type_t::crs) {
  
      resize(tot_cols, max_bandwidth, ncols, ind, -1);
      
      for (size_t i=0; i<num_rows; ++i) {
        auto n = ncols[i];
        cnts[i] = n;
        off[i+1] = off[i] + ncols[i];
      }
      
      row_offsets = off.data();
      row_counts = cnts.data();
      col_ids = ind.data();
    }
    
    // set some general stuff
    num_columns = tot_cols;
    bandwidth = max_bandwidth;
  
  
  }

  //============================================================================
  /// Re-Setup sparse representation
  //============================================================================
  void reconfigure(std::vector<char> & active)
  {
  
    auto tot_cols = active.size() / num_rows;
    col_t num_cols = num_columns;
  
    //----------------------------------
    // Itpack
    if (type == sparse_layout_type_t::itpack) {
  
      for (size_t i=0; i<num_rows; ++i) {
        auto start = bandwidth*i;
        size_t cnt = 0;
        for (col_t j=0; j<num_cols; ++j) {
          if (active[i*tot_cols + j]) {
            assert(cnt<bandwidth && "itpack out of bounds");
            col_ids[start + cnt] = j;
            cnt++;
          }
        }
      }
  
    }
    //----------------------------------
    // crs
    else if (type == sparse_layout_type_t::crs) {
      
      for (size_t i=0; i<num_rows; ++i) {
        auto start = row_offsets[i];
        row_t cnt = 0;
        for (col_t j=0; j<num_cols; ++j) {
          if (active[i*tot_cols + j]) {
            assert(cnt<row_offsets[i+1]-row_offsets[i] && "crs out of bounds");
            col_ids[start + cnt] = j;
            cnt++;
          }
        }
      }
  
    }
  
  }


  //============================================================================
  /// Compress a vector of vectors
  //============================================================================
  template<typename T>
  void compress(
      std::vector<std::vector<col_t>> & ids,
      std::vector<std::vector<T>> & vals,
      std::vector<T> & nonzeros,
      T zero_val = T())
  {
    nonzeros.clear();

    size_t nvals = 1;
    for (size_t i=0; i<num_rows; ++i) {
      auto ncols = ids[i].size();
      auto len = vals[i].size();
      if (ncols > 0) {
        nvals = len / ncols;
        break;
      }
    }

    //----------------------------------
    // Dense
    if (type == sparse_layout_type_t::dense) {

      nonzeros.resize(num_columns*num_rows*nvals, zero_val);  

      for (size_t i=0; i<num_rows; ++i) {
        auto & row_ids = ids[i];
        auto ncols = row_ids.size();
        auto & row_vals = vals[i];
        for (size_t j=0; j<ncols; ++j) {
          auto start = (i*num_columns + row_ids[j]) * nvals;
          for (size_t k=0; k<nvals; ++k)
            nonzeros[start+k] = row_vals[j*nvals + k];
        }
      }

    }
    //----------------------------------
    // Itpack
    else if (type == sparse_layout_type_t::itpack) {

      nonzeros.resize(bandwidth*num_rows*nvals, zero_val);  

      for (size_t i=0; i<num_rows; ++i) {
        auto & row_ids = ids[i];
        auto ncols = row_ids.size();
        auto & row_vals = vals[i];
        for (size_t j=0; j<ncols; ++j) {
          auto start = (i*bandwidth + j) * nvals;
          for (size_t k=0; k<nvals; ++k)
            nonzeros[start + k] = row_vals[j*nvals + k];
        }
      }

    }
    //----------------------------------
    // crs
    else {
      
      nonzeros.resize(row_offsets[num_rows]*nvals); 

      for (size_t i=0; i<num_rows; ++i) {
        auto row_start = row_offsets[i] * nvals;
        auto & row_vals = vals[i];
        auto ncols = row_vals.size();
        for (size_t j=0; j<ncols; ++j)
          nonzeros[row_start + j] = row_vals[j];
      }
    }

  } // make

  //============================================================================
  /// Accessors
  //============================================================================
  col_t row_size(row_t i) const {
    switch (type) {
      case (sparse_layout_type_t::crs):
        return row_counts[i];
      case (sparse_layout_type_t::itpack):
        return row_counts[i];
      default:
        return num_columns;
    };
  }
  
  size_t operator()(row_t i, col_t j) const {
    switch (type) {
      case (sparse_layout_type_t::crs):
        return row_offsets[i] + j;
      case (sparse_layout_type_t::itpack):
        return bandwidth*i + j;
      default:
        return num_columns*i + j;
    }
  }

  col_t column(row_t i, col_t j) const {
    switch (type) {
      case (sparse_layout_type_t::crs):
        return col_ids[row_offsets[i] + j];
      case (sparse_layout_type_t::itpack):
        return col_ids[bandwidth*i + j];
      default:
        return j;
    }
  }

  col_t find_column(row_t i, col_t ind) const
  {
    switch (type) {
      case (sparse_layout_type_t::crs): {
        auto start = row_offsets[i];
        for (auto j=0; j<row_counts[i]; ++j)
          if (ind == col_ids[start+j])
            return j; 
        return -1;
      }
      case (sparse_layout_type_t::itpack): {
        auto start = bandwidth*i;
        auto end = start + row_counts[i];
        for (auto j=start; j<end; ++j)
          if (ind == col_ids[j])
            return j; 
        return -1;
      }
      default:
        return ind<static_cast<col_t>(num_columns) ? ind : -1;
    }
  }

  template<typename T>
  std::vector<T> columns() {
    std::set<T> cols;
    switch (type) {
      case (sparse_layout_type_t::crs):
        for (size_t i=0; i<num_rows; ++i) {
          auto start = row_offsets[i];
          auto end = start + row_counts[i];
          for (auto j=start; j<end; ++j)
            cols.emplace(col_ids[j]);
        }
        break;
      case (sparse_layout_type_t::itpack):
        for (size_t i=0; i<num_rows; ++i) {
          auto start = bandwidth*i;
          auto end = start + row_counts[i];
          for (auto j=start; j<end; ++j)
            cols.emplace(col_ids[j]);
        }
        break;
      default:
        for (size_t j=0; j<num_columns; ++j)
          cols.emplace(j);
    }
    return std::vector<T>(cols.begin(), cols.end());
  }

  bool is_resizeable() const { return type != sparse_layout_type_t::dense; }

  //============================================================================
  /// Check if a resize is needed
  //============================================================================
  bool needs_resize(
      size_t tot_cols,
      size_t max_bandwidth,
      const std::vector<col_t> & ncols,
      bool shrink = false) const
  {
    //--- dense storage
    if (type == sparse_layout_type_t::dense) {
      return shrink ? tot_cols!=num_columns : tot_cols>num_columns;
    }
    //--- max bandwidth
    else if (type == sparse_layout_type_t::itpack) {
        return shrink ? max_bandwidth!=bandwidth : max_bandwidth>bandwidth;
    }
    //--- compressed row
    else if (type == sparse_layout_type_t::crs) {
      if (shrink) {
        for (size_t i=0; i<num_rows; ++i)
          if (ncols[i]!=row_counts[i])
            return true;
      }
      else {
        for (size_t i=0; i<num_rows; ++i) {
          col_t sz = row_offsets[i+1] - row_offsets[i];
          if (ncols[i] > sz)
            return true;
        }
      }
    }
    //--- no resize needed
    return false;
  }


  //============================================================================
  /// Blindly resize an array without checking if needed
  //============================================================================
  template<typename T>
  void resize(
      size_t num_mats,
      size_t max_mats,
      const std::vector<col_t> & ncols,
      std::vector<T> & nonzeros,
      T zero_val = T()) const
  {

    auto nnz = nonzeros.size();

    //----------------------------------
    // Dense or itpack
    if (type == sparse_layout_type_t::dense || type == sparse_layout_type_t::itpack) {

      size_t ncol = (type == sparse_layout_type_t::dense ? num_columns : bandwidth);
      auto nvals = nnz / (num_rows*ncol);
      auto ncol_new = nvals * (type == sparse_layout_type_t::dense ? num_mats : max_mats);
      auto ncol_old = nvals*ncol;
      auto ncol_mid = std::min(ncol_new, ncol_old);
          
      std::vector<T> new_nonzeros(num_rows * ncol_new);

      for (size_t i=0; i<num_rows; ++i) {
        for (size_t j=0; j<ncol_mid; ++j)
          new_nonzeros[i*ncol_new + j] = nonzeros[i*ncol_old + j];
        for (size_t j=ncol_mid; j<ncol_new; ++j)
          new_nonzeros[i*ncol_new + j] = zero_val;
      }

      std::swap(nonzeros, new_nonzeros);

    }
    
    //----------------------------------
    // crs
    else if (type == sparse_layout_type_t::crs) {
      
      auto nvals = nnz / row_offsets[num_rows];
      size_t new_nnz = std::accumulate(ncols.begin(), ncols.end(), 0);

      std::vector<T> new_nonzeros(new_nnz*nvals);
      
      // shrink
      for (size_t i=0, new_start=0; i<num_rows; ++i) {
        auto old_start = row_offsets[i];
        auto nold = row_counts[i];
        auto nnew = ncols[i];
        auto nmid = std::min<col_t>(nold, nnew);
        for (col_t j=0; j<nmid; ++j) {
          auto new_off = (new_start + j) * nvals;
          auto old_off = (old_start + j) * nvals;
          for (size_t k=0; k<nvals; ++k)
            new_nonzeros[new_off + k] = nonzeros[old_off + k];
        }
        for (auto j=nmid; j<nnew; ++j) {
          auto new_off = (new_start + j) * nvals;
          for (size_t k=0; k<nvals; ++k)
            new_nonzeros[new_off + k] = zero_val;
        }
        new_start += nnew;
      }
      
      std::swap(nonzeros, new_nonzeros);
    } // type

  }
  
  //============================================================================
  /// Shuffle 
  //============================================================================
  template<typename T>
  void shuffle(
      const std::vector<char> & active,
      std::vector<T> & nonzeros,
      T zero_val = T()) const
  {

    auto nnz = nonzeros.size();
    size_t tot_cols = active.size() / num_rows;

    //----------------------------------
    // itpack
    if (type == sparse_layout_type_t::itpack) {

      size_t nvals = nnz / (num_rows*bandwidth);

      std::vector<T> row(tot_cols * nvals);
          
      for (size_t i=0; i<num_rows; ++i) {
        auto start = bandwidth*i;
        auto end = start + row_counts[i];

        // expand this row and copy non-zeros
        std::fill(row.begin(), row.end(), zero_val);

        for (auto j=start; j<end && col_ids[j]!=-1; ++j) {
          auto colpos = col_ids[j];
          for (size_t k=0; k<nvals; ++k)
            row[colpos*nvals + k] = nonzeros[j*nvals + k];
        }

        // copy into their new location
        for (size_t j=0, cnt=0; j<tot_cols; ++j) {
          if (active[i*tot_cols + j]) {
            auto pos = start + cnt;
            for (size_t k=0; k<nvals; ++k)
              nonzeros[pos*nvals + k] = row[j*nvals + k];
            cnt++;
          }
        }
      }

    }
    
    //----------------------------------
    // crs
    else if (type == sparse_layout_type_t::crs) {
      
      size_t nvals = nnz / row_offsets[num_rows];
      
      std::vector<T> row(tot_cols * nvals);
      
      // shrink
      for (size_t i=0; i<num_rows; ++i) {
        auto start = row_offsets[i];
        auto end = start + row_counts[i];
        
        // expand this row and copy non-zeros
        std::fill(row.begin(), row.end(), zero_val);

        for (auto j=start; j<end && col_ids[j]!=-1; ++j) {
          auto colpos = col_ids[j];
          for (size_t k=0; k<nvals; ++k)
            row[colpos*nvals + k] = nonzeros[j*nvals + k];
        }

        // copy into their new location
        for (size_t j=0, cnt=0; j<tot_cols; ++j) {
          if (active[i*tot_cols + j]) {
            auto pos = start + cnt;
            for (size_t k=0; k<nvals; ++k)
              nonzeros[pos*nvals + k] = row[j*nvals + k];
            cnt++;
          }
        }

      }

    } // type

  }
  

};

#endif
