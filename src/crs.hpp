#ifndef CRS_HPP
#define CRS_HPP

#include "array_ref.hpp"

#include <algorithm>
#include <iostream>
#include <numeric>

////////////////////////////////////////////////////////////////////////////////
/// This type is a container for compressed-storage of sparse data.
///
/// @var offsets The offset at which each index begins and ends in the
///              list of indices.
/// @var indices The indices of the sparse structure of the data resolved
///              by this storage member.
///
/// @ingroup coloring
////////////////////////////////////////////////////////////////////////////////
template< typename T, typename U >
struct crs_t {

  using value_type = T;
  using index_type = U;

  std::vector<index_type> offsets;
  std::vector<value_type> indices;
  
  bool empty() const { return offsets.empty(); }

  index_type size() const {
    if(offsets.empty())
      return 0;
    else
      return offsets.size() - 1;
  } // size

  index_type size(index_type i) const
  {
    assert(i < size() && "index out of range");
    return offsets[i+1] - offsets[i];
  }

  //============================================================================
  /// \brief erase a bunch of ids
  //============================================================================
  template<typename V>
  void erase(const std::vector<V> & ids) {

    if(ids.empty())
      return;

    // assume sorted
    assert(std::is_sorted(ids.begin(), ids.end()) &&
           "entries to delete are not sorted");

    auto num_offsets = size();
    index_type num_erase = 0;
    
    auto delete_it = ids.begin();

    for ( index_type iold=0, inew=0; iold<num_offsets; ++iold ) {
    
      // skip deleted items
      if ( delete_it != ids.end() )
      {
        if ( static_cast<index_type>(*delete_it) == iold )
        {
          delete_it++;
          num_erase++;
          continue;
        }
      }

      // keep otherwise
      auto old_start = offsets[iold];
      auto old_end = offsets[iold+1];
      auto n = old_end - old_start;
      auto new_start = offsets[inew];
      for ( index_type j=0; j<n; ++j )
        indices[new_start + j] = indices[old_start + j];
      offsets[inew+1] = new_start + n;
      inew++;
    }
  
    num_offsets -= num_erase; 
    offsets.resize( num_offsets+1 );
    indices.resize( offsets.back() );
  }
  
  //============================================================================
  /// \brief reorder according to the provided 
  //============================================================================
  template<typename RandomIt>
  void reorder(RandomIt order) {

    if (offsets.empty()) return;

    index_type num_offsets = offsets.size();
    index_type num_indices = indices.size();

    std::vector<index_type> new_offsets(num_offsets);
    std::vector<value_type> new_indices(num_indices);

    new_offsets[0] = 0;
    for (index_type i=0; i<num_offsets-1; ++i) {
      auto iold = order[i];
      auto n = offsets[iold+1] - offsets[iold];
      new_offsets[i+1] = n;
    }

    std::partial_sum(new_offsets.begin(), new_offsets.end(), new_offsets.begin());

    for (index_type i=0; i<num_offsets-1; ++i) {
      auto iold = order[i];
      auto old_start = offsets[iold];
      auto n = offsets[iold+1] - old_start;
      auto new_start = new_offsets[i];
      for (index_type j=0; j<n; ++j)
        new_indices[new_start + j] = indices[old_start + j];
    }

    std::swap(indices, new_indices);
    std::swap(offsets, new_offsets);
  }


  /// \brief clears the current storage
  void clear() {
    offsets.clear();
    indices.clear();
  }

  void reserve(index_type num_offsets, index_type num_indices)
  {
    offsets.reserve(num_offsets+1);
    indices.reserve(num_indices);
  }

  //============================================================================
  /// Non-const iterator
  //============================================================================
  class iterator
  {
    crs_t & data_;
    index_type pos_;

  public:
    iterator(crs_t & data, index_type pos = 0) : data_(data), pos_(pos) {}
    auto operator*() {
      auto i = data_.offsets[pos_];
      auto n = data_.offsets[pos_ + 1] - i;
      return make_array_ref(&data_.indices[i], n);
    }
    auto & operator++() {
      pos_++;
      return *this;
    } // prefix
    auto operator++(int) {
      auto i = *this;
      pos_++;
      return i;
    } // postfix
    bool operator!=(const iterator & it) const {
      return (pos_ != it.pos_ || &data_ != &it.data_);
    }
  };

  //============================================================================
  /// Const iterator
  //============================================================================
  class const_iterator
  {
    const crs_t & data_;
    index_type pos_;

  public:
    const_iterator(const crs_t & data, index_type pos = 0)
      : data_(data), pos_(pos) {}

    auto operator*() {
      auto i = data_.offsets[pos_];
      auto n = data_.offsets[pos_ + 1] - i;
      return make_array_ref(&data_.indices[i], n);
    }
    auto & operator++() {
      pos_++;
      return *this;
    } // prefix
    auto operator++(int) {
      auto i = *this;
      pos_++;
      return i;
    } // postfix
    bool operator!=(const const_iterator & it) const {
      return (pos_ != it.pos_ || &data_ != &it.data_);
    }
  };

  //============================================================================
  // Accessors
  //============================================================================
  auto begin() {
    return iterator(*this);
  }
  auto end() {
    return iterator(*this, size());
  }

  auto begin() const {
    return const_iterator(*this);
  }
  auto end() const {
    return const_iterator(*this, size());
  }

  auto at(index_type i) {
    assert(i < size() && "index out of range");
    return *iterator(*this, i);
  }
  auto at(index_type i) const {
    assert(i < size() && "index out of range");
    return *const_iterator(*this, i);
  }
  auto operator[](index_type i) {
    return *iterator(*this, i);
  }
  auto operator[](index_type i) const {
    return *const_iterator(*this, i);
  }

  //============================================================================
  // append data
  //============================================================================
  template<typename InputIt>
  void push_back(InputIt first, InputIt last) {
    if(first == last)
      return;
    if(offsets.empty())
      offsets.emplace_back(0);
    offsets.emplace_back(offsets.back() + std::distance(first, last));
    indices.insert(indices.end(), first, last);
  }

  template<typename V>
  void push_back(const V & value) {
    auto ptr = &value;
    push_back(ptr, ptr + 1);
  }

  template<typename V>
  void push_back(std::initializer_list<V> init) {
    push_back(init.begin(), init.end());
  }

  template<typename... Args, template<typename...> class Vector>
  void push_back(const Vector<Args...> & init) {
    push_back(init.begin(), init.end());
  }

  template<typename V, typename W>
  void assign(const crs_t<V, W> & other)
  {
    offsets.assign(other.offsets.begin(), other.offsets.end());
    indices.assign(other.indices.begin(), other.indices.end());
  }

}; // struct crs

//============================================================================
/// Helper function to print a crs instance.
//============================================================================

template< typename T, typename U >
inline std::ostream &
operator<<(std::ostream & stream, const crs_t<T, U> & c) {
  stream << "offsets: ";
  for(auto i : c.offsets) {
    stream << i << " ";
  } // for
  stream << std::endl;

  stream << "indices: ";
  for(auto i : c.indices) {
    stream << i << " ";
  } // for

  return stream;
} // operator <<

#endif
