#ifndef ARRAY_REF_HPP
#define ARRAY_REF_HPP

#include <assert.h>
#include <array>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
/// An \c array_ref<T> represents an immutable array of \c size()
/// elements of type T.  The storage for the array is *not* owned by
/// the \c array_ref object, and clients must arrange for the backing
/// store to remain live while the \c array_ref object is in use.
///
/// Implicit conversion operations are provided from types with
/// contiguous iterators like \c std::vector, \c std::string, \c
/// std::array, and primitive arrays.  \c array_ref objects are
/// invalidated by any operation that invalidates their underlying
/// pointers.
///
/// One common use for \c array_ref is when passing arguments to a
/// routine where you want to be able to accept a variety of array
/// types.  The usual approach here is to have the client explicitly
/// pass in a pointer and a length.
///
/// \todo The existing \c array_ref classes make the view const. It
/// may be useful to extend that to allow modifications of the
/// referenced array elements, and use \c array_ref<const T> for
/// immutable views.
////////////////////////////////////////////////////////////////////////////////
template <typename T>
class array_ref {

  //using DecayedType_ = typename std::remove_const<T>::type;
  //constexpr bool IsSame_ = std::is_same<T, DecayedType_>::value;
  //using Type_ = typename std::conditional< IsSame_, 

public:
  /// \name types
  /// @{
  using element_type = T;
  using value_type = typename std::remove_cv<T>::type;
  /// \todo Should the pointer type be configurable as a template argument?
  using pointer = element_type *;
  using const_pointer = const value_type *;
  using reference = element_type &;
  using const_reference = const element_type &;
  /// random-access, contiguous iterator type
  using iterator = pointer;
  using const_iterator = const_pointer;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  /// @}
        
  /// \name construct/copy
  /// @{
        
  /// \post <code>empty() == true</code>
  constexpr array_ref() : ptr_(nullptr), length_(0) { }
  constexpr array_ref(const array_ref&) = default;
  array_ref& operator=(const array_ref&) = default;
        
  constexpr array_ref(pointer array, size_type length)
    : ptr_(array), length_(length) { }
  
  constexpr array_ref(pointer first, pointer last)
    : ptr_(first), length_(std::distance(first, last)) { }
        
  // Implicit conversion constructors
  template<class C,
    class = std::enable_if_t<
      std::is_convertible<
        std::remove_pointer_t<
          decltype(
              void(std::declval<C &>().size()),
              std::declval<C &>().data()
          )
        > (*)[],
        T (*)[]
      >::value
    >
  >
  //constexpr array_ref(C & c) : array_ref(std::declval<C&>().data(), std::declval<C&>().size()) {}
  constexpr array_ref(C & c) : array_ref(c.data(), c.size()) {}
  
  template<typename U>
  constexpr array_ref(array_ref<U> a) : array_ref(a.ptr_, a.length_) {}

  /// Friend class
  friend class array_ref<value_type>;
  /// @}
        
  /// \name iterators
  /// @{
  constexpr iterator begin() { return ptr_; }
  constexpr iterator end() { return ptr_ + length_; }

  constexpr const_iterator begin() const { return ptr_; }
  constexpr const_iterator end() const { return ptr_ + length_; }

  constexpr const_iterator cbegin() const { return begin(); }
  constexpr const_iterator cend() const { return end(); }

  reverse_iterator rbegin() { return reverse_iterator(end()); }
  reverse_iterator rend() { return reverse_iterator(begin()); }
  
  const_reverse_iterator rbegin() const { return const_reverse_iterator(end()); }
  const_reverse_iterator rend() const { return const_reverse_iterator(begin()); }
  
  const_reverse_iterator crbegin() const { return rbegin(); }
  const_reverse_iterator crend() const { return rend(); }
  /// @}
        
  /// \name capacity
  /// @{
  constexpr size_type size() const { return length_; }
  constexpr size_type max_size() const {
    return std::numeric_limits<size_type>::max() / sizeof(T);
  }
  constexpr bool empty() const { return length_ == 0; }
  /// @}
        
  /// \name element access
  /// @{
  constexpr const T& operator[](size_type i) const { return ptr_[i]; }
  constexpr T& operator[](size_type i) { return ptr_[i]; }
  constexpr const T& at(size_type i) const {
    // This makes at() constexpr as long as the argument is within the
    // bounds of the array_ref.
    return i >= size() ? throw std::out_of_range("at() argument out of range")
      : ptr_[i];
  }
        
  constexpr const T& front() const { return ptr_[0]; }
  constexpr const T& back() const { return ptr_[length_-1]; }
        
  /// \returns A pointer such that [<code>data()</code>,<code>data() +
  /// size()</code>) is a valid range. For a non-empty array_ref,
  /// <code>data() == &front()</code>.
  constexpr const T* data() const { return ptr_; }
  constexpr T* data() { return ptr_; }
  /// @}
        
  /// \name Outgoing conversion operators
  ///
  /// These functions provide explicit conversions to selected other
  /// contiguous sequence types using those types' iterator-range
  /// constructors.  We provide both explicit conversion operators for
  /// use in variable initialization and short member functions for
  /// use in function calls.
  ///
  /// The operators are \c explicit to avoid accidental O(N)
  /// operations on type mismatches.
  ///
  /// @{
        
  /// \todo Arguably, this conversion should be a std::vector
  /// constructor.
  explicit operator std::vector<value_type>() const {
    return std::vector<value_type>(begin(), end());
  }
  std::vector<value_type> vec() const {
    return std::vector<value_type>(*this);
  }
     
  /// \todo Arguably, this conversion should be a std::basic_string
  /// constructor.
  template<typename traits, typename Allocator>
  explicit operator std::basic_string<T, traits, Allocator>() const {
    return std::basic_string<T, traits, Allocator>(begin(), end());
  }
  std::basic_string<T> str() const {
    return std::basic_string<T>(*this);
  }
        
  /// @}
        
  /// \name mutators
  /// @{
        
  /// \par Effects:
  /// Resets *this to its default-constructed state.
  void clear() { *this = array_ref(); }
        
  /// \par Effects:
  /// Advances the start pointer of this array_ref past \p n elements
  /// without moving the end pointer.
  void remove_prefix(size_type n) {
    assert(length_ >= n);
    ptr_ += n;
    length_ -= n;
  }
        
  /// \par Effects:
  /// Moves the end pointer of this array_ref earlier by \p n elements
  /// without moving the start pointer.
  void remove_suffix(size_type n) {
    assert(length_ >= n);
    length_ -= n;
  }
  /// \par Effects:
  /// <code>remove_suffix(1)</code>
  void pop_back() {
    remove_suffix(1);
  }
  /// \par Effects:
  /// <code>remove_prefix(1)</code>
  void pop_front() {
    remove_prefix(1);
  }
  /// @}


private:
  T* ptr_;
  size_type length_;
        
};

////////////////////////////////////////////////////////////////////////////////
/// \name deducing constructor wrappers
/// \relates std::array_ref
/// \xmlonly <nonmember/> \endxmlonly
///
/// These functions do the same thing as the constructor with the same
/// signature. They just allow users to avoid writing the iterator
/// type.
////////////////////////////////////////////////////////////////////////////////
/// @{
    
template<typename T>
constexpr array_ref<const T> make_array_ref(const T* array, std::size_t length) {
  return array_ref<const T>(array, length);
}

template<typename T>
constexpr array_ref<T> make_array_ref(T* array, std::size_t length) {
  return array_ref<T>(array, length);
}
    
    
template<typename T, std::size_t N>
constexpr array_ref<T> make_array_ref(const T (&a)[N]) {
  return array_ref<const T>(a);
}

template<typename T, std::size_t N>
constexpr array_ref<T> make_array_ref(T (&a)[N]) {
  return array_ref<T>(a);
}
    
    
template<typename T>
array_ref<T> make_array_ref(const std::vector<T>& v) {
  return array_ref<const T>(v);
}

template<typename T>
array_ref<T> make_array_ref(std::vector<T>& v) {
  return array_ref<T>(v);
}

template<typename T, std::size_t N>
array_ref<T> make_array_ref(std::array<T,N>& a) {
  return array_ref<T>(a);
}
    
template<typename T, std::size_t N>
array_ref<T> make_array_ref(const std::array<T,N>& a) {
  return array_ref<const T>(a);
}

    
/// @}
    
#endif
