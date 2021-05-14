#------------------------------------------------------------------------------#
# Set the minimum CMake version
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.2)

#------------------------------------------------------------------------------#
# Setup the project
#------------------------------------------------------------------------------#

project(SPARSEHUB VERSION 0.1.0)

set(default_build_type "Debug")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()


# cmake module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")

# We need C++ 17
add_executable( sparsehub )

include(CTest) # note: this adds a BUILD_TESTING which defaults to ON


#------------------------------------------------------------------------------#
# configure library
#------------------------------------------------------------------------------#

add_subdirectory(src)

