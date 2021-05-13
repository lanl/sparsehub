#------------------------------------------------------------------------------#
# Set the minimum CMake version
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.2)
cmake_policy(SET CMP0074 NEW)
include(ExternalProject)

# a macro for gather options matching a pattern
macro(GatherVariables VARNAME PATTERN)
  set(${VARNAME} "")
  get_cmake_property(_variableNames VARIABLES)
  foreach (_variableName ${_variableNames})
    if(_variableName MATCHES ${PATTERN})
      list(APPEND ${VARNAME} "-D${_variableName}=${${_variableName}} ")
    endif()
  endforeach()
endmacro()

#------------------------------------------------------------------------------#
# Setup the project
#------------------------------------------------------------------------------#

project(FLECSIM-TPL VERSION 0.1.0)

# cmake module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/../cmake")

option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# copy some pertinent cmake flags
set(COMPILE_FLAGS
  -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
  -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE})

# For non-cmake projects, try to figure out how many processors i can use
include(ProcessorCount)
ProcessorCount(PROCESSOR_COUNT)
if(NOT PROCESSOR_COUNT EQUAL 0)
  math(EXPR PROCESSOR_COUNT "${PROCESSOR_COUNT} / 4")
endif()

if (PROCESSOR_COUNT AND ${PROCESSOR_COUNT} GREATER 0 )
  set(${CMAKE_PROJECT_NAME}_MAKE_ARGS -j${PROCESSOR_COUNT})
endif()


#------------------------------------------------------------------------------#
# LUA
#------------------------------------------------------------------------------#
option(BUILD_LUA "Build lua." ON)
if (BUILD_LUA)
  MESSAGE(STATUS "Building Lua.")

  set( LUA_ARGS )
  if ( UNIX AND NOT APPLE )
    list(APPEND LUA_ARGS "linux")
  elseif ( APPLE )
    list(APPEND LUA_ARGS "macosx")
  else()
    MESSAGE(FATAL_ERROR "Unknown lua platform")
  endif()

  ExternalProject_Add(lua
    URL https://www.lua.org/ftp/lua-5.4.2.tar.gz
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND 
      make "${LUA_ARGS}"
    INSTALL_COMMAND
      make install INSTALL_TOP=<INSTALL_DIR>
    BUILD_IN_SOURCE 1
  )
endif()

#------------------------------------------------------------------------------#
# MPICH
#------------------------------------------------------------------------------#
option(BUILD_MPICH "Build MPICH." ON)
if (BUILD_MPICH)
  
  MESSAGE(STATUS "Building MPICH.")

  # mpich configure doesnt like having F90 and F90FLAGS set
  #set(F90_OLD $ENV{F90})
  #set(F90FLAGS_OLD $ENV{F90FLAGS})
  #unset(ENV{F90})
  #unset(ENV{F90FLAGS})

  set( MPICH_ARGS "--enable-fortran=no" )
  #set( MPICH_ARGS )
  if (BUILD_SHARED_LIBS)
    list( APPEND MPICH_ARGS "--enable-shared=yes" )
  else()
    list( APPEND MPICH_ARGS "--enable-shared=no" )
  endif()

  ExternalProject_Add(mpich
    URL http://www.mpich.org/static/downloads/3.3.2/mpich-3.3.2.tar.gz
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND 
      CC=${CMAKE_C_COMPILER}
      CCFLAGS=${CMAKE_C_FLAGS}
      CXX=${CMAKE_CXX_COMPILER}
      CXXFLAGS=${CMAKE_CXX_FLAGS}
      ./configure --prefix=<INSTALL_DIR> ${MPICH_ARGS}
    BUILD_COMMAND make ${${CMAKE_PROJECT_NAME}_MAKE_ARGS}
    INSTALL_COMMAND make install
    BUILD_IN_SOURCE 1
    LOG_BUILD 1
  )

  # Set F90 and F90FLAGS back
  #set(ENV{F90} ${F90_OLD})
  #set(ENV{F90FLAGS} ${F90FLAGS_OLD})
    
endif()


