#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

# required includes
include(ProcessorCount)

# get the scripts directory
set(REGRESSION_CMAKE_DIR ${CMAKE_CURRENT_LIST_DIR} )

#-------------------------------------------------------------------------------
# This macro creates a name
#-------------------------------------------------------------------------------
function(name_from_path path out)
  file(RELATIVE_PATH ${out} ${PROJECT_SOURCE_DIR} ${path})
  file(TO_NATIVE_PATH ${${out}} ${out})
  string(REPLACE "/" "_" ${out} ${${out}})
  set(${out} "${${out}}" PARENT_SCOPE)
endfunction()
  
#-------------------------------------------------------------------------------
# A helper to make lists of ffiles
#-------------------------------------------------------------------------------
function(append_file_list prefix num out)
  math(EXPR _num "${num}-1")
  
  foreach(_i RANGE ${_num})
    string(FIND ${prefix} "%" _first)
    string(FIND ${prefix} "%" _last  REVERSE)
    math(EXPR _last "${_last}+1")
    math(EXPR _size "${_last}-${_first}")
    if (_size LESS_EQUAL 0)
      MESSAGE(FATAL_ERROR "Could not find pattern character '%' to replace")
    endif()
    string(REPEAT "%" ${_size} _orig)
    string(LENGTH "${_i}" _ilen)
    math(EXPR _istart "${_size}-${_ilen}")
    if (_istart LESS 0)
      set(_istart "0")
    endif()
    string(REPEAT "0" ${_istart} _new)
    string(CONCAT _new ${_new} ${_i})
    string(REPLACE ${_orig} ${_new} _name ${prefix})
    list(APPEND ${out} ${_name})
  endforeach()
  set(${out} "${${out}}" PARENT_SCOPE)
endfunction()

#-------------------------------------------------------------------------------
# This macro creates a regression test

function(create_test)

  # the command to run to compare outputs
  set (TEST_COMMAND "${Python_EXECUTABLE} ${PRL_TOOL_DIR}/numdiff.py --check-text --absolute ${PRL_TEST_TOLERANCE}")

  # parse the arguments
  set(options)
  set(oneValueArgs NAME THREADS WORKING_DIRECTORY)
  set(multiValueArgs COMPARE STANDARD COMMAND INPUTS)
  cmake_parse_arguments(args "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  
  # check the preconditions
  if( NOT args_NAME )
    message( FATAL_ERROR "You must specify a test name using NAME." )
  endif()
  
  if( NOT args_COMMAND )
    message( FATAL_ERROR "You must specify a test command using COMMAND." )
  endif()
  
  list(LENGTH args_COMPARE  _compare_len)
  list(LENGTH args_STANDARD _standard_len)
  if ( NOT _compare_len EQUAL _standard_len )
    message( FATAL_ERROR "The number of files to compare with does not "
      "match the number of standards.\nLength of STANDARD must match that "
      "of COMPARE" )
  endif()

  # use at least one thread
  if( NOT args_THREADS )
    set( args_THREADS 1 )
  endif()

  if (args_WORKING_DIRECTORY)
    file(MAKE_DIRECTORY ${args_WORKING_DIRECTORY})
  endif()

  
  # add the test
  add_test( 
    NAME ${args_NAME}
    COMMAND ${CMAKE_COMMAND}
      "-Dtest_name=${args_NAME}"
      "-Dtest_cmd=${args_COMMAND}"
      -Dcompare_cmd=${TEST_COMMAND}
      "-Doutput_blessed=${args_STANDARD}"
      "-Doutput_test=${args_COMPARE}"
      -P ${REGRESSION_CMAKE_DIR}/run_test.cmake
    WORKING_DIRECTORY ${args_WORKING_DIRECTORY}
  )
  
  # for openmp
  SET_TESTS_PROPERTIES( ${args_NAME}
    PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=${args_THREADS}")

endfunction()

#-------------------------------------------------------------------------------
# This macro creates a regression test
#-------------------------------------------------------------------------------
function(create_regression)

  # parse the arguments
  set(options)
  set(oneValueArgs NAME)
  set(multiValueArgs ARGS PROCS COMPARE STANDARD PARTS)
  cmake_parse_arguments(args "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  
  # check the preconditions
  if( NOT args_NAME )
    message( FATAL_ERROR "You must specify a regression name using NAME." )
  endif()
  
  if( NOT args_ARGS )
    message( FATAL_ERROR "You must specify arguments." )
  endif()
  
  string(REPLACE ";" " " _args "${args_ARGS}")
  
  list(LENGTH args_COMPARE  _compare_len)
  list(LENGTH args_STANDARD _standard_len)
  if ( NOT _compare_len EQUAL _standard_len )
    message( FATAL_ERROR "The number of files to compare with does not "
      "match the number of standards.\nLength of STANDARD must match that "
      "of COMPARE" )
  endif()
  
  if( NOT args_PROCS )
    create_test(
      NAME ${args_NAME}
      COMMAND
        $<TARGET_FILE:prl> -i ${args_INPUT} ${_part_args}
      COMPARE ${args_COMPARE}
      STANDARD ${args_STANDARD})
  else()
    foreach(_procs ${args_PROCS})
      create_test(
        NAME ${args_NAME}-p${_procs}
        COMMAND
          ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${_procs}
          ${MPIEXEC_PREFLAGS} $<TARGET_FILE:prl>
          ${_args}
          ${MPIEXEC_POSTFLAGS}
        COMPARE ${args_COMPARE}
        STANDARD ${args_STANDARD}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/p${_procs})
    endforeach()
  endif()

endfunction()
 
 
