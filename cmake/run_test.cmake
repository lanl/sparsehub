# some argument checking:

# test_cmd is the command to run with all its arguments
if( NOT test_name )
   message( FATAL_ERROR "Variable test_name not defined" )
endif()

# test_cmd is the command to run with all its arguments
if( NOT test_cmd )
   message( FATAL_ERROR "Variable test_cmd not defined" )
endif()

# compare_cmd is the command to compare files
if( NOT compare_cmd )
   message( FATAL_ERROR "Variable compare_cmd not defined" )
endif()
    
list(LENGTH args_test  _test_len)
list(LENGTH args_blessed _blessed_len)
if ( NOT _blessed_len EQUAL _test_len )
  message( FATAL_ERROR "The number of files to compare with does not "
    "match the number of standards.\nLength of STANDARD must match that "
    "of COMPARE" )
endif()

# how many threads are there
message(STATUS "Using $ENV{OMP_NUM_THREADS} threads")

# blow away the compare-to-file in case it is already there
foreach(_file IN LISTS output_test )
  if (NOT _file STREQUAL "stdout")
    file(REMOVE ${_file})
  endif()
endforeach()

# run the test
separate_arguments( test_cmd ) 

string(REPLACE ";" " " test_cmd_string "${test_cmd}")
message(STATUS "Executing '${test_cmd_string}'")

set(output_file ${test_name}.txt)
execute_process(
   COMMAND ${test_cmd}
   OUTPUT_FILE ${output_file}
   ERROR_FILE ${output_file}
   RESULT_VARIABLE test_not_successful
)

if( test_not_successful )
   message( SEND_ERROR "Error running ${test_cmd}" )
endif()

# dump output
file(READ ${output_file} output)
if (output)
  MESSAGE(${output})
endif()

# need to fix the spaces in the passed command for some reason
separate_arguments( compare_cmd ) 

# run the diff
list(LENGTH output_blessed _blessed_len)
list(LENGTH output_test _test_len)

math(EXPR _iter_max "${_blessed_len}-1")
foreach( _iter RANGE ${_iter_max} )

  list(GET output_blessed ${_iter} _blessed )
  list(GET output_test ${_iter} _test )
  MESSAGE( STATUS "Checking file ${_blessed} against ${_test}" )

  if (_test STREQUAL "stdout")
    set(_test ${output_file})
  endif()

  if(DEFINED ENV{BUILD_STANDARDS})
    set(new_std $ENV{BUILD_STANDARDS})
  else()
    file(TOUCH ${_test})
  endif()

  # NEW STANDARDS
  if(new_std)

    MESSAGE(STATUS "Rebuilding standards.")
    get_filename_component(_dirname ${_blessed} DIRECTORY)
    file(MAKE_DIRECTORY ${_dirname})
    message(STATUS "${_test} vs ${_blessed}")
    configure_file(${_test} ${_blessed} COPYONLY)

  # COMPARE
  else()

    string(REPLACE ";" " " test_cmd_string "${compare_cmd} ${_blessed} ${_test}")
    message(STATUS "Executing '${test_cmd_string}'")

    execute_process(
      COMMAND ${compare_cmd} ${_blessed} ${_test}
      RESULT_VARIABLE test_not_successful
    )

    if( test_not_successful )
      message( SEND_ERROR "${_test} does not match ${_blessed}!" )
    endif()

  endif()

endforeach()
