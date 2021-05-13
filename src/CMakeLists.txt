target_include_directories( sparsehub PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> )

#------------------------------------------------------------------------------#
# configure header
#------------------------------------------------------------------------------#

configure_file(config.hpp.in config.hpp @ONLY)

target_sources( sparsehub PRIVATE  ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp )
