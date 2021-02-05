cmake_minimum_required(VERSION 3.5.1)

project(AggrProject C CXX)

if(NOT MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
  add_definitions(-D_WIN32_WINNT=0x600)
endif()

# Find flint
if(FLINT_INCLUDE_DIR AND FLINT_LIBRARIES)
    # Already in cache, be silent
    set(FLINT_FIND_QUIETLY TRUE)
endif()
find_path(FLINT_INCLUDE_DIR flint/flint.h)
find_library(FLINT_LIBRARIES NAMES flint)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FLINT DEFAULT_MSG FLINT_INCLUDE_DIR FLINT_LIBRARIES)
mark_as_advanced(FLINT_INCLUDE_DIR FLINT_LIBRARIES)

# Palisade
find_package(Palisade)
set( CMAKE_CXX_FLAGS ${PALISADE_CXX_FLAGS} )
set( CMAKE_EXE_LINKER_FLAGS ${PALISADE_EXE_LINKER_FLAGS} )
include_directories( ${OPENMP_INCLUDES} )
include_directories( ${PALISADE_INCLUDE} )
include_directories( ${PALISADE_INCLUDE}/third-party/include )
include_directories( ${PALISADE_INCLUDE}/core )
include_directories( ${PALISADE_INCLUDE}/pke )
### add directories for other PALISADE modules as needed for your project
link_directories( ${PALISADE_LIBDIR} )
link_directories( ${OPENMP_LIBRARIES} )
link_libraries( ${PALISADE_LIBRARIES} )

include(cmake/common.cmake)
find_package(EMP-TOOL REQUIRED)
include_directories(${EMP-TOOL_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

foreach(_target
  server client
)
  add_executable(${_target} "${_target}.cpp" 
                 "constants.cpp" "proto.cpp" "fmpz_utils.cpp" "share.cpp" "net_share.cpp" "edabit.cpp" "he_triples.cpp"
                 "poly/fft.c" "poly/poly_once.c" "poly/poly_batch.c"
                 )
  target_link_libraries(${_target}
    ${OPENSSL_LIBRARIES}
    ${Boost_LIBRARIES}
    ${GMP_LIBRARIES}
    ${EMP-TOOL_LIBRARIES}
    ${FLINT_LIBRARIES}
  )
endforeach()

foreach(_target
  test_net_share
  test_circuit
  test_linreg
  test_proto
  test_bits
  test_he_triples
)
  add_executable(${_target} "test/${_target}.cpp" 
                 "constants.cpp" "proto.cpp" "fmpz_utils.cpp" "share.cpp" "net_share.cpp" "edabit.cpp" "he_triples.cpp"
                 "poly/fft.c" "poly/poly_once.c" "poly/poly_batch.c"
                 )
  target_link_libraries(${_target}
    ${OPENSSL_LIBRARIES}
    ${Boost_LIBRARIES}
    ${GMP_LIBRARIES}
    ${EMP-TOOL_LIBRARIES}
    ${FLINT_LIBRARIES}
  )
endforeach()
