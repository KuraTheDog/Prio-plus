cmake_minimum_required(VERSION 3.5.1)

project(AggrProject C CXX)
set(CMAKE_CXX_STANDARD 11)

# Needed to get openmp to work
if(APPLE)
      set(CMAKE_C_COMPILER "/usr/local/bin/gcc-11")
      set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-11")
endif()

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

include(cmake/common.cmake)
find_package(EMP-TOOL REQUIRED)
include_directories(${EMP-TOOL_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

# Palisade
find_package(OpenMP REQUIRED)
find_package(Palisade)
include_directories( ${OPENMP_INCLUDES} )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PALISADE_CXX_FLAGS}" )
include_directories( ${PALISADE_INCLUDE} )
include_directories( ${PALISADE_INCLUDE}/third-party/include )
include_directories( ${PALISADE_INCLUDE}/core )
include_directories( ${PALISADE_INCLUDE}/pke )
link_directories( ${PALISADE_LIBDIR} )
link_directories( ${OPENMP_LIBRARIES} )
# Non-static
set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${PALISADE_EXE_LINKER_FLAGS}" )
link_libraries( ${PALISADE_SHARED_LIBRARIES} )

foreach(_target
  server client
)
  add_executable(${_target} "${_target}.cpp" 
                 "constants.cpp" "ot.cpp" "fmpz_utils.cpp" "share.cpp" "net_share.cpp" "correlated.cpp" "hash.cpp" "he_triples.cpp"
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

set(test_poly "test_circuit" "test_linreg")
set(test_correlated "test_ot" "test_bits" "test_heavy")
set(test_hash "test_hash")
set(test_he_triples "test_he_triples" ${test_correlated})
# stuff that sends shares
set(test_net_share "test_net_share" ${test_poly} ${test_correlated} ${test_he_triples})
set(test_share "test_share" ${test_net_share})
foreach(_target
  test_net_share
  test_circuit
  test_linreg
  test_ot
  test_bits
  test_hash
  test_he_triples
  test_heavy
)
  set (test_SOURCE_FILES "test/${_target}.cpp")
  set (test_SOURCE_FILES ${test_SOURCE_FILES} "constants.cpp" "fmpz_utils.cpp")
  if (_target IN_LIST test_poly)
    set (test_SOURCE_FILES ${test_SOURCE_FILES} 
         "poly/fft.c" "poly/poly_once.c" "poly/poly_batch.c")
  endif()
  if (_target IN_LIST test_share)
    set (test_SOURCE_FILES ${test_SOURCE_FILES} "share.cpp")
  endif()
  if (_target IN_LIST test_net_share)
    set (test_SOURCE_FILES ${test_SOURCE_FILES} "net_share.cpp")
  endif()
  if (_target IN_LIST test_correlated)
    set (test_SOURCE_FILES ${test_SOURCE_FILES} "correlated.cpp" "ot.cpp")
  endif()
  if (_target IN_LIST test_hash)
    set (test_SOURCE_FILES ${test_SOURCE_FILES} "hash.cpp")
  endif()
  if (_target IN_LIST test_he_triples)
    set (test_SOURCE_FILES ${test_SOURCE_FILES} "he_triples.cpp")
  endif()
  list(REMOVE_DUPLICATES test_SOURCE_FILES)
  # message(STATUS "${_target}: ${test_SOURCE_FILES}")
  add_executable(${_target} ${test_SOURCE_FILES})
  target_link_libraries(${_target}
    ${OPENSSL_LIBRARIES}
    ${Boost_LIBRARIES}
    ${GMP_LIBRARIES}
    ${EMP-TOOL_LIBRARIES}
    ${FLINT_LIBRARIES}
  )
endforeach()
