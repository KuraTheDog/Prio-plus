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

include(cmake/FindGMP.cmake)

find_path(CMAKE_FOLDER NAMES cmake/emp-base.cmake)
message(STATUS "CMAKE_FOLDER: ${CMAKE_FOLDER}")
include(${CMAKE_FOLDER}/cmake/emp-base.cmake)
include(${CMAKE_FOLDER}/cmake/emp-tool-config.cmake)
include(${CMAKE_FOLDER}/cmake/emp-ot-config.cmake)
include_directories(${EMP-OT_INCLUDE_DIRS})

# # For debugging
set(VALGRIND_TEST FALSE)
if(VALGRIND_TEST)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g")
  message(STATUS "DEBUG: Running with debug flags for valgrind testing")
endif()


foreach(_target
  server client
)
  add_executable(${_target} "${_target}.cpp" 
                 "constants.cpp" "ot.cpp" "fmpz_utils.cpp" "share.cpp" "net_share.cpp" "correlated.cpp" "hash.cpp"
                 "poly/fft.c" "poly/poly_once.c" "poly/poly_batch.c"
                 )
  target_link_libraries(${_target}
    ${OPENSSL_LIBRARIES}
    ${Boost_LIBRARIES}
    ${GMP_LIBRARIES}
    ${FLINT_LIBRARIES}
  )
endforeach()

set(test_poly "test_circuit" "test_linreg")
set(test_correlated "test_ot" "test_bits")
set(test_hash "test_hash")
# stuff that sends shares
set(test_net_share "test_net_share" ${test_poly} ${test_correlated})
set(test_share "test_share" ${test_net_share})
foreach(_target
  test_net_share
  test_circuit
  test_linreg
  test_ot
  test_bits
  test_hash
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
  list(REMOVE_DUPLICATES test_SOURCE_FILES)
  # message(STATUS "${_target}: ${test_SOURCE_FILES}")
  add_executable(${_target} ${test_SOURCE_FILES})
  target_link_libraries(${_target}
    ${OPENSSL_LIBRARIES}
    ${Boost_LIBRARIES}
    ${GMP_LIBRARIES}
    ${FLINT_LIBRARIES}
  )
endforeach()
