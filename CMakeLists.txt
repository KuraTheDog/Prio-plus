cmake_minimum_required(VERSION 3.5.1)

project(AggrProject C CXX)
set(CMAKE_CXX_STANDARD 11)

# Needed to get openmp to work
if(APPLE)
  # set(CMAKE_C_COMPILER "/opt/homebrew/Cellar/gcc/12.2.0/bin/g-12")
  # set(CMAKE_CXX_COMPILER "/opt/homebrew/Cellar/gcc/12.2.0/bin/g++-12")
  # set(CMAKE_C_COMPILER "/usr/local/bin/gcc-11")
  # set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-11")
  # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor")
  # include_directories("/usr/local/include" "/usr/local/opt/llvm/include")
  # link_directories("/usr/local/lib" "/usr/local/opt/llvm/lib")
endif()
link_directories("/usr/local/lib")

# if(NOT MSVC)
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
# else()
#   add_definitions(-D_WIN32_WINNT=0x600)
# endif()

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
include(${CMAKE_FOLDER}/cmake/emp-sh2pc-config.cmake)
include_directories(${EMP-SH2PC_INCLUDE_DIRS})

# Palisade
# find_package(Palisade)
# set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PALISADE_CXX_FLAGS}" )
# message("PALISADE_INCLUDE: ${PALISADE_INCLUDE}")
# include_directories( ${OPENMP_INCLUDES} )
# include_directories( ${PALISADE_INCLUDE} )
# include_directories( ${PALISADE_INCLUDE}/third-party/include )
# include_directories( ${PALISADE_INCLUDE}/core )
# include_directories( ${PALISADE_INCLUDE}/pke )
# link_directories( ${PALISADE_LIBDIR} )
# link_directories( ${OPENMP_LIBRARIES} )
# # Non-static
# set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${PALISADE_EXE_LINKER_FLAGS}" )
# link_libraries( ${PALISADE_SHARED_LIBRARIES} )

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
                 "constants.cpp" "ot.cpp" "fmpz_utils.cpp" "share.cpp" "net_share.cpp"
                 "correlated.cpp" "hash.cpp" "eval_heavy.cpp"
                 # "heavy_triples.cpp"
                 "poly/fft.c" "poly/poly_once.c" "poly/poly_batch.c"
                 )
  target_link_libraries(${_target}
    ${OPENSSL_LIBRARIES}
    ${Boost_LIBRARIES}
    ${GMP_LIBRARIES}
    ${FLINT_LIBRARIES}
    ${EMP-SH2PC_LIBRARIES}
  )
endforeach()

set(test_poly "test_circuit" "test_linreg")
set(test_correlated "test_ot" "test_bits" "test_heavy")
set(test_garble "test_garble")
set(test_hash "test_hash" ${test_garble})
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
  # test_he_triples
  test_heavy
  test_garble
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
    # set (test_SOURCE_FILES ${test_SOURCE_FILES} "he_triples.cpp")
  endif()
  if (_target IN_LIST test_garble)
    set (test_SOURCE_FILES ${test_SOURCE_FILES} "eval_heavy.cpp")
  endif()
  list(REMOVE_DUPLICATES test_SOURCE_FILES)
  # message(STATUS "${_target}: ${test_SOURCE_FILES}")
  add_executable(${_target} ${test_SOURCE_FILES})
  target_link_libraries(${_target}
    ${OPENSSL_LIBRARIES}
    ${Boost_LIBRARIES}
    ${GMP_LIBRARIES}
    ${FLINT_LIBRARIES}
    ${EMP-SH2PC_LIBRARIES}
  )
endforeach()
