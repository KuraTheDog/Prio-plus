cmake_minimum_required(VERSION 3.5.1)

project(AggrProject C CXX)

if(NOT MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
  add_definitions(-D_WIN32_WINNT=0x600)
endif()

# find_path(CMAKE_FOLDER NAMES cmake/emp-tool-config.cmake)
include(cmake/common.cmake)
# include(cmake/source_of_randomness.cmake)
# include(cmake/threading.cmake)
find_package(emp-tool REQUIRED)
include_directories(${EMP-TOOL_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

# Targets greeter_[async_](client|server)
foreach(_target
  server client)
  add_executable(${_target} "${_target}.cpp" "proto.cpp")
  target_link_libraries(${_target}
    ${OPENSSL_LIBRARIES}
    ${Boost_LIBRARIES}
    ${GMP_LIBRARIES}
    ${EMP-TOOL_LIBRARIES})
endforeach()