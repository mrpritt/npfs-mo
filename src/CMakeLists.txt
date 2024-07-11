cmake_minimum_required(VERSION 3.25)

project(npfsmo)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
message(STATUS "Build configuration: " ${CMAKE_BUILD_TYPE})

set(CMAKE_CXX_FLAGS "-g -ggdb -Wall -Wno-unknown-pragmas -fPIC -DIL_STD -fexceptions -std=c++17")

find_package(cplex)
set(Boost_NO_SYSTEM_PATHS ON)
find_package(Boost 1.75.0 COMPONENTS program_options timer system chrono filesystem)

include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
if (CPLEX_FOUND)
  include_directories(SYSTEM ${CPLEX_INCLUDE_DIRS})
endif()
include_directories(fmt/include)
link_directories(${Boost_LIBRARY_DIRS})

set(std_sources options.cpp random.cpp logging.cpp instance.cpp helpers.cpp)

add_executable(npfsmo npfsmo.cpp instance.cpp solution.cpp heuristics.cpp holes.cpp ${std_sources})
target_link_libraries(npfsmo ${CPLEX_LIBRARIES} ${Boost_LIBRARIES} rt stdc++fs dl)

if (CPLEX_FOUND)
  add_executable(exact exact.cpp instance.cpp solution.cpp holes.cpp models.cpp heuristics.cpp ${std_sources})
  target_link_libraries(exact ${CPLEX_LIBRARIES} ${Boost_LIBRARIES} rt stdc++fs dl)
else()
  message(STATUS "CPLEX model won't be built.")
endif()
