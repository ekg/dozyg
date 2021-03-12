# Specify the minimum version for CMake

cmake_minimum_required(VERSION 3.1)

# Project's name
project(dozyg)
# We build using c++14
set(CMAKE_CXX_STANDARD 14)

# Use all standard-compliant optimizations
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -O3 -mcx16 -msse4.1 -g")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -O3 -mcx16 -msse4.1 -g")

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  # assumes clang build
  # we can't reliably detect when we're using clang, so for the time being we assume
  # TODO: can't we though?
  # adapted from https://stackoverflow.com/questions/46414660/macos-cmake-and-openmp
  # find_package(OpenMP) does not work reliably on macOS, so we do its work ourselves
  set (OpenMP_C "${CMAKE_C_COMPILER}")
  set (OpenMP_C_FLAGS " -Xpreprocessor -fopenmp -I/opt/local/include/libomp -I/usr/local/include -L/opt/local/lib/libomp -L/usr/local/lib")
  set (OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
  set (OpenMP_CXX "${CMAKE_CXX_COMPILER}")
  set (OpenMP_CXX_FLAGS " -Xpreprocessor -fopenmp -I/opt/local/include/libomp -I/usr/local/include -L/opt/local/lib/libomp -L/usr/local/lib")
  set (OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
  set (OpenMP_libomp_LIBRARY "omp")
  set (OpenMP_libgomp_LIBRARY "gomp")
  set (OpenMP_libiomp5_LIBRARY "iomp5")
  # and now add the OpenMP parameters to the compile flags
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} -lomp")
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  find_package(OpenMP REQUIRED)
  # add the flags it detects to the compile flags
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fopenmp")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

# The following folder will be included
include_directories("${PROJECT_SOURCE_DIR}")

# Add external projects
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

# TODO: We're using INSTALL_DIR very wrong. We *should* be actually installing
# the external projects into their prefixes and working with the installed
# files. Instead we're building but not installing them and trying to work with
# the non-installed build trees.
# 
# Hence the blanked out INSTALL_COMMANDs to suppress the install step.
#
# We need to NOT blank out UPDATE_COMMAND or we can never change the Git revision we point to.
# The cost of this is that we have to re-configure on every build if we do update.

# sdsl-lite (full build using its cmake config)
ExternalProject_Add(sdsl-lite
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/sdsl-lite"
  CMAKE_ARGS "${CMAKE_ARGS};-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(sdsl-lite INSTALL_DIR)
set(sdsl-lite_INCLUDE "${INSTALL_DIR}/src/sdsl-lite-build/include")
set(sdsl-lite-divsufsort_INCLUDE "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/include")
set(sdsl-lite_LIB "${INSTALL_DIR}/src/sdsl-lite-build/lib")
set(sdsl-lite-divsufsort_LIB "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/lib")

# DYNAMIC (full build using its cmake config)
ExternalProject_Add(dynamic
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/DYNAMIC"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(dynamic SOURCE_DIR)
set(dynamic_INCLUDE "${SOURCE_DIR}/include")

# hopscotch_map (required by DYNAMIC)
ExternalProject_Add(hopscotch_map
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/hopscotch-map"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(hopscotch_map SOURCE_DIR)
set(hopscotch_map_INCLUDE "${SOURCE_DIR}/include")

# gfakluge (now header only)
ExternalProject_Add(gfakluge
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/gfakluge"
  CMAKE_ARGS "${CMAKE_ARGS};-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(gfakluge SOURCE_DIR)
set(gfakluge_INCLUDE "${SOURCE_DIR}/src")
set(gfakluge_tinyFA_INCLUDE "${SOURCE_DIR}/src/tinyFA")
ExternalProject_Get_property(gfakluge INSTALL_DIR)
set(gfakluge_LIB "${INSTALL_DIR}/src/gfakluge")

# libhandlegraph (full build using its cmake config)
ExternalProject_Add(handlegraph
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/libhandlegraph"
  CMAKE_ARGS "${CMAKE_ARGS};-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>")
ExternalProject_Get_property(handlegraph INSTALL_DIR)
set(handlegraph_INCLUDE "${INSTALL_DIR}/include")
set(handlegraph_LIB "${INSTALL_DIR}/lib")

# taywee's C++ args library, header only
ExternalProject_Add(tayweeargs
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/args"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(tayweeargs SOURCE_DIR)
set(tayweeargs_INCLUDE "${SOURCE_DIR}")

# BBHash perfect hasher
ExternalProject_Add(bbhash
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/BBHash"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(bbhash SOURCE_DIR)
set(bbhash_INCLUDE "${SOURCE_DIR}")

# ska
ExternalProject_Add(ska
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/flat_hash_map"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(ska SOURCE_DIR)
set(ska_INCLUDE "${SOURCE_DIR}")

# structures
ExternalProject_Add(structures
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/structures"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(structures SOURCE_DIR)
set(structures_INCLUDE "${SOURCE_DIR}/src/include")

ExternalProject_Add(picosha256
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/PicoSHA2"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(picosha256 SOURCE_DIR)
set(picosha256_INCLUDE "${SOURCE_DIR}")

# mmmultimap (memory mapped multimap)
ExternalProject_Add(mmmultimap
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/mmmultimap"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(mmmultimap SOURCE_DIR)
set(mmmultimap_INCLUDE "${SOURCE_DIR}/src")

# mmap allocator
ExternalProject_Add(mmap_allocator
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/mmap_allocator"
  CMAKE_ARGS "${CMAKE_ARGS};-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>"
  UPDATE_COMMAND "")
ExternalProject_Get_property(mmap_allocator INSTALL_DIR)
set(mmap_allocator_INCLUDE "${INSTALL_DIR}/include/mmap_allocator")
set(mmap_allocator_LIB "${INSTALL_DIR}/lib")

# In-place Parallel Super Scalar Samplesort (IPS⁴o), header only
ExternalProject_Add(ips4o
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/ips4o"
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(ips4o SOURCE_DIR)
set(ips4o_INCLUDE "${SOURCE_DIR}")

# gzip file reader, header only
ExternalProject_Add(gzipreader
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/gzip_reader"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(gzipreader SOURCE_DIR)
set(gzipreader_INCLUDE "${SOURCE_DIR}")

# gzip file reader, header only
ExternalProject_Add(cgranges
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/cgranges"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(cgranges SOURCE_DIR)
set(cgranges_INCLUDE "${SOURCE_DIR}/cpp")

# edlib
add_subdirectory(
  deps/edlib EXCLUDE_FROM_ALL)
set(edlib_INCLUDE deps/edlib/edlib/include)

# atomic queue library
ExternalProject_Add(atomicqueue
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/atomic_queue"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(atomicqueue SOURCE_DIR)
set(atomicqueue_INCLUDE "${SOURCE_DIR}")

# dozeu
ExternalProject_Add(dozeu
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/dozeu"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(dozeu SOURCE_DIR)
set(dozeu_INCLUDE "${SOURCE_DIR}")

#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_BUILD_TYPE Debug)

# set up our target executable and specify its dependencies and includes
add_library(dozyg_objs OBJECT
  ${CMAKE_SOURCE_DIR}/src/odgi.cpp
  ${CMAKE_SOURCE_DIR}/src/node.cpp # for odgi
#  ${CMAKE_SOURCE_DIR}/src/xg.cpp
  ${CMAKE_SOURCE_DIR}/src/index.cpp
  ${CMAKE_SOURCE_DIR}/src/chain.cpp
  ${CMAKE_SOURCE_DIR}/src/align.cpp
  ${CMAKE_SOURCE_DIR}/src/mapper.cpp
  ${CMAKE_SOURCE_DIR}/src/position.cpp
  ${CMAKE_SOURCE_DIR}/src/threads.cpp
  ${CMAKE_SOURCE_DIR}/src/utility.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/subcommand.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/index_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/query_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/map_main.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/dfs.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/kmer.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/hash.cpp
  )

set(dozyg_DEPS 
    sdsl-lite
    dynamic
    hopscotch_map
    gfakluge
    handlegraph
    tayweeargs
    bbhash
    ska
    structures
    picosha256
    gzipreader
    mmmultimap
    mmap_allocator
    ips4o
    cgranges
    edlib
    atomicqueue)
add_dependencies(dozyg_objs ${dozyg_DEPS})

set(dozyg_INCLUDES
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/algorithms
  ${sdsl-lite_INCLUDE}
  ${sdsl-lite-divsufsort_INCLUDE}
  ${dynamic_INCLUDE}
  ${hopscotch_map_INCLUDE}
  ${gfakluge_INCLUDE}
  ${gfakluge_tinyFA_INCLUDE}
  ${handlegraph_INCLUDE}
  ${tayweeargs_INCLUDE}
  ${ska_INCLUDE}
  ${bbhash_INCLUDE}
  ${picosha256_INCLUDE}
  ${gzipreader_INCLUDE}
  ${mmmultimap_INCLUDE}
  ${mmap_allocator_INCLUDE}
  ${ips4o_INCLUDE}
  ${cgranges_INCLUDE}
  ${edlib_INCLUDE}
  ${atomicqueue_INCLUDE}
  ${dozeu_INCLUDE})

set(dozyg_LIBS
  ${sdsl-lite_LIB}/libsdsl.a
  ${sdsl-lite-divsufsort_LIB}/libdivsufsort.a
  ${sdsl-lite-divsufsort_LIB}/libdivsufsort64.a
  ${handlegraph_LIB}/libhandlegraph.a
  ${mmap_allocator_LIB}/libmmap_allocator.a
  edlib
#  dl
  atomic
  z)

#set(dozyg_HEADERS
#  ${CMAKE_SOURCE_DIR}/src/xg.hpp
#  ${CMAKE_SOURCE_DIR}/src/dozyg.hpp
#  ${CMAKE_SOURCE_DIR}/src/subcommand/subcommand.hpp
#  ${CMAKE_SOURCE_DIR}/src/algorithms/remove_high_degree.hpp
#  ${CMAKE_SOURCE_DIR}/src/algorithms/hash.hpp
#  ${CMAKE_SOURCE_DIR}/src/algorithms/dfs.hpp
#  ${CMAKE_SOURCE_DIR}/src/algorithms/kmer.hpp
#  ${CMAKE_SOURCE_DIR}/src/algorithms/prune.hpp)

target_include_directories(dozyg_objs PUBLIC ${dozyg_INCLUDES})
set_target_properties(dozyg_objs PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

#add_library(libdozyg_static STATIC $<TARGET_OBJECTS:dozyg_objs>)
#set_target_properties(libdozyg_static PROPERTIES OUTPUT_NAME "dozyg")
#set_target_properties(libdozyg_static PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
#add_library(libdozyg_shared SHARED $<TARGET_OBJECTS:dozyg_objs>)
#set_target_properties(libdozyg_shared PROPERTIES OUTPUT_NAME "dozyg")
#set_target_properties(libdozyg_shared PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

#set_target_properties(libdozyg_static PROPERTIES PUBLIC_HEADER "${dozyg_HEADERS}")
#set_target_properties(libdozyg_shared PROPERTIES PUBLIC_HEADER "${dozyg_HEADERS}")

add_executable(dozyg
  $<TARGET_OBJECTS:dozyg_objs>
  ${CMAKE_SOURCE_DIR}/src/main.cpp)
add_dependencies(dozyg dozyg_objs)
target_link_libraries(dozyg ${dozyg_LIBS})
set_target_properties(dozyg PROPERTIES OUTPUT_NAME "dozyg")


#if (NOT EXISTS ${CMAKE_SOURCE_DIR}/deps/pybind11/CMakeLists.txt)
#  execute_process(COMMAND git submodule update --init --recursive)
#endif()

#add_subdirectory(deps/pybind11)
#pybind11_add_module(dozyg_pybind11 "${CMAKE_SOURCE_DIR}/src/pythonmodule.cpp")
#add_dependencies(dozyg_pybind11 ${dozyg_DEPS} libdozyg_static)
#target_include_directories(dozyg_pybind11 PUBLIC ${dozyg_INCLUDES})
#target_link_libraries(dozyg_pybind11 PRIVATE "${CMAKE_SOURCE_DIR}/lib/libdozyg.a" "${dozyg_LIBS}")
#set_target_properties(dozyg_pybind11 PROPERTIES OUTPUT_NAME "dozyg")

install(TARGETS dozyg DESTINATION bin)
#install(TARGETS libdozyg_static ARCHIVE DESTINATION lib LIBRARY DESTINATION lib PUBLIC_HEADER DESTINATION include/dozyg)
#install(TARGETS libdozyg_shared ARCHIVE DESTINATION lib LIBRARY DESTINATION lib PUBLIC_HEADER DESTINATION include/dozyg)
#install(TARGETS dozyg_pybind11 LIBRARY DESTINATION lib)

if (APPLE)
elseif (TRUE)
  if (BUILD_STATIC)
    set(CMAKE_EXE_LINKER_FLAGS "-static")
  endif()
endif()
