cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(SWIG REQUIRED)
find_package(Boost REQUIRED)

set(GRAPHEXT_WRAP ${CMAKE_CURRENT_BINARY_DIR}/clipGraphExtractor_wrap.cc)
set(GRAPHEXT_TCL_INIT ${CMAKE_CURRENT_BINARY_DIR}/clipGraphExtTclInitVar.cc)

add_custom_command(OUTPUT ${GRAPHEXT_WRAP}
  COMMAND ${SWIG_EXECUTABLE} -tcl8 -c++ -I${CLIP_GRAPHEXT_HOME}/include -I${OPENDB_HOME}/src/swig/tcl -o ${GRAPHEXT_WRAP} ${CLIP_GRAPHEXT_HOME}/src/clipGraphExtractor.i
  WORKING_DIRECTORY ${CLIP_GRAPHEXT_HOME}/src
  DEPENDS ${CLIP_GRAPHEXT_HOME}/src/clipGraphExtractor.i ${CLIP_GRAPHEXT_HOME}/include/clip_graph_ext/clipGraphExtractor.h
)

add_custom_command(OUTPUT ${GRAPHEXT_TCL_INIT}
  COMMAND ${OPENSTA_HOME}/etc/TclEncode.tcl ${GRAPHEXT_TCL_INIT} graph_extractor_tcl_inits ${CLIP_GRAPHEXT_HOME}/src/clipGraphExtractor.tcl
  WORKING_DIRECTORY ${CLIP_GRAPHEXT_HOME}/src
  DEPENDS ${CLIP_GRAPHEXT_HOME}/src/clipGraphExtractor.tcl ${OPENSTA_HOME}/etc/TclEncode.tcl
  )

############################################################
# Libraries setting 
############################################################

add_library( ClipGraphExtractor
  src/clipGraphExtractor.cpp
  src/MakeClipGraphExtractor.cpp
  src/instGraph.cpp
  src/binGraph.cpp

  ${GRAPHEXT_WRAP}
  ${GRAPHEXT_TCL_INIT}
)

target_include_directories( ClipGraphExtractor
  PUBLIC ${CLIP_GRAPHEXT_HOME}/include
  PRIVATE
  ${CLIP_GRAPHEXT_HOME}/src
  ${OPENROAD_HOME}/include
  ${OPENDB_HOME}/include
  ${DBSTA_HOME}/include
  ${OPENSTA_HOME}/include

  ${Boost_INCLUDE_DIR}
 )
