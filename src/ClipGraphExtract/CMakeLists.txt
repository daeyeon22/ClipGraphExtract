cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)



option(USE_CIMG_LIB "Use CImg drawing library for image drawing" ON)


if( USE_CIMG_LIB )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Dcimg_use_jpeg=1 -D_GLIBCXX_USE_CXX11_ABI=0")
    #set(CIMG_HOME ${CLIP_GRAPHEXT_HOME}/module/CImg)
  set(CIMG_HOME ${PROJECT_SOURCE_DIR}/src/CImg)
  set(COLORMAP_HOME ${PROJECT_SOURCE_DIR}/src/tinycolormap)
    add_definitions(-DENABLE_CIMG_LIB=1)
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")
endif()

   
    
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

if( USE_CIMG_LIB )
  find_package(X11 REQUIRED)
  find_package(JPEG REQUIRED)
  set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${X11_LIBRARIES} ${JPEG_LIBRARIES}")
else()
endif()

add_library( ClipGraphExtractor
  src/clipGraphExtractor.cpp
  src/MakeClipGraphExtractor.cpp
#src/instGraph.cpp
#src/vertex.cpp
  src/grid.cpp
  src/gcell.cpp
  src/rsmt.cpp
  src/plot.cpp
  src/marker.cpp
  src/label.cpp
  src/writer.cpp  
  src/feaExt.cpp
  src/mymeasure.cpp
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
  ${FLUTE_HOME}
    #${CIMG_HOME}
  ${Boost_INCLUDE_DIR}
 )


# CIMG_USE case
if( USE_CIMG_LIB )
  target_link_libraries( ClipGraphExtractor
    PRIVATE
    ${X11_LIBRARIES} 
    ${JPEG_LIBRARIES} 
  )
  target_include_directories( ClipGraphExtractor
    PUBLIC
    # For CImg
    ${CIMG_HOME}
    ${COLORMAP_HOME}/include
    ${X11_INCLUDE_DIR}
    ${JPEG_INCLUDE_DIR}
  )
endif()
