cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(LatticeYM)

include_directories(include)
file(GLOB sources src/*.cc)
file(GLOB header include/*.hh)

#To debug the lattice/simulation and SU3/utils sources, respectively
#comment out all those which you don't want to build

#add_executable(test2 test2.cc ${sources} ${headers})
#add_custom_target(LatticeYM DEPENDS test2)
#install(TARGETS test2 DESTINATION build)

#add_executable(test test.cc ${sources} ${headers})
#add_custom_target(LatticeYM DEPENDS test)
#install(TARGETS test DESTINATION build)

#Actual simulation code 

add_executable(latticeYM latticeYM.cc ${sources} ${headers})
add_custom_target(LatticeYM DEPENDS latticeYM)
install(TARGETS latticeYM DESTINATION build)
