cmake_minimum_required(VERSION 2.6)
set(BOOST_ROOT "/home/micha/")
set(Boost_DEBUG ON)
project(idresnet)
find_package(Boost COMPONENTS program_options)
MESSAGE("${BOOST_VERSION}")

include_directories("/usr/include/suitesparse")

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    add_executable(idresnet.bin main.cpp)
    add_executable(exp_cc.bin main_topo.cpp)
    target_link_libraries(idresnet.bin ${Boost_LIBRARIES})                                                                                                                                                                                                                            
    target_link_libraries(exp_cc.bin ${Boost_LIBRARIES})                                                                                                                                                                                                                            
endif()

target_link_libraries(idresnet.bin cholmod)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Ofast -fopenmp -march=native -Wall -Wno-unused-local-typedefs")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=c++11 -O0 -march=native -Wall -Wno-unused-local-typedefs")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g3 -Wall -Wno-unused-local-typedefs")
install(TARGETS idresnet.bin RUNTIME DESTINATION bin)
install(TARGETS exp_cc.bin RUNTIME DESTINATION bin)
