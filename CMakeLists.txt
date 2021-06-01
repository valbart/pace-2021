cmake_minimum_required(VERSION 3.0)

project(mainoptil)

include_directories(${CMAKE_SOURCE_DIR}/exact)
include_directories(${CMAKE_SOURCE_DIR}/shared)
include_directories(${CMAKE_SOURCE_DIR}/heuristic/src)

link_directories(${CMAKE_SOURCE_DIR}/lib)


set(CMAKE_CXX_STANDARD 17)

set(SRCS
     exact/main_optil.cpp 
     exact/cluster_graph.cpp 
     exact/instance.cpp 
     shared/solution.cpp 
     shared/graph.cpp 
     heuristic/src/low_mem_graph.cpp 
     heuristic/src/local_search.cpp  
     heuristic/src/simulated_annealing.cpp 
     heuristic/src/instance.cpp 
     heuristic/src/multi_cc_instance.cpp 
     heuristic/src/profiler.cpp 
     heuristic/src/union_find.cpp 
     exact/lp_solve.cpp
    )
    
set(HEADERS
     exact/cluster_graph.hpp 
     exact/instance.hpp
     shared/solution.hpp 
     shared/graph.hpp 
     heuristic/src/low_mem_graph.h 
     heuristic/src/local_search.h  
     heuristic/src/simulated_annealing.h 
     heuristic/src/instance.h 
     heuristic/src/multi_cc_instance.h 
     heuristic/src/profiler.h 
     heuristic/src/union_find.h 
     exact/lp_solve.hpp
    )

add_executable(mainoptil ${SRCS} ${HEADERS} )

target_link_libraries(mainoptil Clp ClpSolver CoinUtils z)

target_compile_options(mainoptil PUBLIC -Ofast -static)

