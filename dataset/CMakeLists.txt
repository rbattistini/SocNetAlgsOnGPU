add_executable(graph_gen graph_gen.cpp
        ../src/common.cpp
        ../src/matio.cpp
        ../src/matstorage.cpp)

target_link_libraries(graph_gen PRIVATE mmio)
target_link_libraries(graph_gen PRIVATE zf_log)

if(OpenMP_CXX_FOUND)
    target_link_libraries(graph_gen PUBLIC OpenMP::OpenMP_CXX libsnap.a)
endif()
