add_executable(graph_gen graph_gen.cpp
        ../src/common.cpp
        ../src/matds.cpp)

target_link_libraries(graph_gen PRIVATE mmio)
target_link_libraries(graph_gen PRIVATE zf_log)
