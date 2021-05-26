add_executable(graph_gen_boost graph_gen_boost.cpp
        ../src/common.cpp
        ../src/matio.cpp
        ../src/matstorage.cpp)

target_link_libraries(graph_gen_boost PRIVATE mmio)
target_link_libraries(graph_gen_boost PRIVATE zf_log)
