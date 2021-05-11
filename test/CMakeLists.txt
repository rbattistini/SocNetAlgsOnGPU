add_executable(test_cc test_cc.cpp
        ../src/utils.cpp
        ../src/matstorage.cpp
        ../src/graphs.cpp)

add_executable(test_serial_bc test_serial_bc.cpp
        ../src/utils.cpp
        ../src/matstorage.cpp
        ../src/graphs.cpp)

add_executable(test_mod_bc test_mod_bc.cpp
        ../src/utils.cpp
        ../src/matstorage.cpp
        ../src/graphs.cpp
        ../src/spvb.cpp)

target_link_libraries(test_mod_bc PRIVATE OpenMP::OpenMP_CXX libsnap.a)

add_executable(test_matrix_io test_matrix_io.cpp
        ../src/utils.cpp
        ../src/matstorage.cpp
        ../src/graphs.cpp
        ../src/matio.cpp)

target_link_libraries(test_matrix_io PRIVATE mmio)
