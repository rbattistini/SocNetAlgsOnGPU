enable_testing()

add_executable(test_cc test_cc.cpp
        ../src/common.cpp
        ../src/spmatops.cpp
        ../src/ecc.cpp
        ../src/matstorage.cpp
        ../src/graphs.cpp)

target_link_libraries(test_cc PRIVATE zf_log)

add_test(NAME test_cc COMMAND test_cc)

add_executable(test_serial_bc test_serial_bc.cpp
        ../src/common.cpp
        ../src/matstorage.cpp
        ../src/bc.cpp)

target_link_libraries(test_serial_bc PRIVATE zf_log)

add_test(NAME test_serial_bc  COMMAND test_serial_bc)

add_executable(test_mod_bc test_mod_bc.cpp
        ../src/common.cpp
        ../src/matstorage.cpp
        ../src/graphs.cpp
        ../src/ecc.cpp
        ../src/spmatops.cpp
        ../src/spvb.cpp)

target_link_libraries(test_mod_bc PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(test_mod_bc PRIVATE zf_log)
target_link_libraries(test_mod_bc PRIVATE libsnap.a)

add_test(NAME test_mod_bc COMMAND test_mod_bc)

add_executable(test_matrix_io test_matrix_io.cpp
        ../src/common.cpp
        ../src/matstorage.cpp
        ../src/matio.cpp)

target_link_libraries(test_matrix_io PRIVATE mmio)
target_link_libraries(test_matrix_io PRIVATE zf_log)

add_test(NAME test_matrix_io COMMAND test_matrix_io)

add_executable(test_spmatops test_spmatops.cpp
        ../src/common.cpp
        ../src/spmatops.cpp
        ../src/matstorage.cpp)

target_link_libraries(test_spmatops PRIVATE zf_log)

add_test(NAME test_spmatops COMMAND test_spmatops)
