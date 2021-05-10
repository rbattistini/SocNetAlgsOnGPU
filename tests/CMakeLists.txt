#set(CMAKE_CXX_FLAGS "-Wall -Werror -Wpedantic")

add_executable(test_cc test_cc.cpp
        ../include/utils.h ../src/utils.cpp
        ../include/matstorage.h ../src/matstorage.cpp
        ../include/graphs.h ../src/graphs.cpp)

add_executable(test_serial_bc test_serial_bc.cpp
        ../include/utils.h ../src/utils.cpp
        ../include/matstorage.h ../src/matstorage.cpp
        ../include/graphs.h ../src/graphs.cpp)

add_executable(test_mod_bc test_mod_bc.cpp
        ../include/utils.h ../src/utils.cpp
        ../include/matstorage.h ../src/matstorage.cpp
        ../include/graphs.h ../src/graphs.cpp
        ../include/spvb.h ../src/spvb.cpp)

target_link_libraries(test_mod_bc PUBLIC OpenMP::OpenMP_CXX libsnap.a)

add_executable(test_matrix_io test_matrix_io.cpp
        ../include/utils.h ../src/utils.cpp
        ../include/matio.h ../src/matio.cpp
        ../include/matstorage.h ../src/matstorage.cpp
        ../include/graphs.h ../src/graphs.cpp)

target_link_libraries(test_matrix_io PRIVATE mmio)
