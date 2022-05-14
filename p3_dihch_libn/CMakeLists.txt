CMAKE_MINIMUM_REQUIRED(VERSION 3.14)
PROJECT(project3)
set(CMAKE_C_FLAGS_DEBUG "-O0 -fopenmp -pthread -mavx2 -mavx")
set(CMAKE_C_FLAGS_RELEASE  "-O2 -fopenmp -pthread -mavx2 -mavx")

link_libraries(m)

add_executable(fast gbfloat_fast.c stb_image.h stb_image_write.h)
add_executable(tmp gbfloat_tmp.c stb_image.h stb_image_write.h)
add_executable(base gbfloat_tmp.c stb_image.h stb_image_write.h)
add_executable(check_base_fast test_accuracy.c)
add_executable(check_base_tmp test_accuracy.c)
add_executable(check_tmp_fast test_accuracy.c)
