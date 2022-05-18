CMAKE_MINIMUM_REQUIRED(VERSION 3.14)
PROJECT(project3)
set(CMAKE_C_FLAGS_DEBUG "-g -O0 -fopenmp -pthread -mavx2 -mavx -mfma")
set(CMAKE_C_FLAGS_RELEASE  "-O2 -fopenmp -pthread -mavx2 -mavx -mfma")

link_libraries(m)

add_executable(fast gbfloat_fast.c stb_image.h stb_image_write.h)
add_executable(tmp gbfloat_m128.c stb_image.h stb_image_write.h)
add_executable(base gbfloat_base.c stb_image.h stb_image_write.h)
add_executable(check_base_fast test_accuracy.c)
add_executable(check_base_tmp test_accuracy.c)
add_executable(check_tmp_fast test_accuracy.c)

