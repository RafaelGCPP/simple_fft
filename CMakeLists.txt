
cmake_minimum_required(VERSION 3.12)
set(CMAKE_VERBOSE_MAKEFILE ON)
project(simple_fft)

add_compile_options(-Wall -Wextra -Wpedantic -Werror)
set(CMAKE_C_STANDARD 11)

# Add include directory to the include path
include_directories(include)

# Add subdirectories
add_subdirectory(floating-point)
add_subdirectory(fixed-point)

if(BUILD_TESTSET)
    add_subdirectory(testset)
endif()

target_link_libraries(fixed-fft)
target_link_libraries(float-fft)