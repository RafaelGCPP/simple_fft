
cmake_minimum_required(VERSION 3.12)
set(CMAKE_VERBOSE_MAKEFILE ON)

project(float-fft)


add_compile_options(-Wall -Wextra -Wpedantic -Werror)
set(CMAKE_C_STANDARD 99)

# Add your source files here
set(SOURCES    
    bitrev.c
    fft.c
    fft_core.c
    fft_init.c
    rfft.c
)

# Add include directory to the include path
#include_directories(../include)


add_library(float-fft STATIC ${SOURCES})
target_include_directories(float-fft PUBLIC ../include)
target_link_libraries(float-fft PRIVATE m)
