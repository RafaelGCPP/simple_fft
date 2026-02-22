
cmake_minimum_required(VERSION 3.12)
set(CMAKE_VERBOSE_MAKEFILE ON)

project(fixed-fft)


add_compile_options(-Wall -Wextra -Wpedantic -Werror)
set(CMAKE_C_STANDARD 99)

# Add your source files here
set(SOURCES  
    fix_bitrev.c
    fix_fft.c
    fix_rfft.c
    fix_fft_core.c
    fix_fft_init.c
    cordic.c
)

# Add include directory to the include path

add_library(fixed-fft STATIC ${SOURCES})
target_include_directories(fixed-fft PUBLIC ../include)