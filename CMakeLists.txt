cmake_minimum_required(VERSION 3.20)
project(Custom-SmoothD_for_AviUtl)

set_property(GLOBAL PROPERTY EXPORT_COMPILE_COMMANDS TRUE)
string(TOLOWER ${CMAKE_BUILD_TYPE} build_type)

# cmake tool
include(CheckIPOSupported)

# set source files
add_library(Custom-SmoothDfA SHARED ./src/SmoothDfA.cpp ./src/idct_int32.c ./src/idct_clip_table.c ./src/SmoothDfA.def)
target_include_directories(Custom-SmoothDfA PUBLIC ./include)

# set output file name
set_target_properties(Custom-SmoothDfA PROPERTIES PREFIX "")
set_target_properties(Custom-SmoothDfA PROPERTIES OUTPUT_NAME "Custom-SmoothDfA")
set_target_properties(Custom-SmoothDfA PROPERTIES SUFFIX ".auf")

# add option library
if(MIMALLOC)
    add_subdirectory(library/mimalloc)
    target_link_libraries(Custom-SmoothDfA PUBLIC mimalloc-static)
    target_compile_definitions(Custom-SmoothDfA PRIVATE MIMALLOC)
endif()

# LTO enabled
if(build_type STREQUAL release)
    check_ipo_supported(RESULT supported OUTPUT error)
    if(supported)
        message(STATUS "IPO / LTO enabled")
        set_property(GLOBAL PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    else()
        message(STATUS "IPO / LTO not supported: <${error}>")
    endif()
endif()

# set build options
target_compile_features(Custom-SmoothDfA PUBLIC PRIVATE cxx_std_11)

if(build_type STREQUAL debug)
    target_compile_definitions(Custom-SmoothDfA PRIVATE DEBUG)
endif()

if(MSVC)
    target_compile_options(Custom-SmoothDfA PUBLIC /arch:AVX2)
    if(build_type STREQUAL debug)
        target_compile_options(Custom-SmoothDfA PUBLIC /Od)
    else()
        target_compile_options(Custom-SmoothDfA PUBLIC /O2 /fp:fast)
    endif()
else()
    target_compile_options(Custom-SmoothDfA PUBLIC -march=x86-64-v3)
    if(NOT WIN32)
        target_link_options(Custom-SmoothDfA PUBLIC -static -s -Wl,--enable-stdcall-fixup -Wl,-kill-at)
    endif()
    if(build_type STREQUAL debug)
        target_compile_options(Custom-SmoothDfA PUBLIC -g3)
    else()
        target_compile_options(Custom-SmoothDfA PUBLIC -Ofast)
    endif()
endif()
