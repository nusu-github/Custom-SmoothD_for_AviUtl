cmake_minimum_required(VERSION 3.23)
project(Custom-SmoothD_for_AviUtl)

set(CMAKE_CXX_STANDARD 20)

# cmake tool
include(CheckIPOSupported)

# in-source build error
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "In-source builds are not allowed.
    Please make a new directory (called a build directory) and run CMake from there.
    You may need to remove CMakeCache.txt.")
endif ()

# set source files
add_library(Custom-SmoothDfA SHARED ./src/SmoothDfA.cpp ./src/ht.cpp ./src/dct_int32.c ./src/idct_clip_table.c ./src/idct_int32.c ./src/SmoothDfA.def)
target_include_directories(Custom-SmoothDfA PUBLIC "./include")
target_include_directories(Custom-SmoothDfA PUBLIC "./include/aviutl_exedit_sdk")

# set output file name
set_target_properties(Custom-SmoothDfA PROPERTIES PREFIX "")
set_target_properties(Custom-SmoothDfA PROPERTIES OUTPUT_NAME "Custom-SmoothDfA")
set_target_properties(Custom-SmoothDfA PROPERTIES SUFFIX ".auf")

# openmp enable
find_package(OpenMP REQUIRED)
if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(Custom-SmoothDfA PRIVATE /openmp:experimental)
else ()
    target_link_libraries(Custom-SmoothDfA PUBLIC OpenMP::OpenMP_CXX)
endif ()

# LTO enabled
if (ENABLE_LTO)
    check_ipo_supported(RESULT result OUTPUT output)
    if (result)
        message(STATUS "IPO / LTO enabled")
        set_property(GLOBAL PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
        if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
            target_compile_options(Custom-SmoothDfA PRIVATE /GL /Gw)
            target_link_options(Custom-SmoothDfA PRIVATE /LTCG /OPT:REF,ICF)
        endif ()
    else ()
        message(WARNING "LTO is not supported: ${output}")
    endif ()
endif ()

# PGO enabled
if (PGO_USE OR PGO_2PASS_USE)
    if (PGO_2PASS_USE)
        file(COPY 2pass_code.profdata DESTINATION ${CMAKE_BINARY_DIR})
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-use=2pass_code.profdata")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-use=2pass_code.profdata")
    else ()
        file(COPY code.profdata DESTINATION ${CMAKE_BINARY_DIR})
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-use=code.profdata")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-use=code.profdata")
    endif ()
    message(STATUS "PGO use enabled")
endif ()

if (PGO_GENERATE OR PGO_2PASS_GENERATE)
    if (PGO_2PASS_GENERATE)
        if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fcs-profile-generate=2pass_code-%m.profraw")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fcs-profile-generate=2pass_code-%m.profraw")
        else ()
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-generate=2pass_code-%m.profraw")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-generate=2pass_code-%m.profraw")
        endif ()
    else ()
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-generate=code-%m.profraw")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-generate=code-%m.profraw")
    endif ()
    message(STATUS "PGO generation enabled")
endif ()

if (CMAKE_BUILD_TYPE STREQUAL debug)
    target_compile_definitions(Custom-SmoothDfA PRIVATE DEBUG)
endif ()

# set build options

if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    add_compile_options(/arch:AVX2)
    if (CMAKE_BUILD_TYPE STREQUAL debug)
        add_compile_options(/Od)
    else ()
        add_compile_options(/O2)
        if (MAX_OPTIMIZE)
            add_compile_options(/fp:fast /Ob3)
        endif ()
    endif ()
endif ()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GUN")
        add_compile_options(-mfpmath=sse)
    endif ()

    add_compile_options(-mvzeroupper)
    target_link_options(Custom-SmoothDfA PUBLIC -Wl,--enable-stdcall-fixup -Wl,-kill-at)

    if (NATIVE)
        add_compile_options(-march=native -mtune=native -mcpu=native)
    else ()
        add_compile_options(-march=x86-64-v3 -mtune=x86-64-v3 -mcpu=x86-64-v3)
    endif ()

    if (CMAKE_BUILD_TYPE STREQUAL debug)
        add_compile_options(-O0)
    else ()
        if (MAX_OPTIMIZE)
            add_compile_options(-Ofast)
        else ()
            add_compile_options(-O3 -ffast-math)
        endif ()
    endif ()
endif ()