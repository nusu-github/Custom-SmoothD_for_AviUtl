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

# LTO enabled
if (ENABLE_LTO)
    check_ipo_supported(RESULT result OUTPUT output)
    if (result)
        message(STATUS "IPO / LTO enabled")
        set_property(GLOBAL PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
        if (MSVC)
            target_compile_options(Custom-SmoothDfA PRIVATE /GL /Gw)
            target_link_options(Custom-SmoothDfA PRIVATE /OPT:REF,ICF)
        endif ()
    else ()
        if (MSVC AND CMAKE_CXX_COMPILER_ID MATCHES "Intel")
            message(STATUS "IPO / LTO enabled")
            target_link_options(Custom-SmoothDfA PRIVATE /Qipo)
        else ()
            message(WARNING "LTO is not supported: ${output}")
        endif ()
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


# set build options
if (MSVC)
    target_compile_options(Custom-SmoothDfA PUBLIC /arch:AVX2)

    if (CMAKE_BUILD_TYPE MATCHES "Debug")
        taget_compile_options(Custom-SmoothDfA PUBLIC /Od /ZI /DEBUG)
    else ()
        target_compile_options(Custom-SmoothDfA PUBLIC /O2)
    endif ()

    if (MAX_OPTIMIZE)
        target_compile_options(Custom-SmoothDfA PUBLIC /fp:fast /Ob3)
    endif ()

    if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        message(STATUS "MSVC-like Intel C++ Compiler detected")

        target_compile_options(Custom-SmoothDfA PUBLIC /Qvec /Qopt-report:3)
        target_link_options(Custom-SmoothDfA PUBLIC /Quse-intel-optimized-headers /Qno-intel-lib)

        if (NATIVE)
            target_compile_options(Custom-SmoothDfA PUBLIC /QxHost)
        endif ()

    endif ()

    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        message(STATUS "MSVC-like Clang detected")

        if (NATIVE)
            target_compile_options(Custom-SmoothDfA PUBLIC /clang:-march=native /clang:-mtune=native)
        else ()
            target_compile_options(Custom-SmoothDfA PUBLIC /clang:-march=x86-64-v3 /clang:-mtune=x86-64)
        endif ()

        if (MAX_OPTIMIZE)
            target_compile_options(Custom-SmoothDfA PUBLIC /clang:-ffast-math)
        endif ()

    endif ()
else ()
    if (STATIC)
        target_link_options(Custom-SmoothDfA PUBLIC -static)
    endif ()

    if (NATIVE)
        target_compile_options(Custom-SmoothDfA PUBLIC -march=native -mtune=native)
    else ()
        target_compile_options(Custom-SmoothDfA PUBLIC -march=x86-64-v3 -mtune=generic)
    endif ()

    if (CMAKE_BUILD_TYPE MATCHES "Debug")
        message(STATUS "Debug build")
        target_compile_options(Custom-SmoothDfA PUBLIC -O0 -g3 -ggdb)
    else ()
        message(STATUS "Release build")
        if (MAX_OPTIMIZE)
            target_link_options(Custom-SmoothDfA PUBLIC -Ofast)
        else ()
            target_link_options(Custom-SmoothDfA PUBLIC -O3 -ffast-math)
        endif ()
    endif ()
endif ()