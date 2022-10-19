cmake_minimum_required(VERSION 3.23.2)

project(
        Custom-SmoothD_for_AviUtl
        VERSION 0.01.1
        LANGUAGES CXX C
)

# in-source build error
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "In-source builds are not allowed.
    Please make a new directory (called a build directory) and run CMake from there.
    You may need to remove CMakeCache.txt.")
endif ()

# Option list
## compile option
option(STATIC-LINK "Enable static compile" ON)
option(ENABLE_LTO "Enable LTO" OFF)
option(ENABLE_PGO_GEN "Enable PGO generate" OFF)
option(ENABLE_PGO_USE "Enable PGO use" OFF)
option(NATIVE "Enable native compile" OFF)
option(ENABLE_AVX512 "Enable AVX512" OFF)
option(ENABLE_AVX2 "Enable AVX2" ON)

# cmake tool
include(CheckIPOSupported)
include(FetchContent)

# Library list
## aviutl_exedit_sdk
FetchContent_Declare(
        aviutl_exedit_sdk
        GIT_REPOSITORY "https://github.com/ePi5131/aviutl_exedit_sdk"
        GIT_SHALLOW true
        GIT_PROGRESS true
)
## eve
FetchContent_Declare(
        eve
        GIT_REPOSITORY "https://github.com/jfalcou/eve"
        GIT_TAG "main"
        GIT_SHALLOW true
        GIT_PROGRESS true
)
FetchContent_MakeAvailable(aviutl_exedit_sdk)

# Project settings
## Set source files
add_library(Custom-SmoothDfA SHARED ./src/SmoothDfA.cpp ./src/ht.cpp src/dct_int32.cpp src/idct_int32.cpp ./src/SmoothDfA.def)

## Set include directories
target_include_directories(Custom-SmoothDfA PRIVATE ./include)
### aviutl_exedit_sdk
target_include_directories(Custom-SmoothDfA PRIVATE ${aviutl_exedit_sdk_SOURCE_DIR})
### eve
set(EVE_BUILD_TEST OFF CACHE BOOL false FORCE)
set(EVE_BUILD_BENCHMARKS OFF CACHE BOOL false FORCE)
set(EVE_BUILD_DOCUMENTATION OFF CACHE BOOL false FORCE)
FetchContent_MakeAvailable(eve)
target_link_libraries(Custom-SmoothDfA PRIVATE eve::eve)

## C++ standard
target_compile_features(Custom-SmoothDfA PRIVATE cxx_std_20)

## set output file name
set_target_properties(Custom-SmoothDfA PROPERTIES PREFIX "")
### Release or Debug
if (CMAKE_BUILD_TYPE STREQUAL "Release")
    set_target_properties(Custom-SmoothDfA PROPERTIES OUTPUT_NAME "Custom-SmoothDfA")
else ()
    set_target_properties(Custom-SmoothDfA PROPERTIES OUTPUT_NAME "Custom-SmoothDfA_D")
endif ()
set_target_properties(Custom-SmoothDfA PROPERTIES SUFFIX ".auf")

### OpenMP
find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found")
    message(STATUS "OpenMP_CXX_VERSION: ${OpenMP_CXX_VERSION}")
    if (OpenMP_CXX_VERSION VERSION_LESS 4.5)
        if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
            target_compile_options(Custom-SmoothDfA PRIVATE -openmp:experimental)
            message(WARNING
                    "OpenMP4.0 in the MSVC compiler is an experimental feature. It may fail to compile.
                    Please use the Clang compiler.")
        else ()
            message(FATAL_ERROR "OpenMP version is too old. Please update your compiler.")
        endif ()
    else ()
        if (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
            target_compile_options(Custom-SmoothDfA PRIVATE -Qopenmp-simd)
        else ()
            target_compile_options(Custom-SmoothDfA PRIVATE -fopenmp-simd)
        endif ()
    endif ()
else ()
    message(FATAL_ERROR "OpenMP not found. Please install OpenMP Library.")
endif ()

### LTO
check_ipo_supported(RESULT LTO_SUPPORTED OUTPUT LTO_OUTPUT)
if (LTO_SUPPORTED)
    message(STATUS "LTO is supported")
    if (ENABLE_LTO)
        message(STATUS "LTO is enabled")
        if (CMAKE_BUILD_TYPE STREQUAL "Release")
            set_target_properties(Custom-SmoothDfA PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
            if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
                target_compile_options(Custom-SmoothDfA PRIVATE /GL)
                target_link_options(Custom-SmoothDfA PRIVATE /LTCG)
            endif ()
        else ()
            message(WARNING "LTO is enabled but CMAKE_BUILD_TYPE is not Release")
            set(ENABLE_LTO OFF)
        endif ()
    endif ()
else ()
    if (ENABLE_LTO)
        message(WARNING "LTO is enabled but not supported")
        set(ENABLE_LTO OFF)
    endif ()
endif ()

if (MSVC)
    # Compile detection
    if (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
        message(STATUS "Compiler is MSVC-Like IntelLLVM")
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        message(STATUS "Compiler is MSVC-Like Clang")
        message(FATAL_ERROR "Do not use MSVC-like Clang(clang-cl).")
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        message(STATUS "Compiler is MSVC-Like MSVC")
    else ()
        message(WARNING "Compiler is MSVC-Like but not IntelLLVM, Clang, MSVC")
    endif ()

    target_link_options(Custom-SmoothDfA PRIVATE /LARGEADDRESSAWARE)

    # Release
    if (CMAKE_BUILD_TYPE STREQUAL "Release")

        target_compile_options(Custom-SmoothDfA PRIVATE /GF /Gw /Gy /J)
        target_link_options(Custom-SmoothDfA PRIVATE /OPT:REF,ICF,LBR /INCREMENTAL:NO)

        if (OPTIMIZE_FOR_SIZE)
            message(STATUS "Optimize for size")
            target_compile_options(Custom-SmoothDfA PRIVATE /Os)
        else ()
            message(STATUS "Optimize for speed")
            target_compile_options(Custom-SmoothDfA PRIVATE /Ot)
        endif ()
    else ()
    endif ()

    target_compile_options(Custom-SmoothDfA PRIVATE /fp:fast)

    if (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")

        if (SAVE_TEMPS)
            target_compile_options(Custom-SmoothDfA PRIVATE /Qsave-temps)
        endif ()

        if (NATIVE)
            target_compile_options(Custom-SmoothDfA PRIVATE /QxHost)
        else ()
            if (ENABLE_AVX512)
                target_compile_options(Custom-SmoothDfA PRIVATE /QxAVX512)
            elseif (ENABLE_AVX2)
                target_compile_options(Custom-SmoothDfA PRIVATE /QxAVX2)
            endif ()
        endif ()

        target_compile_options(Custom-SmoothDfA PRIVATE /Qopt-report:2)
        target_compile_options(Custom-SmoothDfA PRIVATE /Qvec /Qopt-multiple-gather-scatter-by-shuffles /Qopt-dynamic-align)

    else ()
        if (ENABLE_AVX512)
            target_compile_options(Custom-SmoothDfA PRIVATE /arch:AVX512)
        elseif (ENABLE_AVX2)
            target_compile_options(Custom-SmoothDfA PRIVATE -arch:AVX2)
        endif ()

        target_compile_options(Custom-SmoothDfA PRIVATE -Qvec-report:2)

    endif ()

else ()
    # Compile detection
    if (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
        message(STATUS "Compiler is IntelLLVM")
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        message(STATUS "Compiler is Clang")
    else ()
        message(STATUS "Compiler is GCC")
    endif ()

    target_link_options(Custom-SmoothDfA PRIVATE -s)

    if (NATIVE)
        target_compile_options(Custom-SmoothDfA PRIVATE -march=native -mtune=native)
    else ()
        if (ENABLE_AVX512)
            target_compile_options(Custom-SmoothDfA PRIVATE -march=x86-64-v4)
        elseif (ENABLE_AVX2)
            target_compile_options(Custom-SmoothDfA PRIVATE -march=x86-64-v3)
        endif ()
    endif ()

    if (ENABLE_PGO_GEN)
        target_compile_options(Custom-SmoothDfA PRIVATE -fprofile-generate -fprofile-update=atomic)
        target_link_options(Custom-SmoothDfA PRIVATE -fprofile-generate -fprofile-update=atomic)
    elseif (ENABLE_PGO_USE)
        target_compile_options(Custom-SmoothDfA PRIVATE -fprofile-use)
    endif ()

    if (CMAKE_BUILD_TYPE STREQUAL "Release")
        target_compile_options(Custom-SmoothDfA PRIVATE -ffast-math -freciprocal-math)
    endif ()

    if (STATIC_LINK)
        target_link_options(Custom-SmoothDfA PRIVATE -static)
    endif ()

    if (SAVE_TEMPS)
        target_compile_options(Custom-SmoothDfA PRIVATE -save-temps)
    endif ()

endif ()
