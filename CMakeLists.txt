cmake_minimum_required(VERSION 3.9)
project(Custom-SmoothD_for_AviUtl)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

include_directories(./include)
add_library(Custom-SmoothDfA SHARED ./src/SmoothDfA.cpp ./src/idct_int32.c ./src/idct_clip_table.c ./src/SmoothDfA.def)
set_target_properties(Custom-SmoothDfA PROPERTIES PREFIX "")
set_target_properties(Custom-SmoothDfA PROPERTIES OUTPUT_NAME "Custom-SmoothDfA")
set_target_properties(Custom-SmoothDfA PROPERTIES SUFFIX ".auf")

include(CheckIPOSupported)
check_ipo_supported(RESULT supported OUTPUT error)

if( supported )
    message(STATUS "IPO / LTO enabled")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(STATUS "IPO / LTO not supported: <${error}>")
endif()

if(MSVC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W4 /fp:fast")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /Od /DDEBUG")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /O2")

    set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -march=x86-64-v3")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g3 -O0 -DDEBUG")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -Ofast")

    set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")

    set(CMAKE_SHARED_LINKER_FLAGS "-static ${CMAKE_SHARED_LINKER_FLAGS} -s -Wl,--enable-stdcall-fixup,-kill-at")
endif()
