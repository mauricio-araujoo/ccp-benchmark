cmake_minimum_required(VERSION 3.20)

project(ccp-benchmark)

set(CMAKE_CXX_STANDARD 20)
set(FETCHCONTENT_QUIET OFF CACHE BOOL "" FORCE)

# Add CPM
set(CPM_DOWNLOAD_VERSION 0.40.2)
set(CPM_DOWNLOAD_LOCATION "${CMAKE_BINARY_DIR}/cmake/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
if(NOT (EXISTS ${CPM_DOWNLOAD_LOCATION}))
    message(STATUS "Downloading CPM.cmake...")
    file(DOWNLOAD https://github.com/cpm-cmake/CPM.cmake/releases/download/v${CPM_DOWNLOAD_VERSION}/CPM.cmake ${CPM_DOWNLOAD_LOCATION})
endif()
include(${CPM_DOWNLOAD_LOCATION})

# Add external libs
CPMAddPackage("gh:d99kris/rapidcsv@8.83")
CPMAddPackage("gh:p-ranav/argparse@3.1")
add_subdirectory(spark)

# Source files
add_executable(ccp-benchmark
        src/main.cpp
        src/simulation.cpp
        src/parameters.cpp
        src/parameters.h
        src/simulation.h
        src/reactions.cpp
        src/simulation_events.cpp
        src/simulation_events.h
)

target_link_libraries(ccp-benchmark PUBLIC spark::spark rapidcsv argparse)