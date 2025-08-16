# Simple app so just simply make a variable of target names
SET(TARGET "transcendental-integrals")
SET(TARGET_SRC ${TARGET}.cc)

# Minimum CMake and requiring boost
CMAKE_MINIMUM_REQUIRED(VERSION 4.0.1)
set(Boost_DEBUG 1)
#find_package(Boost REQUIRED)

# Including stuff
include_directories(${Boost_INCLUDE_DIRS})
add_executable(${TARGET} ${TARGET_SRC})

# Project templating
project(${TARGET})