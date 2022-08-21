cmake_minimum_required(VERSION 3.23)
project(AtomProc)

set(ATOMPROC_SRC ${CMAKE_CURRENT_SOURCE_DIR}/AtomProc)
set(CMAKE_CXX_STANDARD 17)
add_executable(AtomProc ${ATOMPROC_SRC}/AtomProc.cpp ${ATOMPROC_SRC}/main.cpp)
if (APPLE)
    target_include_directories(AtomProc PUBLIC /usr/local/Cellar/eigen/3.4.0_1/include/eigen3)
endif()

