cmake_minimum_required(VERSION 3.6)
project(sturdins)

string(ASCII 27 Esc)
set(Reset       "${Esc}[m")
set(Bold        "${Esc}[1m")
set(Red         "${Esc}[31m")
set(Green       "${Esc}[32m")
set(Yellow      "${Esc}[33m")
set(Blue        "${Esc}[34m")
set(Magenta     "${Esc}[35m")
set(Cyan        "${Esc}[36m")
set(White       "${Esc}[37m")
set(BoldRed     "${Esc}[1;31m")
set(BoldGreen   "${Esc}[1;32m")
set(BoldYellow  "${Esc}[1;33m")
set(BoldBlue    "${Esc}[1;34m")
set(BoldMagenta "${Esc}[1;35m")
set(BoldCyan    "${Esc}[1;36m")
set(BoldWhite   "${Esc}[1;37m")

# --- Version Control ---
set(STURDR_MAJOR_VERSION 1)
set(STURDR_MINOR_VERSION 0)
set(STURDR_PATCH_VERSION 0)
set(PROJECT_VERSION ${STURDR_MAJOR_VERSION}.${STURDR_MINOR_VERSION}.${STURDR_PATCH_VERSION})
message(STATUS "${BoldCyan}${PROJECT_NAME} version ${PROJECT_VERSION}${Reset}")

# --- Compiler Flags ---
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -W")
# set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Ofast -ffast-math -fomit-frame-pointer -fstrict-aliasing -flto -DNDEBUG -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")

# --- Add Dependencies ---
find_package(Eigen3 REQUIRED)
find_package(navtools REQUIRED)
find_package(satutils REQUIRED)

set(STURDINS_HDRS
    include/sturdins/inertial-nav.hpp
    include/sturdins/kinematic-nav.hpp
    include/sturdins/least-squares.hpp
    include/sturdins/nav-clock.hpp
    include/sturdins/nav-imu.hpp
    include/sturdins/strapdown.hpp
)

set(STURDINS_SRCS
    src/inertial-nav.cpp
    src/kinematic-nav.cpp
    src/least-squares.cpp
    src/nav-clock.cpp
    src/nav-imu.cpp
    src/strapdown.cpp
)

# --- Create the C++ Library ---
add_library(${PROJECT_NAME} ${STURDINS_HDRS} ${STURDINS_SRCS})
target_include_directories(
    ${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${Eigen3_INCLUDE_DIRS}>
    $<INSTALL_INTERFACE:${Eigen3_INCLUDE_DIRS}>
)
target_link_libraries(
    ${PROJECT_NAME} PUBLIC
    Eigen3::Eigen
    navtools
    satutils
)
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

# --- Add Executables ---
if (NOT DEFINED INSTALL_STURDINS_TESTS OR NOT INSTALL_STURDINS_TESTS)
else()
    file(GLOB testfiles "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp")
    foreach(file ${testfiles})
        # message(${file})
        get_filename_component(testname ${file} NAME_WE)
        add_executable(${testname} ${file})
        target_link_libraries( ${testname} PUBLIC ${PROJECT_NAME})
    endforeach(file ${testfiles})
endif()

# --- Make Library 'Findable' for other CMake Packages ---
include(CMakePackageConfigHelpers)

install(
    DIRECTORY include/ 
    DESTINATION include
)

install(
    TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}Targets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    INCLUDES DESTINATION include
)

export(
    TARGETS ${PROJECT_NAME} 
    FILE ${PROJECT_NAME}Targets.cmake
)

install(
    EXPORT ${PROJECT_NAME}Targets 
    FILE ${PROJECT_NAME}Targets.cmake 
    DESTINATION lib/cmake/${PROJECT_NAME}
)

configure_package_config_file(
    "${PROJECT_NAME}-config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
    INSTALL_DESTINATION lib/cmake/${PROJECT_NAME}
    PATH_VARS PROJECT_NAME
)

install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
    DESTINATION lib/cmake/${PROJECT_NAME}
)

# --- Python ---
if (NOT DEFINED INSTALL_PYTHON OR INSTALL_PYTHON)
    find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
    find_package(pybind11 CONFIG REQUIRED)
    python_add_library(_sturdins_core MODULE src/sturdins-python.cpp WITH_SOABI)
    target_link_libraries(_sturdins_core PRIVATE pybind11::headers ${PROJECT_NAME} navtools satutils)
    target_compile_definitions(_sturdins_core PRIVATE VERSION_INFO=${PROJECT_VERSION})
    install(TARGETS _sturdins_core DESTINATION sturdins)
endif()
