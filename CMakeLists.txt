cmake_minimum_required(VERSION 3.17)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

project(sgrid)

find_package(pybind11 REQUIRED)

find_package(Eigen3 REQUIRED)

find_package(VTK REQUIRED)

find_package(pugixml)

set(SOURCE_CODE Sgrid.cpp SgridVtk.cpp SgridGetSet.cpp)

add_library(${PROJECT_NAME} ${SOURCE_CODE})

target_link_libraries(${PROJECT_NAME} PUBLIC
        ${VTK_LIBRARIES}
        pugixml
        Eigen3::Eigen)

pybind11_add_module(${PROJECT_NAME}_bind wrapper.cpp)

target_link_libraries(${PROJECT_NAME}_bind PRIVATE ${PROJECT_NAME})

set_target_properties(${PROJECT_NAME}_bind PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
