cmake_minimum_required(VERSION 3.5.0)

project(linking_example)

find_library(MATRICES Matrix_Utilities)

file(GLOB SOURCES src/*.c)

include_directories(include)

add_executable(${PROJECT_NAME} ${SOURCES})

install(TARGETS ${PROJECT_NAME} DESTINATION lib/${PROJECT_NAME})

target_link_libraries(linking_example PRIVATE ${MATRICES})