cmake_minimum_required(VERSION 3.5.0)

project(Matrix_Utilities)

file(GLOB SOURCES src/*.c)

include_directories(include)

add_library(${PROJECT_NAME} SHARED ${SOURCES})

install(TARGETS ${PROJECT_NAME} DESTINATION lib)

file(GLOB HEADERS include/*.h)
install(FILES ${HEADERS} DESTINATION include/${PROJECT_NAME})
