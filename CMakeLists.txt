cmake_minimum_required(VERSION 3.10)
project(executable VERSION 1.0 LANGUAGES CXX)

set(EIGEN_INCLUDE_DIR "/usr/local/include/eigen3")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include_directories(${EIGEN_INCLUDE_DIR})

# Установка директории для исполняемых файлов
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(${PROJECT_NAME}
    src/main.cpp
    src/init/init.cpp
    src/curve/curve.cpp
    src/preprocess/preprocess.cpp
    src/system/system.cpp
)

target_include_directories(${PROJECT_NAME}
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)

install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)
