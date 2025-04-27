cmake_minimum_required(VERSION 3.24)
project(main LANGUAGES C)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

file(GLOB SRC_LIST CONFIGURE_DEPENDS src/*.c)

add_executable(main ${SRC_LIST})

target_include_directories(main PRIVATE inc)

set_target_properties(main PROPERTIES
    C_STANDARD 99
    C_STANDARD_REQUIRED ON
    C_EXTENSIONS OFF
)