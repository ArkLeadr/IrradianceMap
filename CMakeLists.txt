project(IrradianceMap)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules/")

aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu11")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} -std=gnu11 -Ofast")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS} -std=gnu11 -g")


find_package(DevIL REQUIRED) # DevIL


include_directories(. ${DEVIL_INCLUDE_DIR})

target_link_libraries(${PROJECT_NAME} ${DEVIL_LIBRARY})


