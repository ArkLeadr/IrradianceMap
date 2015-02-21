project(IrradianceMap)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu11")

set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} -std=gnu11 -Ofast")

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS} -std=gnu11 -g")

if(APPLE)
	set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules/")
	find_package(FreeImage REQUIRED)
	include_directories( . ${FREE_IMAGE_INCLUDE_DIR})
	target_link_libraries(${PROJECT_NAME} ${FREE_IMAGE_LIBRARY})
else(APPLE)
	target_link_libraries(${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/FreeImage.lib)
endif(APPLE)