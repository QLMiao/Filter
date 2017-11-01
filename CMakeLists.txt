cmake_minimum_required(VERSION 2.8)
project(fft_test)

set(PROJECT ${CMAKE_CURRENT_LIST_DIR})
set(FFTW_PATH ${PROJECT}/fftw)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")

file(GLOB SRC_FILE
	${PROJECT}/*.c
	${PROJECT}/src/*.c
	)


include_directories(${FFTW_PATH}/include)
link_directories(${FFTW_PATH}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT}/bin)
add_executable(project ${SRC_FILE})
target_link_libraries(project fftw3 m)
