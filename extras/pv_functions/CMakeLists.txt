cmake_minimum_required(VERSION 2.8)
project(main)
message(status "${CMAKE_CURRENT_SOURCE_DIR}")
include_directories(src include)
set(HEADER_FILES /Users/jamesmckay/Documents/Programs/tsil-1.3/tsil_cpp.h  include/pv.hpp include/mdm.hpp include/figures.hpp include/decays.hpp include/data.hpp)
set(SRC_FILES  src/main.cpp src/pv.cpp src/mdm.cpp src/figures.cpp src/decays.cpp)


set(Boost_INCLUDE_DIR /usr/local/Cellar/boost/1.59.0/include)
set(Boost_LIBRARY_DIR /usr/local/Cellar/boost/1.59.0/lib)

find_package(Boost COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

LINK_DIRECTORIES(/Users/jamesmckay/Documents/Programs/tsil-1.3)


set(CMAKE_CXX_FLAGS -L.)

add_executable(main ${SRC_FILES} ${HEADER_FILES})

target_link_libraries( main ${Boost_LIBRARIES} )
TARGET_LINK_LIBRARIES(main libtsil.a)

if (CMAKE_COMPILER_IS_GNUCXX)
    message(status "Setting GCC flags")

    # Remove the -g flag when not debugging
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O2 -Wall -Wextra -W -pedantic -std=c++0x -O3")
else()
    message(status "Setting MSVC flags")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHc-")
endif()
install(TARGETS main DESTINATION ${CMAKE_BINARY_DIR}/lib)