cmake_minimum_required(VERSION 2.8)
project(mass_builder)
message(status "${CMAKE_CURRENT_SOURCE_DIR}")

include_directories(src include)
set(HEADER_FILES /Users/jamesmckay/Documents/Programs/tsil-1.3/tsil_cpp.h  include/data.hpp include/self_energy.hpp include/supplements.hpp)
set(SRC_FILES src/self_energy.cpp src/supplements.cpp )


set(HEADER_FILES_GEN include/calc_amplitudes.hpp include/utils.hpp)
set(SRC_FILES_GEN src/mass_builder.cpp src/calc_amplitudes.cpp src/generate_code.cpp src/utils.cpp)



set(Boost_INCLUDE_DIR /usr/local/Cellar/boost/1.59.0/include)
set(Boost_LIBRARY_DIR /usr/local/Cellar/boost/1.59.0/lib)

find_package(Boost COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

LINK_DIRECTORIES(/Users/jamesmckay/Documents/Programs/tsil-1.3)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_CXX_FLAGS -L.)

add_executable(mass_builder ${SRC_FILES_GEN} ${HEADER_FILES_GEN} ${SRC_FILES} ${HEADER_FILES})

target_link_libraries(mass_builder ${Boost_LIBRARIES} )
TARGET_LINK_LIBRARIES(mass_builder libtsil.a)

if (CMAKE_COMPILER_IS_GNUCXX)
    message(status "Setting GCC flags")

    # Remove the -g flag when not debugging
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O2 -Wall -Wextra -W -pedantic -std=c++0x -O3")
else()
    message(status "Setting MSVC flags")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHc-")
endif()
install(TARGETS mass_builder DESTINATION ${CMAKE_BINARY_DIR}/lib)