cmake_minimum_required(VERSION 2.8)
project(mass_builder)
message(status "${CMAKE_CURRENT_SOURCE_DIR}")

include_directories(src include)
set(TSIL_HEADER_FILE /Users/jamesmckay/Documents/Programs/tsil-1.3/tsil_cpp.h ) 
set(HEADER_FILES include/options.hpp include/bases.hpp include/data.hpp include/self_energy.hpp)
set(SRC_FILES src/options.cpp src/bases.cpp src/self_energy.cpp )

set(SRC_example examples/example_mass_splitting.cpp)
set(SRC_VDM examples/VDM.cpp)

set(SRC_MB src/mass_builder.cpp)

set(HEADER_TRIPLET include/supplements.hpp)
set(SRC_TRIPLET examples/example_triplet.cpp examples/supplements.cpp)


set(HEADER_FILES_GEN include/calc_amplitudes.hpp include/calc_counter_terms.hpp include/utils.hpp include/write_tsil_ini.hpp include/print_vertices.hpp)
set(SRC_FILES_GEN src/calc_amplitudes.cpp src/calc_counter_terms.cpp src/generate_code.cpp src/utils.cpp src/write_tsil_ini.cpp src/print_vertices.cpp)

include(src/CMakeLists.txt)

LINK_DIRECTORIES(/Users/jamesmckay/Documents/Programs/tsil-1.3) 

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_CXX_FLAGS -L.)

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} ${MPI_LINK_FLAGS}")
endif()

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

message("${MPI_LIBRARIES}")

add_executable(mass_builder ${SRC_MB} ${SRC_FILES_GEN} ${HEADER_FILES_GEN} ${TSIL_HEADER_FILE} ${SRC_FILES} ${HEADER_FILES} ${GENERATED_SRC_FILES})

add_executable(example EXCLUDE_FROM_ALL ${SRC_example} ${SRC_FILES_GEN} ${HEADER_FILES_GEN} ${TSIL_HEADER_FILE} ${SRC_FILES} ${HEADER_FILES} ${GENERATED_SRC_FILES})

add_executable(VDM EXCLUDE_FROM_ALL ${SRC_VDM} ${SRC_FILES_GEN} ${HEADER_FILES_GEN} ${TSIL_HEADER_FILE} ${SRC_FILES} ${HEADER_FILES} ${GENERATED_SRC_FILES})

add_executable(triplet EXCLUDE_FROM_ALL ${SRC_TRIPLET} ${HEADER_TRIPLET} ${SRC_FILES_GEN} ${HEADER_FILES_GEN} ${TSIL_HEADER_FILE} ${SRC_FILES} ${HEADER_FILES} ${GENERATED_SRC_FILES})

TARGET_LINK_LIBRARIES(mass_builder libtsil.a ${MPI_LIBRARIES})


TARGET_LINK_LIBRARIES(example libtsil.a)

TARGET_LINK_LIBRARIES(triplet libtsil.a)

TARGET_LINK_LIBRARIES(VDM libtsil.a)

if (CMAKE_COMPILER_IS_GNUCXX)
    message(status "Setting GCC flags")

    # Remove the -g flag when not debugging
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O2 -Wall -Wextra -W -pedantic -std=c++0x -O3")
else()
    message(status "Setting MSVC flags")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHc-")
endif()
install(TARGETS mass_builder DESTINATION ${CMAKE_BINARY_DIR}/lib)
