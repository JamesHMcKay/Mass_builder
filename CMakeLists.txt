cmake_minimum_required(VERSION 2.8)
project(mass_builder)
message(status "${CMAKE_CURRENT_SOURCE_DIR}")


# Check for Mathematica
include(FindMathematica.cmake)
if(Mathematica_FOUND)
  message("${BoldYellow}   Found Mathematica")
  if(Mathematica_WSTP_FOUND)
    message("${BoldYellow}   Found Wolfram Symbolic Transfer Protocol. Mathematica backends enabled")
    set(HAVE_MATHEMATICA 1)
    set(MATHEMATICA_WSTP_H "${Mathematica_WSTP_INCLUDE_DIR}/wstp.h")
    set(MATHEMATICA_KERNEL "${Mathematica_KERNEL_EXECUTABLE}")
  else()
    message("${BoldRed}  WSTP not found. Please make sure it is installed before attempting to use Mathematica backends")
    set(HAVE_MATHEMATICA 0)
  endif()
else()
  message("${BoldRed}   Mathematica not found. Backends using Mathematica will be disabled")
  set(HAVE_MATHEMATICA 0)
endif()


if(Mathematica_FOUND AND Mathematica_WSTP_FOUND)
  set(LIBRARIES ${LIBRARIES} ${Mathematica_WSTP_LIBRARIES} ${Mathematica_MathLink_LIBRARY})
endif()

set(TSIL_PATH "" CACHE STRING "The path to the TSIL file")

# Configure cmake variables file
set(outhdr "${PROJECT_SOURCE_DIR}/include/cmake_variables.hpp")
configure_file("${PROJECT_SOURCE_DIR}/include/cmake_variables.hpp.in" ${outhdr})

set(outhdr "${PROJECT_SOURCE_DIR}/include/self_energy.hpp")
if(NOT EXISTS ${outhdr})
  # We set up skeleton files that will be replaced by generated files by Mass Builder later
  # Configure skeleton self_energy.hpp
  set(outhdr "${PROJECT_SOURCE_DIR}/include/self_energy.hpp")
  configure_file("${PROJECT_SOURCE_DIR}/include/self_energy.hpp.in" ${outhdr})

  # Configure skeleton data.hpp
  set(outhdr "${PROJECT_SOURCE_DIR}/include/data.hpp")
  configure_file("${PROJECT_SOURCE_DIR}/include/data.hpp.in" ${outhdr})

  # Configure skeleton self_energy.cpp
  set(outhdr "${PROJECT_SOURCE_DIR}/src/self_energy.cpp")
  configure_file("${PROJECT_SOURCE_DIR}/src/self_energy.cpp.in" ${outhdr})
endif()

include_directories(src include ${Mathematica_WSTP_INCLUDE_DIR})

# specify location of TSIL header file
set(TSIL_HEADER_FILE ${TSIL_PATH}/tsil_cpp.h )

# set source and header files used in all executables
set(HEADER_FILES  include/cmake_variables.hpp include/templates.hpp include/options.hpp include/bases.hpp include/data.hpp include/self_energy.hpp include/utils.hpp include/write_tsil_ini.hpp)

set(SRC_FILES src/templates.cpp src/options.cpp src/bases.cpp src/self_energy.cpp src/generate_code.cpp src/utils.cpp src/write_tsil_ini.cpp src/print_vertices.cpp)

# set source file for Mass Builder executable
set(SRC_MB src/mass_builder.cpp  src/compute_amp.cpp)
set(HEADER_MB include/compute_amp.hpp)


# set example source and header files here
set(SRC_Scalar examples/Scalar.cpp)
set(SRC_MSSM examples/MSSM.cpp)
set(SRC_VDM examples/VDM.cpp)
set(SRC_TRIPLET examples/EW_triplet.cpp)
set(HEADER_TRIPLET examples/EW_triplet.hpp)

# include generated source files
include(src/CMakeLists.txt)


LINK_DIRECTORIES(${TSIL_PATH})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_CXX_FLAGS -L.)

FIND_PACKAGE(MPI REQUIRED)
if(MPI_FOUND)
message("MPI FOUND")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_LINK_FLAGS}")
include_directories(${MPI_INCLUDE_PATH})
endif()

# The main Mass Builder executable
add_executable(mass_builder ${SRC_MB} ${SRC_FILES_GEN} ${HEADER_FILES_GEN} ${HEADER_MB} ${TSIL_HEADER_FILE} ${SRC_FILES} ${HEADER_FILES} ${GENERATED_SRC_FILES})

# Add example executables for each of the models
add_executable(Scalar EXCLUDE_FROM_ALL ${SRC_Scalar} ${SRC_FILES_GEN} ${HEADER_FILES_GEN} ${TSIL_HEADER_FILE} ${SRC_FILES} ${HEADER_FILES} ${GENERATED_SRC_FILES})

add_executable(MSSM EXCLUDE_FROM_ALL ${SRC_MSSM} ${SRC_FILES_GEN} ${HEADER_FILES_GEN} ${TSIL_HEADER_FILE} ${SRC_FILES} ${HEADER_FILES} ${GENERATED_SRC_FILES})

add_executable(VDM EXCLUDE_FROM_ALL ${SRC_VDM} ${SRC_FILES_GEN} ${HEADER_FILES_GEN} ${TSIL_HEADER_FILE} ${SRC_FILES} ${HEADER_FILES} ${GENERATED_SRC_FILES})

add_executable(EW_triplet EXCLUDE_FROM_ALL ${SRC_TRIPLET} ${HEADER_TRIPLET} ${SRC_FILES_GEN} ${HEADER_FILES_GEN} ${TSIL_HEADER_FILE} ${SRC_FILES} ${HEADER_FILES} ${GENERATED_SRC_FILES})

# Mass Builder required the WSTP and MPI libraries
TARGET_LINK_LIBRARIES(mass_builder libtsil.a ${MPI_LIBRARIES} ${LIBRARIES})

Mathematica_ABSOLUTIZE_LIBRARY_DEPENDENCIES(mass_builder)

# Examples only require the TSIL libraries
TARGET_LINK_LIBRARIES(MSSM libtsil.a)

TARGET_LINK_LIBRARIES(EW_triplet libtsil.a)

TARGET_LINK_LIBRARIES(VDM libtsil.a)

TARGET_LINK_LIBRARIES(Scalar libtsil.a)

if (CMAKE_COMPILER_IS_GNUCXX)
    message(status "Setting GCC flags")

    # Remove the -g flag when not debugging
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O2 -Wall -Wextra -W -pedantic -std=c++0x -O3")
else()
    message(status "Setting MSVC flags")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHc-")
endif()
install(TARGETS mass_builder DESTINATION ${CMAKE_BINARY_DIR}/lib)
