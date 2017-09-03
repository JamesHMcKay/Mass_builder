cmake_minimum_required(VERSION 2.8)
project(mass_builder)
message(status "${CMAKE_CURRENT_SOURCE_DIR}")


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# Check for Mathematica
include(cmake/FindMathematica.cmake)
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
get_filename_component(TSIL_PATH ${TSIL_PATH} REALPATH)

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

###### the below is relevant for a build using FlexibleSUSY

# Check for Eigen
find_package(Eigen3 3.1.0)
if(EIGEN3_FOUND)
  include_directories("${EIGEN3_INCLUDE_DIR}")
  message("-- Eigen version: ${EIGEN3_VERSION}")
else()
  message("${BoldRed}   Eigen v3.1.0 or greater not found.")
  message("${BoldYellow} Eigen not found -- this is required for FlexibleSUSY functionality")
endif()

message("${BoldYellow} Eigen found in ${EIGEN3_INCLUDE_DIR}")


# find boost (currently hardcoded to my system)

find_package(Boost 1.41)
if(Boost_FOUND)
  include_directories("${Boost_INCLUDE_DIR}")
	#link_directories("${Boost_LIBRARY_DIR}")
else()
  message("${BoldYellow} Boost not found -- required for mass splittings plotting functions")
endif()

# Check for Gnu Scientific Library (GSL)
include(cmake/FindGSL.cmake)
if(GSL_FOUND)
  if (NOT GSL_INCLUDE_DIRS STREQUAL "")
    include_directories("${GSL_INCLUDE_DIRS}")
  endif()
else()
  message("FlexibleSUSY requires the GSL libraries.")
endif()

# if all required packages are found then link the directories

if (GSL_FOUND AND Boost_FOUND AND EIGEN3_FOUND)
	include_directories(flexiblesusy/src/)
	set(FS_INCLUDE include/flexiblesusy.hpp flexiblesusy/src/utils.h)
	# include flexiblesusy
	include(cmake/flexiblesusy.cmake)
	set(MS mass_splittings)
	set(SRC_MASS_SPLITTINGS ${MS}/decays.hpp ${MS}/decays.cpp ${MS}/splittings.cpp ${MS}/mssm_spectrum.cpp ${MS}/mssm_spectrum.hpp ${MS}/ew_triplet_spectrum.cpp ${MS}/ew_triplet_spectrum.hpp ${MS}/mdm_spectrum.cpp ${MS}/mdm_spectrum.hpp ${MS}/figures.cpp ${MS}/figures.hpp)
endif()


include_directories(src include ${Mathematica_WSTP_INCLUDE_DIR})

# specify location of TSIL header file
set(TSIL_HEADER_FILE ${TSIL_PATH}/tsil_cpp.h )

# set source and header files used in all executables
set(HEADER_FILES  include/cmake_variables.hpp include/templates.hpp include/options.hpp include/bases.hpp include/data.hpp include/self_energy.hpp include/utils.hpp include/write_tsil_ini.hpp)

set(SRC_FILES src/templates.cpp src/options.cpp src/bases.cpp src/self_energy.cpp src/generate_code.cpp src/utils.cpp src/write_tsil_ini.cpp )

# set source file for Mass Builder executable
set(SRC_MB src/mass_builder.cpp  src/compute_amp.cpp src/print_vertices.cpp)
set(HEADER_MB include/compute_amp.hpp)

# set example source and header files here
set(SRC_Scalar examples/Scalar.cpp)
set(SRC_MSSM examples/MSSM.cpp)
set(SRC_itMSSM examples/itMSSM.cpp)
set(SRC_VDM examples/VDM.cpp)
set(SRC_QED examples/QED.cpp)
set(SRC_QED2 examples/QED2.cpp)
set(SRC_MDM examples/MDM.cpp)


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

#set(LIBRARIES ${LIBRARIES} ${GSL_LIBRARIES})


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FS_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${flexiblesusy_LDFLAGS}")

# The main Mass Builder executable
add_executable(mass_builder ${SRC_MB} ${SRC_FILES_GEN} ${HEADER_FILES_GEN} ${HEADER_MB} ${TSIL_HEADER_FILE} ${SRC_FILES} ${HEADER_FILES} ${GENERATED_SRC_FILES})

# Add example executables for each of the models
add_executable(Scalar EXCLUDE_FROM_ALL ${SRC_Scalar} ${SRC_FILES_GEN} ${HEADER_FILES_GEN} ${TSIL_HEADER_FILE} ${SRC_FILES} ${HEADER_FILES} ${GENERATED_SRC_FILES})

add_executable(VDM EXCLUDE_FROM_ALL ${SRC_VDM} ${SRC_FILES_GEN} ${HEADER_FILES_GEN} ${TSIL_HEADER_FILE} ${SRC_FILES} ${HEADER_FILES} ${GENERATED_SRC_FILES})

add_executable(QED EXCLUDE_FROM_ALL ${SRC_QED} ${SRC_FILES_GEN} ${HEADER_FILES_GEN} ${TSIL_HEADER_FILE} ${SRC_FILES} ${HEADER_FILES} ${GENERATED_SRC_FILES})

add_executable(QED2 EXCLUDE_FROM_ALL ${SRC_QED2} ${SRC_FILES_GEN} ${HEADER_FILES_GEN} ${TSIL_HEADER_FILE} ${SRC_FILES} ${HEADER_FILES} ${GENERATED_SRC_FILES})

add_executable(MDM EXCLUDE_FROM_ALL ${SRC_MDM} ${SRC_FILES_GEN} ${HEADER_FILES_GEN} ${TSIL_HEADER_FILE} ${SRC_FILES} ${HEADER_FILES} ${GENERATED_SRC_FILES})


# Mass Builder required the WSTP and MPI libraries
TARGET_LINK_LIBRARIES(mass_builder libtsil.a ${MPI_LIBRARIES} ${LIBRARIES})

Mathematica_ABSOLUTIZE_LIBRARY_DEPENDENCIES(mass_builder)


# Examples only require the TSIL libraries

TARGET_LINK_LIBRARIES(MDM libtsil.a)

TARGET_LINK_LIBRARIES(VDM libtsil.a)

TARGET_LINK_LIBRARIES(Scalar libtsil.a)

TARGET_LINK_LIBRARIES(QED libtsil.a)

TARGET_LINK_LIBRARIES(QED2 libtsil.a)


# mass splittings target
LINK_DIRECTORIES(${PROJECT_SOURCE_DIR}/flexiblesusy/src/ ${PROJECT_SOURCE_DIR}/flexiblesusy/models/EW_triplet/ ${PROJECT_SOURCE_DIR}/flexiblesusy/models/MSSM/ ${EIGEN3_INCLUDE_DIR})
add_executable(splittings EXCLUDE_FROM_ALL ${FS_INCLUDE} ${SRC_MASS_SPLITTINGS} ${SRC_FILES_GEN} ${HEADER_FILES_GEN} ${TSIL_HEADER_FILE} ${SRC_FILES} ${HEADER_FILES} ${GENERATED_SRC_FILES})
Mathematica_ABSOLUTIZE_LIBRARY_DEPENDENCIES(splittings)
TARGET_LINK_LIBRARIES(splittings libflexisusy.a libMSSM.a libEW_triplet.a libtsil.a ${GSL_LIBRARIES} ${LIBRARIES} ${Boost_LIBRARIES} )


if (CMAKE_COMPILER_IS_GNUCXX)
    message(status "Setting GCC flags")

    # Remove the -g flag when not debugging
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O2 -Wall -Wextra -W -pedantic -std=c++0x -O3")
else()
    message(status "Setting MSVC flags")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHc-")
endif()
install(TARGETS mass_builder DESTINATION ${CMAKE_BINARY_DIR}/lib)
