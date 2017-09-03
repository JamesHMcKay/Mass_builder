#Copyright (c) 2017, The GAMBIT Collaboration
#All rights reserved.

#Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

#1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#interfaces with external packages (via a "backend" system), a complete
#2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

#3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


include(ExternalProject)

set(FS_DIR "${PROJECT_SOURCE_DIR}/flexiblesusy/")
set(FS_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2 -fPIC")
set(FS_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")

# Determine compiler libraries needed by flexiblesusy.
if(CMAKE_Fortran_COMPILER MATCHES "gfortran*")
	set(flexiblesusy_extralibs "${flexiblesusy_extralibs} -lgfortran -lm")
elseif(CMAKE_Fortran_COMPILER MATCHES "g77" OR CMAKE_Fortran_COMPILER MATCHES "f77")
	set(flexiblesusy_extralibs "${flexiblesusy_extralibs} -lg2c -lm")
elseif(CMAKE_Fortran_COMPILER MATCHES "ifort")
	set(flexiblesusy_extralibs "${flexiblesusy_extralibs} -lifcore -limf -ldl -lintlc -lsvml")
endif()
message("${BoldYellow}-- Determined FlexibleSUSY compiler library dependencies: ${flexiblesusy_extralibs}${ColourReset}")
set(flexiblesusy_LDFLAGS "${flexiblesusy_LDFLAGS} ${flexiblesusy_extralibs}")

# Silence the deprecated-declarations warnings comming from Eigen3
#set_compiler_warning("no-deprecated-declarations" FS_CXX_FLAGS)

# Silence the unused parameter and variable warnings comming from FlexibleSUSY
#set_compiler_warning("no-unused-parameter" FS_CXX_FLAGS)
#set_compiler_warning("no-unused-variable" FS_CXX_FLAGS)

# FlexibleSUSY configure options
set(FS_OPTIONS ${FS_OPTIONS}
		 --with-cxx=${CMAKE_CXX_COMPILER}
		 --with-cxxflags=${FS_CXX_FLAGS}
		 --with-shared-ldflags=${OpenMP_CXX_FLAGS}
		 --with-fc=${CMAKE_Fortran_COMPILER}
		 --with-fflags=${FS_Fortran_FLAGS}
		 --with-eigen-incdir=${EIGEN3_INCLUDE_DIR}
		 --with-boost-libdir=${Boost_LIBRARY_DIR}
		 --with-boost-incdir=${Boost_INCLUDE_DIR}
		 --with-lapack-libs=${LAPACK_LINKLIBS}
		 --with-blas-libs=${LAPACK_LINKLIBS}
		 --disable-librarylink
		 --with-fc=gfortran
		 --disable-threads
		#--enable-verbose flag causes verbose output at runtime as well. Maybe set it dynamically somehow in future.
	 )

# Set the models (spectrum generators) existing in flexiblesusy (could autogen this, but that would build some things we don't need)
set(BUILT_FS_MODELS EW_triplet MSSM MDM)

# Explain how to build each of the flexiblesusy spectrum generators we need.  Configure now, serially, to prevent parallel build issues.
string (REPLACE ";" "," BUILT_FS_MODELS_COMMAS "${BUILT_FS_MODELS}")
set(config_command ./configure ${FS_OPTIONS} --with-models=${BUILT_FS_MODELS_COMMAS})
add_custom_target(configure-flexiblesusy COMMAND cd ${FS_DIR} && ${config_command})
message("${Yellow}-- Configuring FlexibleSUSY for models: ${BoldYellow}${BUILT_FS_MODELS_COMMAS}${ColourReset}")
execute_process(COMMAND ${config_command}
								WORKING_DIRECTORY ${FS_DIR}
								RESULT_VARIABLE result
								OUTPUT_VARIABLE output
							 )
if (NOT "${result}" STREQUAL "0")
	message("${BoldRed}-- Configuring FlexibleSUSY failed.  Here's what I tried to do:\n${config_command}\n${output}${ColourReset}" )
	message(FATAL_ERROR "Configuring FlexibleSUSY failed." )
endif()
set(rmstring "${CMAKE_BINARY_DIR}/flexiblesusy-prefix/src/flexiblesusy-stamp/flexiblesusy")
execute_process(COMMAND ${CMAKE_COMMAND} -E touch ${rmstring}-configure)

message("${Yellow}-- Configuring FlexibleSUSY - done.${ColourReset}")

# Add FlexibleSUSY as an external project
ExternalProject_Add(flexiblesusy
	SOURCE_DIR ${FS_DIR}
	BUILD_IN_SOURCE 1
	BUILD_COMMAND $(MAKE) alllib
	CONFIGURE_COMMAND ${config_command}
	INSTALL_COMMAND ""
)

# Set linking commands.  Link order matters! The core flexiblesusy libraries need to come after the model libraries but before the other link flags.
set(flexiblesusy_LDFLAGS "-L${FS_DIR}/src -lflexisusy -L${FS_DIR}/legacy -llegacy ${flexiblesusy_LDFLAGS}")
foreach(_MODEL ${BUILT_FS_MODELS})
	set(flexiblesusy_LDFLAGS "-L${FS_DIR}/models/${_MODEL} -l${_MODEL} ${flexiblesusy_LDFLAGS}")
endforeach()

# Set up include paths
include_directories("${FS_DIR}/..")
include_directories("${FS_DIR}/src")
include_directories("${FS_DIR}/legacy")
include_directories("${FS_DIR}/config")
include_directories("${FS_DIR}/slhaea")
# Dig through flexiblesusy "models" directory and add all subdirectories to the include list
# (these contain the headers for the generated spectrum generators)
foreach(_MODEL ${BUILT_FS_MODELS})
	include_directories("${FS_DIR}/models/${_MODEL}")
endforeach()

# Strip out leading and trailing whitespace
string(STRIP "${flexiblesusy_LDFLAGS}" flexiblesusy_LDFLAGS)
