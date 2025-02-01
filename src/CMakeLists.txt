cmake_minimum_required(VERSION 3.2)

# Dolphinx specifics
set(PROJECT_NAME photoelasticity)
project(${PROJECT_NAME} LANGUAGES C CXX)

# Set C++20 standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT TARGET dolfinx)
  find_package(DOLFINX REQUIRED)
endif()

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(MPI REQUIRED)

include(CheckSymbolExists)
set(CMAKE_REQUIRED_INCLUDES ${PETSC_INCLUDE_DIRS})
check_symbol_exists(PETSC_USE_COMPLEX petscsystypes.h PETSC_SCALAR_COMPLEX)
check_symbol_exists(PETSC_USE_REAL_DOUBLE petscsystypes.h PETSC_REAL_DOUBLE)

# Add target to compile UFL files
if(PETSC_SCALAR_COMPLEX EQUAL 1)
  if(PETSC_REAL_DOUBLE EQUAL 1)
    set(SCALAR_TYPE "--scalar_type=complex128")
  else()
    set(SCALAR_TYPE "--scalar_type=complex64")
  endif()
else()
  if(PETSC_REAL_DOUBLE EQUAL 1)
    set(SCALAR_TYPE "--scalar_type=float64")
  else()
    set(SCALAR_TYPE "--scalar_type=float32")
  endif()
endif()

add_custom_command(
  OUTPUT sphere.c
  COMMAND ffcx ${CMAKE_CURRENT_SOURCE_DIR}/sphere.py ${SCALAR_TYPE}
  VERBATIM
  DEPENDS sphere.py
  COMMENT "Compile sphere.py using FFCx"
)

# Generate the fem_template files
set (max_contacts "7" CACHE STRING "Maximum numbre of contact per grain")
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/fem_template_generator.sh ${CMAKE_CURRENT_BINARY_DIR}/fem_template_generator.sh USE_SOURCE_PERMISSIONS)
add_custom_command(
  OUTPUT fem_template
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/fem_template_generator.sh ${max_contacts} ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM
  DEPENDS fem_impl_2.cpp fem_template_generator.sh
  COMMENT "Generate the FEM templates"
)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wextra -pedantic  ")
# -L/home/fguillard/miniconda3/envs/fenicsx-env2/lib -L/usr/lib/x86_64-linux-gnu
add_executable(${PROJECT_NAME} Main.cpp Image.cpp Geometry.cpp Fem.cpp Cells.cpp ${CMAKE_CURRENT_BINARY_DIR}/sphere.c fem_template)
target_link_libraries(${PROJECT_NAME} dolfinx SDL2 pthread Eigen3::Eigen mpi)




