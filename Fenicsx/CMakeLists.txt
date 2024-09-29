cmake_minimum_required(VERSION 3.5)
set(PROJECT_NAME sphere)
project(${PROJECT_NAME})
# Set CMake behavior
cmake_policy(SET CMP0004 NEW)
# Get DOLFIN configuration data
find_package(DOLFINX REQUIRED)

# Add executable
add_executable(${PROJECT_NAME} sphere.c)
# Target libraries
target_link_libraries(${PROJECT_NAME} ${DOLFIN_LIBRARIES})

