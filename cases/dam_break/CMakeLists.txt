cmake_minimum_required(VERSION 3.16)
string(REGEX REPLACE ".*/(.*)" "\\1" CURRENT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR})
project(${CURRENT_FOLDER} LANGUAGES CXX)

find_package(glm REQUIRED)

# set
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

# file
file(GLOB_RECURSE CPP_FILES "particle_generate.cpp")

# include
include_directories(${PROJECT_SOURCE_DIR}/../../include)

include_directories(${GLM_INCLUDE_DIRS})

# target
add_executable(${PROJECT_NAME} ${CPP_FILES})
# link
target_link_libraries(${PROJECT_NAME} ${GLM_LIBRARIES})