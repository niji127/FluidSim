cmake_minimum_required(VERSION 3.26)
string(REGEX REPLACE ".*/(.*)" "\\1" CURRENT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR})
project(${CURRENT_FOLDER} LANGUAGES CUDA CXX)

find_package(CUDA REQUIRED)
find_package(Thrust REQUIRED CONFIG)
find_package(glm REQUIRED)
find_package(Boost REQUIRED COMPONENTS filesystem system)

include(CheckLanguage)
check_language(CUDA)

# file
file(GLOB_RECURSE CUDA_FILES "src/*.cu")
file(GLOB_RECURSE CPP_FILES "src/*.cpp")

# Set compiler options
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_90")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored")
set(CMAKE_CUDA_FLAGS_RELEASE "${CMAKE_CUDA_FLAGS} -Xptxas -O3,-dlcm=cg")
set(CMAKE_CUDA_FLAGS_RELEASE "${CMAKE_CUDA_FLAGS_RELEASE} -use_fast_math")
set(CMAKE_CUDA_FLAGS_RELEASE "${CMAKE_CUDA_FLAGS_RELEASE} -Xcompiler /O2")

# include
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${CUDA_INCLUDE_DIRS})
include_directories(${GLM_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

link_directories(${Boost_LIBRARY_DIRS})

# Add the executable
add_executable(${PROJECT_NAME} ${CPP_FILES} ${CUDA_FILES})

thrust_create_target(Thrust)
target_link_libraries(${PROJECT_NAME} Thrust)
target_link_libraries(${PROJECT_NAME} ${GLM_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})

# set executable path
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin) 
# Set CUDA specific properties
set_target_properties(${PROJECT_NAME} PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_STANDARD 17
    CUDA_STANDARD_REQUIRED ON
    CXX_STANDARD 17
)