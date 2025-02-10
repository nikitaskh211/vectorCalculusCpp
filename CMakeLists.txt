# writing minimal requirements
CMAKE_MINIMUM_REQUIRED(VERSION 3.30)

# making a project
PROJECT(MoiseenkoMaxwell VERSION 0.0.1 LANGUAGES CXX)

# message about project's metadata
MESSAGE(STATUS "[R/T] Project name : ${PROJECT_NAME}")
MESSAGE(STATUS "[R/T] Project version : ${PROJECT_VERSION}")

# message about the standard
MESSAGE(STATUS "[R/T] Setting the minimal C++ standard.")

# setting a standard
SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_REQUIRED ON)
SET(CMAKE_CXX_EXTENSIONS OFF)

# message about C++17 standard
MESSAGE(STATUS "[R/T] C++ standard : " ${CMAKE_CXX_STANDARD})

# message about llvm
MESSAGE(STATUS "[R/T] Setting the LLVM toolchain.")

# llvm toolchain setup
SET(CMAKE_C_COMPILER clang)
SET(CMAKE_CXX_COMPILER clang++)
SET(CMAKE_LINKER lld-link)

# message about compiler
MESSAGE(STATUS "[R/T] C compiler : " ${CMAKE_C_COMPILER})
MESSAGE(STATUS "[R/T] C++ compiler : " ${CMAKE_CXX_COMPILER})
MESSAGE(STATUS "[R/T] Linker : " ${CMAKE_LINKER})

# including subdirectories
ADD_SUBDIRECTORY(lib)
ADD_SUBDIRECTORY(src)

# message about libraries
MESSAGE(STATUS "[R/T] Subdirectories for libraries : " ${CMAKE_SOURCE_DIR}/lib)
MESSAGE(STATUS "[R/T] Subdirectories for sources : " ${CMAKE_SOURCE_DIR}/src)

# setting a path for output executable
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

# message about output binary
MESSAGE(STATUS "[R/T] Returning output to : " ${CMAKE_SOURCE_DIR}/bin)
