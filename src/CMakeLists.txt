# message about current directory instructions
MESSAGE(STATUS "[SRC] Adding the target executable for unit testing.")

# creating a target executable
ADD_EXECUTABLE(unitTest unitTesting.cpp)

# message about executable name
MESSAGE(STATUS "[SRC] Target executable for unit testing : unitTest")

# message about linking the libraries to the executable
MESSAGE(STATUS "[SRC] Linking libraries to the executable.")

# linking the libraries to the executable
TARGET_LINK_LIBRARIES(unitTest vecLib)
