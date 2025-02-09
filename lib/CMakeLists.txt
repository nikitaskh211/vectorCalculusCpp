# message about the begin of the linking process
MESSAGE(STATUS "[LIB] Started library linking.")

# library linking
ADD_LIBRARY(vecLib vec.cpp alias.cpp)

# message about the end of the linking process
MESSAGE(STATUS "[LIB] Ended library linking.")

# message about the begin of the linking process
MESSAGE(STATUS "[LIB] Started header linking.")

# header->library linking
TARGET_INCLUDE_DIRECTORIES(vecLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/inc)

# message about the end of the linking process
MESSAGE(STATUS "[LIB] Ended header linking.")
