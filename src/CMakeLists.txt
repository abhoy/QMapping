# cmake-lint: disable=C0103
set(JSON_BuildTests
    OFF
    CACHE INTERNAL "")
# cmake-lint: disable=C0103
set(JSON_MultipleHeaders
    OFF
    CACHE INTERNAL "")
# add nlohmann::json library
add_subdirectory("${PROJECT_SOURCE_DIR}/extern/json" "extern/json" EXCLUDE_FROM_ALL)

# main project library
add_library(
  ${PROJECT_NAME} SHARED
  ${PROJECT_SOURCE_DIR}/include/graph.hpp
  ${PROJECT_SOURCE_DIR}/include/gate.hpp
  ${PROJECT_SOURCE_DIR}/include/circuit.hpp
  ${PROJECT_SOURCE_DIR}/include/costmetrics.hpp
  ${PROJECT_SOURCE_DIR}/include/initialmapping.hpp
  ${PROJECT_SOURCE_DIR}/include/localmapping.hpp
  circuit.cpp
  gate.cpp
  graph.cpp
  costmetrics.cpp
  initialmapping.cpp
  localmapping.cpp
)

set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

# set include directories
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include
                                                  ${PROJECT_BINARY_DIR}/include)

target_link_libraries(${PROJECT_NAME} PUBLIC nlohmann_json)

SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)