cmake_minimum_required(VERSION 3.10)

project(
    QMapping 
    VERSION 1.0
    LANGUAGES CXX
    DESCRIPTION "QCORE: A library for core quantum gate and circuit representation")

# check whether the submodule ``modulename`` is correctly cloned in the ``/extern`` directory.
macro(CHECK_SUBMODULE_PRESENT modulename)
  if(NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/${modulename}/CMakeLists.txt")
    message(
      FATAL_ERROR
        "${modulename} submodule not cloned properly. \
        Please run `git submodule update --init --recursive` \
        from the main project directory")
  endif()
endmacro()

check_submodule_present(json)

# add main library code
add_subdirectory(src)

#add executable 
add_subdirectory(test)

# add Python binding code
option(BINDINGS "Configure for building Python bindings" OFF)
if(BINDINGS)
    check_submodule_present(pybind11)
    check_submodule_present(pybind11_json)
    add_subdirectory(QMapping)
endif()
