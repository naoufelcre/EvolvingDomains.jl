cmake_minimum_required(VERSION 3.21)
project(EvolvingDomainsCutCells LANGUAGES C CXX)

include(FetchContent)

# The package ships one shared C ABI library with CutCells linked statically.
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(MSVC)
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

FetchContent_Declare(cutcells
  GIT_REPOSITORY https://github.com/sclaus2/CutCells.git
  GIT_TAG ed5232ecb6858ae31729a907d66bb8f2cc8b4c56
  GIT_SHALLOW FALSE
  SOURCE_SUBDIR cpp)
FetchContent_MakeAvailable(cutcells)
set_target_properties(cutcells PROPERTIES
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN YES)

add_library(evolvingdomains_cut_cells SHARED cutcells_compat.cpp)
target_link_libraries(evolvingdomains_cut_cells PRIVATE cutcells)
target_compile_features(evolvingdomains_cut_cells PRIVATE cxx_std_20)
if(UNIX AND NOT APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  target_link_options(evolvingdomains_cut_cells PRIVATE
    -static-libstdc++ -static-libgcc -Wl,--exclude-libs,ALL)
endif()
set_target_properties(evolvingdomains_cut_cells PROPERTIES
  OUTPUT_NAME evolvingdomains_cut_cells
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN YES)
if(WIN32)
  set_target_properties(evolvingdomains_cut_cells PROPERTIES PREFIX "")
endif()

install(TARGETS evolvingdomains_cut_cells
  RUNTIME DESTINATION lib COMPONENT RuntimeLibraries
  LIBRARY DESTINATION lib COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION lib COMPONENT RuntimeLibraries)
install(FILES
  "${cutcells_SOURCE_DIR}/LICENSE"
  "${cutcells_SOURCE_DIR}/third_party/THIRD_PARTY_NOTICES.md"
  "${cutcells_SOURCE_DIR}/third_party/VTK-Copyright.txt"
  DESTINATION share/licenses/cutcells COMPONENT RuntimeLibraries)
