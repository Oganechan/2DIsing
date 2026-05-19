add_library(optimization_flags INTERFACE)
target_compile_options(optimization_flags INTERFACE
    $<$<CONFIG:Release>:-O3 -march=native -ffast-math -funroll-loops>
)

add_library(quadmath INTERFACE)
target_link_libraries(quadmath INTERFACE -lquadmath)

add_library(lattice STATIC)
target_sources(lattice
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/lattice/topology/topology.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/lattice/interaction/interaction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/lattice/atoms.cpp
)
target_include_directories(lattice
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/lattice
    ${CMAKE_CURRENT_SOURCE_DIR}/lattice/topology
    ${CMAKE_CURRENT_SOURCE_DIR}/lattice/interaction
)
target_link_libraries(lattice
  PUBLIC
    config
    random
  PRIVATE
    optimization_flags
)

add_library(physics STATIC)
target_sources(physics
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/physics/calculator.cpp
)
target_include_directories(physics
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/physics
)
target_link_libraries(physics
  PUBLIC
    lattice
  PRIVATE
    optimization_flags
)

add_library(algorithm STATIC)
target_sources(algorithm
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/algorithm/wanglandau.cpp
)
target_include_directories(algorithm
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/algorithm
)
target_link_libraries(algorithm
  PUBLIC
    physics
    quadmath
  PRIVATE
    optimization_flags
)

add_library(simulation STATIC)
target_sources(simulation
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/simulation.cpp
)
target_include_directories(simulation
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(simulation
    PUBLIC
      algorithm
    PRIVATE
      optimization_flags
      OpenMP::OpenMP_CXX
)
