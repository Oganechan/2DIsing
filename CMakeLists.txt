cmake_minimum_required(VERSION 4.0)
project(Rectangular2DIsing)

find_package(OpenMP REQUIRED)

add_subdirectory(external)
add_subdirectory(src)
add_subdirectory(tests)

add_executable(${PROJECT_NAME})
target_sources(${PROJECT_NAME}
  PRIVATE
    src/main.cpp
)
target_link_libraries(${PROJECT_NAME}
  PUBLIC
    simulation
  PRIVATE
    optimization_flags
    quadmath
    OpenMP::OpenMP_CXX
 
)

install(TARGETS Rectangular2DIsing
  RUNTIME DESTINATION bin
)

install(TARGETS simulation algorithm physics lattice
  ARCHIVE DESTINATION lib
)

install(TARGETS test_geometry test_atoms test_calculator test_wanglandau
  RUNTIME DESTINATION tests
)

install(PROGRAMS 
  scripts/visualize_lattice.py
  scripts/visualize_ln_g.py
  scripts/visualize_plots.py
  scripts/visualize_Tc.py
  scripts/auto_Tc.sh
  scripts/run.sh

  DESTINATION scripts
)

install(FILES config/config.json
  DESTINATION bin
)