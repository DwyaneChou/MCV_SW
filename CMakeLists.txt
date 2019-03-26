cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

project(tmcore LANGUAGES Fortran)

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Ofast")
  if (CMAKE_BUILD_TYPE MATCHES DEBUG)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -gdwarf -Og -traceback -check all")
  endif ()
elseif (CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  if (APPLE)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Ofast -floop-parallelize-all -fdefault-real-8 -ffree-line-length-none")
  else ()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Ofast -floop-parallelize-all -fdefault-real-8 -ftree-parallelize-loops=8 -ffree-line-length-none")
  endif ()
  if (CMAKE_BUILD_TYPE MATCHES DEBUG)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -gdwarf -Og -fbacktrace -fcheck=all")
  endif ()
endif ()

if (DEFINED ENV{NETCDF} AND (NOT DEFINED ENV{NETCDF_ROOT}))
  set(ENV{NETCDF_ROOT} $ENV{NETCDF})
endif ()
if (DEFINED ENV{NETCDF_ROOT})
  include_directories("$ENV{NETCDF_ROOT}/include")
  link_directories("$ENV{NETCDF_ROOT}/lib")
else ()
  find_package(PkgConfig)
  if (PKG_CONFIG_FOUND)
    pkg_search_module(NETCDF REQUIRED netcdf)
    if (NETCDF_FOUND)
      include_directories(${NETCDF_INCLUDEDIR})
      link_directories(${NETCDF_LIBDIR})
    else ()
      message(FATAL_ERROR "Pkg-config could not find netcdf library!")
    endif ()
  else ()
    message(FATAL_ERROR "Unable to find pkg-config library!")
  endif ()
endif ()

if (DEFINED ENV{LAPACK_ROOT})
  include_directories("$ENV{LAPACK_ROOT}")
  link_directories("$ENV{LAPACK_ROOT}")
  set(LAPACK_LINKER_FLAGS "-llapack")
else ()
  find_package(LAPACK REQUIRED)
  if (NOT LAPACK_FOUND)
    message(DEBUG "Unable to find LAPACK library, use internal routines.")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DNOLAPACK")
  endif ()
endif ()

add_subdirectory(lib/container)
include_directories(${CMAKE_BINARY_DIR}/fortran_container)

add_subdirectory(lib/datetime)

set(source_files
  src/utils/const_mod.F90
  src/utils/string_mod.F90
  src/utils/time_mod.F90
  src/utils/timer_mod.F90
  src/utils/log_mod.F90
  src/utils/io_mod.F90
  src/utils/history_mod.F90
  src/utils/sort_mod.F90
  src/utils/sphere_geometry_mod.F90
  src/utils/math_mod.F90
  src/utils/poly_fit_mod.F90
  src/utils/random_number_mod.F90
  src/utils/delaunay_voronoi_mod.F90
  src/core/params_mod.F90
  src/core/mesh_mod.F90
  src/core/static_mod.F90
  src/core/state_mod.F90
  src/core/tend_mod.F90
  src/core/operators_mod.F90
  src/core/adv_schemes/adv_scheme_sg11_mod.F90
  src/core/adv_schemes/adv_scheme_mod.F90
  src/core/diag_mod.F90
  src/core/time_scheme_mod.F90
  src/core/tmcore_adv_mod.F90
  src/core/tmcore_swm_mod.F90
  src/test_cases/test_cases_mod.F90
  src/test_cases/swm/mountain_zonal_flow_test_mod.F90
  src/test_cases/swm/steady_geostrophic_flow_test_mod.F90
  src/test_cases/swm/rossby_haurwitz_wave_test_mod.F90
  src/test_cases/swm/linearized_rossby_wave_test_mod.F90
  src/test_cases/adv/cosine_bell_rotation_test_mod.F90
  src/test_cases/adv/deformation_flow_test_mod.F90
)

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/.mods")

add_library(tmcore ${source_files})
target_link_libraries(tmcore fortran_datetime fortran_container netcdf netcdff ${LAPACK_LIBRARIES})

add_executable(test_tmcore_swm.exe src/core/test_tmcore_swm.F90)
target_link_libraries(test_tmcore_swm.exe tmcore)

add_executable(test_tmcore_adv.exe src/core/test_tmcore_adv.F90)
target_link_libraries(test_tmcore_adv.exe tmcore)

add_executable(gen_mesh.exe src/tools/gen_mesh.F90)
target_link_libraries(gen_mesh.exe tmcore)
