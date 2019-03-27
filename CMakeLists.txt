cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

project(MCV_SW LANGUAGES Fortran)

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

add_subdirectory(lib/fortran-container)

add_subdirectory(lib/fortran-datetime)

set(source_files
  src/MCV_SW.f90
  src/dynamics/initialize_mod.f90
  src/share/constants_mod.f90
  src/share/domain_mod.f90
  src/share/parameters_mod.f90
  src/share/projection_mod.f90
)

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/.mods")

add_library(MCV_SW ${source_files})
target_link_libraries(MCV_SW fortran_datetime fortran_container netcdf netcdff ${LAPACK_LIBRARIES})

add_executable(MCV_SW.exe src/MCV_SW.f90)
target_link_libraries(MCV_SW.exe MCV_SW)

