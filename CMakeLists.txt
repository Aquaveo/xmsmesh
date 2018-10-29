set(CMAKE_CXX_STANDARD 11)
cmake_minimum_required(VERSION 3.1.2)
cmake_policy(SET CMP0015 NEW) # Link Directory Pathing
set(CMAKE_DEBUG_POSTFIX _d)

project(xmsmesh C CXX)

if (APPLE)  
  set(CMAKE_POSITION_INDEPENDENT_CODE False)
else()
  set(CMAKE_POSITION_INDEPENDENT_CODE True)
endif()

set(BUILD_TESTING NO CACHE BOOL "Enable/Disable testing")
set(IS_CONDA_BUILD NO CACHE BOOL "Set this if you want to make a conda package.")
set(CONDA_PREFIX "" CACHE PATH "Path to the conda environment used to build.")
set(IS_PYTHON_BUILD NO CACHE BOOL "Set this if you want to build the python bindings.")
set(XMS_TEST_PATH ${PROJECT_SOURCE_DIR}/test_files/ CACHE PATH "Path to test files for testing")
set(PYTHON_TARGET_VERSION 3.6 CACHE STRING "Version of python to link to for python wrapping.")
set(XMS_VERSION "\"99.99.99\"" CACHE STRING "Library Version")

add_definitions(-DXMS_VERSION=\"${XMS_VERSION}\")

if(WIN32)
    if(XMS_BUILD)
        add_definitions(/D _WIN32_WINNT=0x0501)  # Windows XP and higher
        add_definitions(/Zc:wchar_t-)  # Treat wchar_t as built-in type
    else(NOT XMS_BUILD)
        add_definitions(/D BOOST_ALL_NO_LIB)
    endif()
endif()

if(IS_CONDA_BUILD)
  include(${CMAKE_CURRENT_LIST_DIR}/condabuildinfo.cmake)
else() # If we are not using conda, we are using conan
  # Conan setup
  include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
  conan_basic_setup(TARGETS)
  set(EXT_INCLUDE_DIRS ${CONAN_INCLUDE_DIRS})
  set(EXT_LIB_DIRS ${CONAN_LIB_DIRS})
  set(EXT_LIBS ${CONAN_LIBS})
endif(IS_CONDA_BUILD)

if(IS_PYTHON_BUILD)
    # linux and mac builds for conan (on TRAVISCI) use a docker that has python
    # 2.7 as system python. We do not have control over that docker image so we
    # can't change this and it is fine for building conan packages and checking
    # that we don't have errors in the python wrapping. We have conda recipes
    # for building python packages that target other versions of python.
    find_package(PythonLibs ${PYTHON_TARGET_VERSION} EXACT REQUIRED)

    message("PYTHON_INCLUDE: ${PYTHON_INCLUDE_DIRS}")
    message("PYTHON_LIBS: ${PYTHON_LIBRARIES}")

    # Pybind11 module
    if(IS_CONDA_BUILD)
        include("${CONDA_PREFIX}/share/cmake/pybind11/pybind11Targets.cmake")
        include("${CONDA_PREFIX}/share/cmake/pybind11/FindPythonLibsNew.cmake")
        include("${CONDA_PREFIX}/share/cmake/pybind11/pybind11Config.cmake")
        include("${CONDA_PREFIX}/share/cmake/pybind11/pybind11ConfigVersion.cmake")
        include("${CONDA_PREFIX}/share/cmake/pybind11/pybind11Tools.cmake")

    else()
        include("${CONAN_PYBIND11_ROOT}/share/cmake/pybind11/pybind11Targets.cmake")
        include("${CONAN_PYBIND11_ROOT}/share/cmake/pybind11/FindPythonLibsNew.cmake")
        include("${CONAN_PYBIND11_ROOT}/share/cmake/pybind11/pybind11Config.cmake")
        include("${CONAN_PYBIND11_ROOT}/share/cmake/pybind11/pybind11ConfigVersion.cmake")
        include("${CONAN_PYBIND11_ROOT}/share/cmake/pybind11/pybind11Tools.cmake")
    endif()
    list(APPEND EXT_INCLUDE_DIRS
      ${PYTHON_INCLUDE_DIRS}
    )
endif()

message("External Include Dirs: ${EXT_INCLUDE_DIRS}")
message("External Lib Dirs: ${EXT_LIB_DIRS}")
message("External Libs: ${EXT_LIBS}")

include_directories(${CMAKE_CURRENT_LIST_DIR})
include_directories(${EXT_INCLUDE_DIRS})
link_directories(${EXT_LIB_DIRS})

# Static library sources
set(xmsmesh_sources
  xmsmesh/meshing/MeMeshUtils.cpp
  xmsmesh/meshing/MeMultiPolyTo2dm.cpp
  xmsmesh/meshing/MeMultiPolyMesher.cpp
  xmsmesh/meshing/detail/MeBadQuadRemover.cpp
  xmsmesh/meshing/detail/MeIntersectPolys.cpp
  xmsmesh/meshing/detail/MePolyPatcher.cpp
  xmsmesh/meshing/detail/MePolyOffsetter.cpp
  xmsmesh/meshing/detail/MePolyPaverToMeshPts.cpp
  xmsmesh/meshing/detail/MePolyCleaner.cpp
  xmsmesh/meshing/detail/MePolyPts.cpp
  xmsmesh/meshing/detail/MePolyRedistributePtsCurvature.cpp
  xmsmesh/meshing/detail/MeQuadBlossom.cpp
  xmsmesh/meshing/detail/MeRefinePtsToPolys.cpp
  xmsmesh/meshing/detail/MeRelaxer.cpp
  xmsmesh/meshing/detail/MeWeightMatcher.cpp
  xmsmesh/meshing/MePolyMesher.cpp
  xmsmesh/meshing/MePolyRedistributePts.cpp
  xmsmesh/tutorial/TutMeshing.cpp
)

set(xmsmesh_headers
  xmsmesh/meshing/MeMeshUtils.h
  xmsmesh/meshing/MePolyMesher.h
  xmsmesh/meshing/MeMultiPolyMesher.h
  xmsmesh/meshing/MeMultiPolyMesherIo.h
  xmsmesh/meshing/MeMultiPolyTo2dm.h
  xmsmesh/meshing/MePolyRedistributePts.h
  xmsmesh/meshing/detail/MeBadQuadRemover.h
  xmsmesh/meshing/detail/MePolyCleaner.h
  xmsmesh/meshing/detail/MePolyOffsetter.h
  xmsmesh/meshing/detail/MePolyPts.h
  xmsmesh/meshing/detail/MePolyPatcher.h
  xmsmesh/meshing/detail/MeIntersectPolys.h
  xmsmesh/meshing/detail/MePolyPaverToMeshPts.h
  xmsmesh/meshing/detail/MePolyRedistributePtsCurvature.h
  xmsmesh/meshing/detail/MeQuadBlossom.h
  xmsmesh/meshing/detail/MeRefinePtsToPolys.h
  xmsmesh/meshing/detail/MeRelaxer.h
  xmsmesh/meshing/detail/MeWeightMatcher.h
)

# Pybind11 sources
set(xmsmesh_py_source
  # xmsmesh/python/cxxtest.cpp
  xmsmesh/python/xmsmesh_py.cpp
  # Meshing
  xmsmesh/python/meshing/meshing_py.h
  xmsmesh/python/meshing/meshing_py.cpp
  xmsmesh/python/meshing/MeBadQuadRemover_py.cpp
  xmsmesh/python/meshing/MeMeshUtils_py.cpp
  xmsmesh/python/meshing/MeMultiPolyMesher_py.cpp
  xmsmesh/python/meshing/MeMultiPolyMesherIo_py.cpp
  xmsmesh/python/meshing/MeMultiPolyTo2dm_py.cpp
  xmsmesh/python/meshing/MePolyMesher_py.cpp
  xmsmesh/python/meshing/MePolyRedistributePts_py.cpp
  xmsmesh/python/meshing/MeQuadBlossom_py.cpp
)

# Tests
if (BUILD_TESTING)
  add_definitions(-DXMS_TEST_PATH="${XMS_TEST_PATH}/")
  add_definitions(-DCXX_TEST -DCXXTEST4)

  list(APPEND xmsmesh_sources
    xmsmesh/meshing/MeMeshUtils.t.h
    xmsmesh/meshing/MeMultiPolyTo2dm.t.h
    xmsmesh/meshing/MePolyMesher.t.h
    xmsmesh/meshing/MeMultiPolyMesher.t.h
    xmsmesh/meshing/MePolyRedistributePts.t.h
    xmsmesh/meshing/detail/MeBadQuadRemover.t.h
    xmsmesh/meshing/detail/MePolyPaverToMeshPts.t.h
    xmsmesh/meshing/detail/MeIntersectPolys.t.h
    xmsmesh/meshing/detail/MePolyPatcher.t.h
    xmsmesh/meshing/detail/MePolyOffsetter.t.h
    xmsmesh/meshing/detail/MePolyCleaner.t.h
    xmsmesh/meshing/detail/MePolyRedistributePtsCurvature.t.h
    xmsmesh/meshing/detail/MeQuadBlossom.t.h
    xmsmesh/meshing/detail/MeRefinePtsToPolys.t.h
    xmsmesh/meshing/detail/MeRelaxer.t.h
    xmsmesh/meshing/detail/MeWeightMatcher.t.h
    xmsmesh/tutorial/TutMeshing.t.h
  )

  find_package(CxxTest)
  if(CXXTEST_FOUND)
    include_directories(${CXXTEST_INCLUDE_DIRS})
    enable_testing()

    set(CXXTEST_TESTGEN_ARGS --xunit-printer --have-eh)
    file(GLOB_RECURSE test_headers ${CMAKE_CURRENT_LIST_DIR}/xmsmesh/*.t.h)
    CXXTEST_ADD_TEST(
      runner runner.cpp ${test_headers}
    )
    target_link_libraries(runner ${PROJECT_NAME})
  endif()
endif ()

# Static library
add_library(${PROJECT_NAME} STATIC
  ${xmsmesh_sources} ${xmsmesh_headers}
)
target_include_directories(${PROJECT_NAME}
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/>
    ${Boost_INCLUDE_DIR}
 )
find_package(Threads REQUIRED)
target_link_libraries(${PROJECT_NAME}
  ${EXT_LIBS}
)
target_link_libraries(${PROJECT_NAME}
  ${CMAKE_THREAD_LIBS_INIT}
)
if(UNIX AND NOT APPLE)
    target_link_libraries(${PROJECT_NAME} rt)
endif()

if(IS_PYTHON_BUILD)
    pybind11_add_module(xmsmesh_py
      ${xmsmesh_py_source}
    )

    target_include_directories(xmsmesh_py
        PRIVATE
          ${EXT_LIB_DIRS}
          ${PYTHON_INCLUDE_DIRS}
    )
    target_link_libraries(xmsmesh_py
      PRIVATE
        ${EXT_LIBS}
        ${PROJECT_NAME}
    )
    set_target_properties(xmsmesh_py PROPERTIES
      LINKER_LANGUAGE CXX
    )

    # Install recipe
    install(
      TARGETS xmsmesh_py
      ARCHIVE DESTINATION "site-packages"
      LIBRARY DESTINATION "site-packages"
    )

endif()


# Install recipe
install(
  TARGETS ${PROJECT_NAME}
  ARCHIVE DESTINATION "lib"
  LIBRARY DESTINATION "lib"
)
foreach (header IN LISTS xmsmesh_headers)
  get_filename_component(subdir "${header}" DIRECTORY)
  install(
    FILES "${header}"
    DESTINATION "include/${subdir}"
  )
endforeach ()
