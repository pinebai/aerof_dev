cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
enable_language(Fortran)

SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_admin")

###################
# Local Variables #
###################
SITE_NAME(HOSTNAME)
FIND_PACKAGE(LOCAL NAMES ${HOSTNAME} PATHS . ./localconfig QUIET)
MARK_AS_ADVANCED(HOSTNAME)
MARK_AS_ADVANCED(LOCAL_DIR)
###################

option(AEROACOUSTIC "Enable aeroacoustic computations" OFF)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-write-strings -Wno-unused-result -Wno-aggressive-loop-optimizations -Wno-format-zero-length")
endif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")

IF(USE_OPENMP)
  FIND_PACKAGE(OpenMP REQUIRED)
  if(OpenMP_FLAG_DETECTED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_CXX_FLAGS}")
    #add_definitions(-D_OPENMP)
  endif(OpenMP_FLAG_DETECTED)
ENDIF(USE_OPENMP)

option(USE_STRONG_FARFIELD_BC "Enable strong farfield bc" OFF)

IF(USE_STRONG_FARFIELD_BC)
  add_definitions(-DSTRONG_FARFIELD_BC)
ENDIF(USE_STRONG_FARFIELD_BC)

IF(YDEBUG EQUAL 1)
  add_definitions(-DYDEBUG)
endif(YDEBUG EQUAL 1)

# Code to determine if we should add some link libraries
IF("${CMAKE_Fortran_COMPILER_ID}" MATCHES GNU)
  EXEC_PROGRAM(${CMAKE_Fortran_COMPILER}
      ARGS "-### CMake/Modules/CMakeTestFortran.f"
      OUTPUT_VARIABLE CMAKE_COMPILER_OUTPUT RETURN_VALUE CMAKE_COMPILER_RETURN)
  IF(NOT CMAKE_COMPILER_RETURN)
#      MESSAGE("Will add -lgfortran")
      IF("${CMAKE_COMPILER_OUTPUT}" MATCHES ".*-lgfortran.*" )
         SET(FORTRAN_HYBRID_LIBRARIES "-lgfortran")
      ENDIF("${CMAKE_COMPILER_OUTPUT}" MATCHES ".*-lgfortran.*" )
#  ELSE(NOT CMAKE_COMPILER_RETURN)
#     MESSAGE("Problem: " ${CMAKE_COMPILER_OUTPUT} )
  ENDIF(NOT CMAKE_COMPILER_RETURN)
ENDIF("${CMAKE_Fortran_COMPILER_ID}" MATCHES GNU)

# SET(CMAKE_VERBOSE_MAKEFILE true)

IF(NOT EXTRALIB)
SET(EXTRALIB "")
ENDIF(NOT EXTRALIB)

IF(AEROACOUSTIC)
  FIND_PACKAGE(GSL REQUIRED)
  FIND_PACKAGE(FFTW REQUIRED)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GSL_CFLAGS} -I${FFTW_INCLUDE_DIR} -DAEROACOUSTIC")
ENDIF(AEROACOUSTIC)

#message("${GSL_LIBRARIES} ${FFTW_LIBRARIES}")

IF(NOT WITH_IPOPT EQUAL 0)
  FIND_PACKAGE(IPOPT)
  if(IPOPT_FOUND)
    include_directories(${IPOPT_INCLUDE_DIRS})
    SET(EXTRALIBRARIES ${IPOPT_LIBRARIES} ${EXTRALIBRARIES})
  endif(IPOPT_FOUND)
ENDIF(NOT WITH_IPOPT EQUAL 0)

SET(EXTRALIB ${EXTRALIB} CACHE STRING "Extra link parameters")

SET(CMAKE_DEBUG_POSTFIX "_debug")

if(NOT LAPACK_FOUND)
 FIND_PACKAGE(LAPACK)
endif(NOT LAPACK_FOUND)
IF(BLAS_LIBRARIES)
  message("BLAS Libraries " ${BLAS_LIBRARIES} " found")
  SET(EXTRALIBRARIES ${BLAS_LIBRARIES})
ELSE(BLAS_LIBRARIES)
  SET(EXTRALIBRARIES "")
ENDIF(BLAS_LIBRARIES)
IF(LAPACK_LIBRARIES)
  message("LAPCK Libraries " ${LAPACK_LIBRARIES} " found")
  SET(EXTRALIBRARIES ${LAPACK_LIBRARIES} ${EXTRALIBRARIES})
ENDIF(LAPACK_LIBRARIES)

#message("${GSL_LIBRARIES} ${FFTW_LIBRARIES}")
if (AEROACOUSTIC)
  SET(EXTRALIB ${EXTRALIB} ${GSL_LIBRARIES} ${FFTW_LIBRARIES})
ENDIF(AEROACOUSTIC)

FIND_PACKAGE(MPI)

SET(TARGETS_INST fluid)

IF(UNIX)
  MESSAGE("-- Building for system type: " ${CMAKE_SYSTEM_NAME} ".")
  IF(${CMAKE_SYSTEM_NAME} STREQUAL Linux)
    add_definitions(-D__LINUX)
  ENDIF(${CMAKE_SYSTEM_NAME} STREQUAL Linux)

  IF(${CMAKE_SYSTEM_NAME} STREQUAL Linux OR ${CMAKE_SYSTEM_NAME} STREQUAL Darwin)
    IF(NOT WITH_ARPACK EQUAL 0)
      FIND_LIBRARY(ARPACKLIB NAMES arpack arpack_LINUX
         PATHS  . /usr/lib /usr/local/lib)
    ELSE(NOT WITH_ARPACK EQUAL 0)
      SET(ARPACKLIB "ARPACKLIB-NOTFOUND")
    ENDIF(NOT WITH_ARPACK EQUAL 0)

    SET(DO_MODAL 0)
    IF(ARPACKLIB AND LAPACK_FOUND)
      MESSAGE("-- ARPACK library " ${ARPACKLIB} " found.")
      add_definitions(-DDO_MODAL)
      SET(EXTRALIBRARIES 
          ${ARPACKLIB}
          ${EXTRALIBRARIES})
      SET(DO_MODAL 1)
    ENDIF(ARPACKLIB AND LAPACK_FOUND)

    SET(CMAKE_FIND_LIBRARY_PREFIXES "lib" "")
    IF(NOT WITH_SCALAPACK EQUAL 0)
      IF(NOT SCALAPACK_FOUND)
        FIND_LIBRARY(SCALAPACKLIB NAMES scalapack scalapack-openmpi
                     PATHS  . /usr/lib /usr/local/lib /home/hg/lib/SCALAPACK ${SCALAPACKDIR})
        IF(SCALAPACKLIB)
          MESSAGE("-- SCALAPACK library " ${SCALAPACKLIB} " found.")
          SET(SCALAPACK_FOUND TRUE)
        ENDIF(SCALAPACKLIB)
      ENDIF(NOT SCALAPACK_FOUND)
    ELSE(NOT WITH_SCALAPACK EQUAL 0)
      SET(SCALAPACKLIB "SCALAPACKLIB-NOTFOUND")
      SET(SCALAPACK_FOUND FALSE)
    ENDIF(NOT WITH_SCALAPACK EQUAL 0)

    IF(NOT WITH_BLACS EQUAL 0)
      IF(NOT BLACS_FOUND)
        FIND_LIBRARY(BLACSLIB NAMES blacs blacs-openmpi
                     PATHS  . /usr/lib /usr/local/lib /home/hg/lib/SCALAPACK ${SCALAPACKDIR})
        FIND_LIBRARY(BLACSCLIB NAMES blacsC blacsCinit-openmpi
                     PATHS  . /usr/lib /usr/local/lib /home/hg/lib/SCALAPACK ${SCALAPACKDIR})
        IF(BLACSLIB AND BLACSCLIB)
          MESSAGE("-- BLACS library " ${BLACSLIB},${BLACSCLIB} " found.")
          SET(BLACS_FOUND TRUE)
        ENDIF(BLACSLIB AND BLACSCLIB)
      ENDIF(NOT BLACS_FOUND)
    ELSE(NOT WITH_BLACS EQUAL 0)
      SET(BLACSLIB "BLACSLIB-NOTFOUND")
      SET(BLACSCLIB "BLACSCLIB-NOTFOUND")
      SET(BLACSFOUND FALSE)
    ENDIF(NOT WITH_BLACS EQUAL 0)

    IF(SCALAPACK_FOUND)
      add_definitions(-DDO_SCALAPACK)
      IF(BLACS_FOUND)
        SET(EXTRALIBRARIES
          ${BLACSCLIB} ${BLACSLIB} 
          ${EXTRALIBRARIES})
      ENDIF(BLACS_FOUND)
      SET(EXTRALIBRARIES
        ${SCALAPACKLIB}
        ${EXTRALIBRARIES})
      SET(USE_SCALAPACK 1)
    ELSE(SCALAPACK_FOUND)
      SET(USE_SCALAPACK 0)
    ENDIF(SCALAPACK_FOUND)
    
  ENDIF(${CMAKE_SYSTEM_NAME} STREQUAL Linux OR ${CMAKE_SYSTEM_NAME} STREQUAL Darwin)
ENDIF(UNIX)


SET(BINDIR ${CMAKE_CURRENT_BINARY_DIR}/bin)

FIND_PROGRAM(MERCURIAL_EXECUTABLE hg)

SET(MAIN_FILES Main.C)

ADD_LIBRARY(FluidBase
             FSI/CrackingSurface.cpp
             FSI/DynamicNodalTransfer.cpp
             NavierStokesSolver.C 
             Solvers/LevelSetSolver5.C
             Solvers/NavierStokesCoupledSolver5.C
             Solvers/NavierStokesEmbeddedCoupledSolver5.C
             Solvers/NavierStokesEmbeddedCoupledSolver6.C
             Solvers/NavierStokesEmbeddedCoupledSolver7.C
             Solvers/NavierStokesEmbeddedSegSolver651.C
             Solvers/NavierStokesEmbeddedSegSolver752.C
             Solvers/NavierStokesMultiPhysicsEmbedded5.C
             Solvers/NavierStokesCoupledSolver6.C
             Solvers/NavierStokesCoupledSolver7.C
             Solvers/NavierStokesSegSolver651.C
             Solvers/NavierStokesSegSolver752.C
             AgglomeratedFace.C
             BcFcnCore.C
             BlockTridiagonalMatrix.C
             CommunicatorCore.C 
             CommunicatorMPI.C 
             ConnectivityCore.C 
             CorotSolver.C 
             CurvatureDetection.C 
             Dev/devtools.C 
             DistGeoState.C
             DistLeastSquareSolver.C
             DistMacroCellCore.C 
             DomainCore.C 
	     Dunavant.C
             EdgeCore.C 
	     EmbeddedCorotSolver.C
	     ExactSolution.C
             FaceCore.C 
             FaceTriaCore.C 
             ElemCore.C 
             ElemTetCore.C 
             FemEquationTermDesc.C 
             FluxFcnDesc.C 
             FluxFcnDescSG.C 
             FluxFcnDescTait.C 
             FluxFcnDescJwl.C
             GeoData.C 
             GeoSource.C 
             GeoState.C 
             GappyPreprocessingCore.C
             HeatTransferHandlerCore.C 
             ImplicitRiemann.C
             InletNodeCore.C 
             IoDataCore.C 
             KspConvCriterion.C 
             MacroCellCore.C 
             MatchNodeCore.C 
             MemoryPool.C 
             MeshMotionHandlerCore.C 
             MeshMotionSolver.C 
             MultiGridPrec.C
             MultiGridDistSVec.C
             MultiGridLevel.C
             MultiGridKernel.C
             MultiGridKspSolver.C
             MultiGridMvpMatrix.C
             MultiGridOperator.C 
             MultiGridSpaceOperator.C 
             MultiGridSmoothingMatrix.C
             MultiGridSmoothingMatrices.C
             MultiGridCoupledTsDesc.C
             MultiGridSegTsDesc.C
             MultiGridEmbeddedTsDesc.C
             ModalSolver.C
             NonlinearRomOfflineSolver.C 
             Node.C 
             ParallelRomCore.C
             PostFcn.C 
             ProgrammedBurnCore.C
             RefVal.C 
             ReinitializeDistanceToWall.C
             SmagorinskyLESTerm.C 
             WaleLESTerm.C 
             DynamicLESTerm.C 
             StructExc.C 
             SubDomainCore.C 
             SparseGridCore.C
             SparseGridClusterCore.C
             TimeData.C 
             Timer.C 
             TsInput.C 
             TsParameters.C 
             TsRestartCore.C 
             VMSLESTerm.C 
             DynamicVMSTerm.C 
             WallFcnCore.C 
             BCApplierCore.C 
	     BCond.C 
             BlockAlloc.C 
             TriangulatedSurface.C
             FluidSelectorCore.C
             OneDimensionalSolver.C
	     Multigrid.C
             PolygonReconstructionData.C
             KirchhoffIntegrator.C
             TriangulatedInterface.C
	     LevelSet/MultiGridLevelSetStructure.C
            # RadialBasisFunction/r8lib.C
            # RadialBasisFunction/RBFInterpND.C
            )

ADD_LIBRARY(Utils utils/Predicate.C utils/Laguerre.C)

SET_TARGET_PROPERTIES(Utils
                PROPERTIES COMPILE_FLAGS -fPIC)

IF(MPI_FOUND)
  ADD_EXECUTABLE(loader tools/loader.C)
  TARGET_LINK_LIBRARIES(loader ${MPI_LIBRARIES} ${CMAKE_DL_LIBS} -lpthread)
  SET_TARGET_PROPERTIES(loader
                      PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BINDIR})
ENDIF(MPI_FOUND)

# Parser
# search bison
MACRO(FIND_BISON)
    IF(NOT BISON_EXECUTABLE)
        FIND_PROGRAM(BISON_EXECUTABLE bison)
        IF (NOT BISON_EXECUTABLE)
          MESSAGE(FATAL_ERROR "bison not found - aborting")
        ENDIF (NOT BISON_EXECUTABLE)
    ENDIF(NOT BISON_EXECUTABLE)
ENDMACRO(FIND_BISON)

MACRO(ADD_BISON_FILES _sources )
    FIND_BISON()

    FOREACH (_current_FILE ${ARGN})
      GET_FILENAME_COMPONENT(_in ${_current_FILE} ABSOLUTE)
      GET_FILENAME_COMPONENT(_basename ${_current_FILE} NAME_WE)

      SET(_out ${CMAKE_CURRENT_BINARY_DIR}/${_basename}.tab.C)

      ADD_CUSTOM_COMMAND(
         OUTPUT ${_out}
         COMMAND ${BISON_EXECUTABLE}
         ARGS
         -d
         -v
         -p yyCmdf
         -o${_out}
         ${_in}
         DEPENDS ${_in}
      )

      SET(${_sources} ${${_sources}} ${_out} )
   ENDFOREACH (_current_FILE)
ENDMACRO(ADD_BISON_FILES)


# flex a .ll file

# search flex
MACRO(FIND_FLEX)
    IF(NOT FLEX_EXECUTABLE)
        FIND_PROGRAM(FLEX_EXECUTABLE flex)
        IF (NOT FLEX_EXECUTABLE)
          MESSAGE(FATAL_ERROR "flex not found - aborting")
        ENDIF (NOT FLEX_EXECUTABLE)
    ENDIF(NOT FLEX_EXECUTABLE)
ENDMACRO(FIND_FLEX)

MACRO(ADD_FLEX_FILES _sources )
    FIND_FLEX()

    FOREACH (_current_FILE ${ARGN})
      GET_FILENAME_COMPONENT(_in ${_current_FILE} ABSOLUTE)
      GET_FILENAME_COMPONENT(_basename ${_current_FILE} NAME_WE)

      SET(_out ${CMAKE_CURRENT_BINARY_DIR}/flex_${_basename}.cpp)

      ADD_CUSTOM_COMMAND(
         OUTPUT ${_out}
         COMMAND ${FLEX_EXECUTABLE}
         ARGS
         -PyyCmdf 
         -o${_out}
         ${_in}
         DEPENDS ${_in}
      )

      SET(${_sources} ${${_sources}} ${_out} )
   ENDFOREACH (_current_FILE)
ENDMACRO(ADD_FLEX_FILES)

ADD_FLEX_FILES(MyLib parser/Lexer.l)
ADD_BISON_FILES(MyLib parser/Parser.y)

ADD_LIBRARY(Parser parser/Dictionary.C parser/StringManager.C parser/AssignerCore.C parser/ParseTree.cpp ${MyLib})
SET_TARGET_PROPERTIES(Parser
                PROPERTIES COMPILE_FLAGS -fPIC)

install(TARGETS FluidBase DESTINATION lib EXPORT fluid-targets)
install(EXPORT fluid-targets DESTINATION lib/base)
export(TARGETS FluidBase FILE ${CMAKE_CURRENT_BINARY_DIR}/base-exports.cmake)

	
IF(MERCURIAL_EXECUTABLE)
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/hgidentity.C 
                   COMMAND echo "const char *THE_VERSION = \\\"" 
                   `${MERCURIAL_EXECUTABLE} identify -in` "\\\"" "\;" > ${CMAKE_CURRENT_BINARY_DIR}/hgidentity.C
                 DEPENDS Main.C FluidBase
                 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  add_definitions(-DPRINT_CHANGESETID)
ELSE (MERCURIAL_EXECUTABLE)
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/hgidentity.C 
                   COMMAND echo "const char *THE_VERSION = \\\"\\\"\;" > ${CMAKE_CURRENT_BINARY_DIR}/hgidentity.C
                 DEPENDS Main.C FluidBase
                 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF (MERCURIAL_EXECUTABLE)

###################
# Local Variables #
###################
#SITE_NAME(HOSTNAME)
#FIND_PACKAGE(LOCAL NAMES ${HOSTNAME} PATHS . ./localconfig QUIET)
#MARK_AS_ADVANCED(HOSTNAME)
#MARK_AS_ADVANCED(LOCAL_DIR)
###################

IF(EXTRALIB)
  SET(EXTRALIBRARIES ${EXTRALIBRARIES} ${EXTRALIB})
ENDIF(EXTRALIB)

SET(INC_DIR .)

FIND_PACKAGE(Boost 1.19.0)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif(Boost_FOUND)

##############################
# For the embedded framework #
##############################
if (NOT WITH_PHYSBAM EQUAL 0)

  if(NOT Boost_FOUND) 
    message(FATAL_ERROR "Compiling with PhysBAM requires boost")
  endif(NOT Boost_FOUND)

  SET(INC_DIR PhysBAM_Interface PhysBAM_Interface/Public_Library)
  INCLUDE_DIRECTORIES(${INC_DIR})

  ADD_LIBRARY(physbam
  PhysBAM_Interface/PHYSBAM_INTERFACE.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Arrays/ARRAY_COLLECTION.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Vectors/VECTOR.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Utilities/EXCEPTIONS.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Utilities/TIMER.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Utilities/PROCESS_UTILITIES.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Data_Structures/KD_TREE.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Parsing/STRING_UTILITIES.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Log/DEBUG_UTILITIES.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Log/LOG.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Log/LOG_ENTRY.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Log/DEBUG_PRINT.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Nonlinear_Equations/ITERATIVE_SOLVER.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Matrices/MATRIX_BASE.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Matrices/MATRIX_3X2.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Matrices/MATRIX_3X3.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Matrices/MATRIX_MXN.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Matrices/ROTATION.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Matrices/MATRIX_4X4.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Polynomials/CUBIC.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Random_Numbers/MT19937.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Basic_Geometry/POLYGON.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Basic_Geometry/SEGMENT_3D.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Basic_Geometry/RAY.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Basic_Geometry/PLANE.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Basic_Geometry/BOX.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Continuous_Collision_Detection/POINT_FACE_COLLISION.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Continuous_Collision_Detection/EDGE_EDGE_COLLISION.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Intersections/BOX_LINE_2D_INTERSECTION.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Intersections/BOX_PLANE_INTERSECTION.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Intersections/BOX_POINT_SIMPLEX_1D_INTERSECTION.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Intersections/RAY_PLANE_INTERSECTION.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Intersections/RAY_BOX_INTERSECTION.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Intersections/RAY_ORIENTED_BOX_INTERSECTION.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Intersections/RAY_SEGMENT_2D_INTERSECTION.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Intersections/RAY_TRIANGLE_3D_INTERSECTION.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Registry/STRUCTURE_REGISTRY.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Topology/SIMPLEX_MESH.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Topology/SEGMENT_MESH.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Topology/TRIANGLE_MESH.cpp
  )
  
  SET_TARGET_PROPERTIES(physbam
                PROPERTIES COMPILE_FLAGS "-fPIC -DCOMPILE_WITHOUT_READ_WRITE_SUPPORT")

  #SET(INC_DIR_TMP . IntersectorFRG IntersectorPhysBAM LevelSet ${PHYSBAM_BASE} ${PHYSBAM_BASE}/Public_Library)
  SET(INC_DIR_TMP . IntersectorFRG IntersectorPhysBAM LevelSet)
  INCLUDE_DIRECTORIES(${INC_DIR_TMP})

  #same definitions as AERO-F main code
  add_definitions(-fPIC -DTEMPLATE_FIX -D__LINUX -DF_NEEDS_UNDSC -DTYPE_PREC=float)
  add_library(Intersector
              IntersectorFRG/IntersectorFRG.C
              IntersectorPhysBAM/IntersectorPhysBAM.C
              IntersectorPhysBAM/FloodFill.C
              IntersectorPhysBAM/Mpi_Utilities.C)
  TARGET_LINK_LIBRARIES(Intersector physbam)

  SET (DO_EMBEDDED 1)
  ADD_DEFINITIONS(-DDO_EMBEDDED)
  SET (EXTRALIBRARIES ${EXTRALIBRARIES} Intersector)
ELSE (NOT WITH_PHYSBAM EQUAL 0)
  SET (DO_EMBEDDED 0)
ENDIF (NOT WITH_PHYSBAM EQUAL 0)
##############################

##############################
# add ALS
##############################

include_directories( AlternatingLeastSquare AlternatingLeastSquare/include)

add_library(als AlternatingLeastSquare/als_lapack.cpp AlternatingLeastSquare/als_util.cpp)

target_link_libraries(als ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

#find_package(MPI REQUIRED)
#include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(als ${MPI_LIBRARIES} -lpthread)
set_target_properties(als PROPERTIES LINK_FLAGS "${MPI_COMPILE_FLAGS}")
set_target_properties(als PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")

set(EXTRALIBRARIES ${EXTRALIBRARIES} als)

##############################

ADD_EXECUTABLE(fluid ${MAIN_FILES} ${CMAKE_CURRENT_BINARY_DIR}/hgidentity.C)

SET_TARGET_PROPERTIES(FluidBase
                PROPERTIES COMPILE_FLAGS -fPIC)

ADD_LIBRARY(FortranLib 
f77src/ballvertex.f
f77src/boundflux5.f
f77src/boundjac2.f
f77src/dsvdc.f
f77src/eriemann_gg.f
f77src/eriemann_gw.f
f77src/eriemann_ww.f
f77src/g_boundflux5.f
f77src/g_boundflux5temp.f
f77src/computedboundfluxoperator.f
f77src/genbcfluxgas.f
f77src/genbcfluxgas_HH.f
f77src/genbcfluxtait.f
f77src/genbcfluxtait_HH.f
f77src/genbcfluxjwl.f
f77src/g_roeflux5.f
f77src/g_roeflux6.f
f77src/g_roeflux6temp.f
f77src/computedfluxoperators.f
f77src/roeflux5jwl.f
f77src/hllcflux.f
f77src/hllcflux1.f
f77src/hllcjac.f
f77src/hlleflux.f
f77src/hlleflux1.f
f77src/hllejac.f
f77src/printmask.f
f77src/roeflux5.f
f77src/roeflux5prim.f
f77src/roeflux5waterdissprec.f
f77src/roeflux5waterdissprec_burnable.f
f77src/roeflux6.f
f77src/roejac2.f
f77src/roejac2water.f
f77src/roejac5.f
f77src/roejac5waterdissprec.f
f77src/roejac6.f
f77src/roejac6jwl.f
f77src/roejacwaterdissprim.f
f77src/roeflux7.f
f77src/torsionspring.f
f77src/dgecp.f
f77src/dgers.f
f77src/diff.f
f77src/g1.f
f77src/h12.f
f77src/spnnls.f)

SET_TARGET_PROPERTIES(FortranLib
                      PROPERTIES
                      COMPILE_FLAGS -fPIC
                      RELEASE_OUTPUT_NAME Fortran
                      DEBUG_OUTPUT_NAME Fortran)

IF (USE_SCALAPACK)
  ADD_LIBRARY(FortranScalapackLib
              f77src/thinsvd.f
              f77src/lworksizesvd.f
              f77src/globalmatrices.f)
  SET_TARGET_PROPERTIES(FortranScalapackLib
                        PROPERTIES
                        COMPILE_FLAGS -fPIC
                        RELEASE_OUTPUT_NAME FortranScalapack
                        DEBUG_OUTPUT_NAME FortranScalapack)
ENDIF(USE_SCALAPACK)


IF(MPI_FOUND)
  MESSAGE("-- Will compile with MPI API " ${MPI_INCLUDE_PATH} ".")
  SET(INC_DIR ${INC_DIR} ${MPI_INCLUDE_PATH})
  add_definitions(-DUSE_MPI -DMPICH_IGNORE_CXX_SEEK)
  
  IF(NOT USE_SCALAPACK)
    ADD_LIBRARY(fluidlib SHARED Main.C ${CMAKE_CURRENT_BINARY_DIR}/hgidentity.C)
    SET_TARGET_PROPERTIES(fluidlib 
                          PROPERTIES COMPILE_FLAGS -DCREATE_DSO)
    TARGET_LINK_LIBRARIES(fluidlib FluidBase Parser FortranLib Utils 
                          ${EXTRALIBRARIES} ${MPI_LIBRARIES} ${FORTRAN_HYBRID_LIBRARIES})
    SET_TARGET_PROPERTIES(fluidlib 
                          PROPERTIES 
                          LIBRARY_OUTPUT_DIRECTORY ${BINDIR}
                          RELEASE_OUTPUT_NAME aerof.opt
                          DEBUG_OUTPUT_NAME aerof.debug
                          PREFIX ""
                         )
    SET(TARGETS_INST ${TARGETS_INST} fluidlib)
  ENDIF(NOT USE_SCALAPACK)

ELSE(MPI_FOUND)
  MESSAGE("-- Only a non-MPI based executable will be compiled.")
  SET(MPI_LIBRARIES "")
ENDIF(MPI_FOUND)

IF(NOT WITH_EIGEN3 EQUAL 0)
  find_package(Eigen3 3.1)
  IF(EIGEN3_FOUND)
    MESSAGE("-- A library with EIGEN3 API found.")
    include_directories(${EIGEN3_INCLUDE_DIR})
    add_definitions(-DUSE_EIGEN3)
    set(CMAKE_MODULE_PATH ${EIGEN3_INCLUDE_DIR}/cmake ${CMAKE_MODULE_PATH})
  ENDIF(EIGEN3_FOUND)
ENDIF(NOT WITH_EIGEN3 EQUAL 0)

# SET(CMAKE_VERBOSE_MAKEFILE true)

#MESSAGE("All include dirs" ${INC_DIR})
INCLUDE_DIRECTORIES( ${INC_DIR} )

#add_subdirectory(parser)
#add_subdirectory(f77src)
#add_subdirectory(utils)
#add_subdirectory(tools)


add_definitions(-DTEMPLATE_FIX -DF_NEEDS_UNDSC -DTYPE_PREC=float)

message( "Extra Libraries are " ${EXTRALIBRARIES} ". MPI libraries are " ${MPI_LIBRARIES} ". CMake DL libraries are " ${CMAKE_DL_LIBS})
IF(NOT USE_SCALAPACK)
  TARGET_LINK_LIBRARIES(fluid FluidBase Parser FortranLib Utils ${EXTRALIBRARIES} ${MPI_LIBRARIES} ${CMAKE_DL_LIBS} -lpthread)
ELSE(NOT USE_SCALAPACK)
  TARGET_LINK_LIBRARIES(fluid FluidBase Parser FortranLib FortranScalapackLib Utils ${EXTRALIBRARIES} ${MPI_LIBRARIES} ${CMAKE_DL_LIBS} -lpthread)
ENDIF(NOT USE_SCALAPACK)


SET_TARGET_PROPERTIES(fluid PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BINDIR}
	RELEASE_OUTPUT_NAME aerof.opt DEBUG_OUTPUT_NAME aerof.debug)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF(NOT CMAKE_BUILD_TYPE)

INSTALL(TARGETS ${TARGETS_INST}
        RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib)

MESSAGE("")
MESSAGE("=================================================")
MESSAGE("           Summary of build options")
MESSAGE("-------------------------------------------------")
IF(MPI_FOUND)
 MESSAGE("Distributed Execution:    YES")
 MESSAGE("Aeroelastic:              YES")
ELSE(MPI_FOUND)
 MESSAGE("Distributed Execution:    NO *******")
 MESSAGE("Aeroelastic:              NO *******")
ENDIF(MPI_FOUND)
IF(DO_EMBEDDED)
 MESSAGE("Embedded framework:       YES")
ELSE(DO_EMBEDDED)
 MESSAGE("Embedded framework:       NO *******")
ENDIF(DO_EMBEDDED)
IF(DO_MODAL)
 MESSAGE("Modal capability:         YES")
ELSE(DO_MODAL)
 MESSAGE("Modal capability:         NO *******")
ENDIF(DO_MODAL)
IF(USE_SCALAPACK)
 MESSAGE("Parallel SVD capability:  YES")
ELSE(USE_SCALAPACK)
 MESSAGE("Parallel SVD capability:  NO *******")
ENDIF(USE_SCALAPACK)
MESSAGE("Aeroacoustic capability:  ${AEROACOUSTIC}")
MESSAGE("Build type:               " ${CMAKE_BUILD_TYPE})
MESSAGE("Extra link flags:         " ${EXTRALIB})
MESSAGE("=================================================")
MESSAGE("")

if(REG_SRC_DIR)
  enable_testing()
  add_subdirectory(${REG_SRC_DIR} ${REG_SRC_DIR})
else(REG_SRC_DIR)
  MESSAGE("To active regression testing, set variable REG_SRC_DIR to location of test suite.")
  MESSAGE("")
endif(REG_SRC_DIR)
