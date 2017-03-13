cmake_minimum_required(VERSION 2.6)

SET(INC_DIR . Public_Library)
INCLUDE_DIRECTORIES(${INC_DIR})

FIND_PACKAGE(Boost 1.19.0 REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif(Boost_FOUND)

ADD_LIBRARY(physbam SHARED PHYSBAM_INTERFACE.cpp
./Public_Library/PhysBAM_Tools/Arrays/ARRAY_COLLECTION.cpp
./Public_Library/PhysBAM_Tools/Vectors/VECTOR.cpp
./Public_Library/PhysBAM_Tools/Utilities/EXCEPTIONS.cpp
./Public_Library/PhysBAM_Tools/Utilities/TIMER.cpp
./Public_Library/PhysBAM_Tools/Utilities/PROCESS_UTILITIES.cpp
./Public_Library/PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.cpp
./Public_Library/PhysBAM_Tools/Data_Structures/KD_TREE.cpp
./Public_Library/PhysBAM_Tools/Parsing/STRING_UTILITIES.cpp
./Public_Library/PhysBAM_Tools/Log/DEBUG_UTILITIES.cpp
./Public_Library/PhysBAM_Tools/Log/LOG.cpp
./Public_Library/PhysBAM_Tools/Log/LOG_ENTRY.cpp
./Public_Library/PhysBAM_Tools/Log/DEBUG_PRINT.cpp
./Public_Library/PhysBAM_Tools/Nonlinear_Equations/ITERATIVE_SOLVER.cpp
./Public_Library/PhysBAM_Tools/Matrices/MATRIX_BASE.cpp
./Public_Library/PhysBAM_Tools/Matrices/MATRIX_3X2.cpp
./Public_Library/PhysBAM_Tools/Matrices/MATRIX_3X3.cpp
./Public_Library/PhysBAM_Tools/Matrices/MATRIX_MXN.cpp
./Public_Library/PhysBAM_Tools/Matrices/ROTATION.cpp
./Public_Library/PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.cpp
./Public_Library/PhysBAM_Tools/Matrices/MATRIX_4X4.cpp
./Public_Library/PhysBAM_Tools/Polynomials/CUBIC.cpp
./Public_Library/PhysBAM_Tools/Random_Numbers/MT19937.cpp
./Public_Library/PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.cpp
./Public_Library/PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.cpp
./Public_Library/PhysBAM_Geometry/Basic_Geometry/POLYGON.cpp
./Public_Library/PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.cpp
./Public_Library/PhysBAM_Geometry/Basic_Geometry/SEGMENT_3D.cpp
./Public_Library/PhysBAM_Geometry/Basic_Geometry/RAY.cpp
./Public_Library/PhysBAM_Geometry/Basic_Geometry/PLANE.cpp
./Public_Library/PhysBAM_Geometry/Basic_Geometry/BOX.cpp
./Public_Library/PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.cpp
./Public_Library/PhysBAM_Geometry/Continuous_Collision_Detection/POINT_FACE_COLLISION.cpp
./Public_Library/PhysBAM_Geometry/Continuous_Collision_Detection/EDGE_EDGE_COLLISION.cpp
./Public_Library/PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.cpp
./Public_Library/PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY.cpp
./Public_Library/PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.cpp
./Public_Library/PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.cpp
./Public_Library/PhysBAM_Geometry/Intersections/BOX_LINE_2D_INTERSECTION.cpp
./Public_Library/PhysBAM_Geometry/Intersections/BOX_PLANE_INTERSECTION.cpp
./Public_Library/PhysBAM_Geometry/Intersections/BOX_POINT_SIMPLEX_1D_INTERSECTION.cpp
./Public_Library/PhysBAM_Geometry/Intersections/RAY_PLANE_INTERSECTION.cpp
./Public_Library/PhysBAM_Geometry/Intersections/RAY_BOX_INTERSECTION.cpp
./Public_Library/PhysBAM_Geometry/Intersections/RAY_ORIENTED_BOX_INTERSECTION.cpp
./Public_Library/PhysBAM_Geometry/Intersections/RAY_SEGMENT_2D_INTERSECTION.cpp
./Public_Library/PhysBAM_Geometry/Intersections/RAY_TRIANGLE_3D_INTERSECTION.cpp
./Public_Library/PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.cpp
./Public_Library/PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE.cpp
./Public_Library/PhysBAM_Geometry/Registry/STRUCTURE_REGISTRY.cpp
./Public_Library/PhysBAM_Geometry/Topology/SIMPLEX_MESH.cpp
./Public_Library/PhysBAM_Geometry/Topology/SEGMENT_MESH.cpp
./Public_Library/PhysBAM_Geometry/Topology/TRIANGLE_MESH.cpp
)

SET_TARGET_PROPERTIES(physbam
                PROPERTIES COMPILE_FLAGS "-fPIC -DCOMPILE_WITHOUT_READ_WRITE_SUPPORT -IPublic_Library")

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF(NOT CMAKE_BUILD_TYPE)
