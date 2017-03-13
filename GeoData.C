#include <GeoData.h>

//------------------------------------------------------------------------------

GeoData::GeoData(IoData &ioData)
{

// Included (MB)
  if (ioData.problem.alltype == ProblemData::_SHAPE_OPTIMIZATION_ ||
      ioData.problem.alltype == ProblemData::_AEROELASTIC_SHAPE_OPTIMIZATION_ ||
      ioData.problem.alltype == ProblemData::_ROM_SHAPE_OPTIMIZATION_ ||
      ioData.problem.alltype == ProblemData::_SENSITIVITY_ANALYSIS_ ) {

    if (ioData.ts.type != TsData::IMPLICIT)
      ioData.ts.type = TsData::IMPLICIT;

    if (ioData.ts.implicit.type != ImplicitData::SPATIAL_ONLY)
      ioData.ts.implicit.type = ImplicitData::BACKWARD_EULER;

    if (ioData.dgcl.normals != DGCLData::IMPLICIT_FIRST_ORDER_GCL)
      ioData.dgcl.normals = DGCLData::IMPLICIT_FIRST_ORDER_GCL;

    if (ioData.dgcl.velocities != DGCLData::IMPLICIT_ZERO)
      typeVelocities = DGCLData::IMPLICIT_ZERO;

  }

  use_n = false;
  use_nm1 = false;
  use_nm2 = false;
  use_save = false;
  typeNormals = DGCLData::IMPLICIT_FIRST_ORDER_GCL;
  typeVelocities = DGCLData::IMPLICIT_BACKWARD_EULER_VEL;

  if (ioData.problem.type[ProblemData::ACCELERATED] ||
      ioData.problem.type[ProblemData::AERO] ||
      ioData.problem.type[ProblemData::FORCED] ||
      ioData.problem.type[ProblemData::ROLL]) {
    //FOR IMPLICIT SCHEMES
    if (ioData.ts.type == TsData::IMPLICIT) {
      use_n = true;
      if (ioData.ts.implicit.type == ImplicitData::THREE_POINT_BDF)
        use_nm1 = true;
      if (ioData.ts.implicit.type == ImplicitData::FOUR_POINT_BDF) {
        use_nm1 = true;
        use_nm2 = true;
      }

      //Choice of Normals
      if (ioData.dgcl.normals == DGCLData::AUTO) {
        if (ioData.ts.implicit.type == ImplicitData::BACKWARD_EULER ||
            ioData.ts.implicit.type == ImplicitData::CRANK_NICOLSON)
          typeNormals = DGCLData::IMPLICIT_FIRST_ORDER_GCL;
        else if (ioData.ts.implicit.type == ImplicitData::THREE_POINT_BDF)
          typeNormals = DGCLData::IMPLICIT_SECOND_ORDER_GCL;
        else if (ioData.ts.implicit.type == ImplicitData::FOUR_POINT_BDF)
          typeNormals = DGCLData::IMPLICIT_THIRD_ORDER_EZGCL;
      }
      else
        typeNormals = ioData.dgcl.normals;

      //Choice of Velocities
      if (ioData.dgcl.velocities == DGCLData::AUTO_VEL)
        typeVelocities = DGCLData::IMPLICIT_BACKWARD_EULER_VEL;
      else
        typeVelocities = ioData.dgcl.velocities;
    }

    //FOR EXPLICIT SCHEMES
    else if (ioData.ts.type == TsData::EXPLICIT) {
      use_n = true;
      if (ioData.ts.expl.type == ExplicitData::RUNGE_KUTTA_2 ||
          ioData.ts.expl.type == ExplicitData::ONE_BLOCK_RK2 ||
          ioData.ts.expl.type == ExplicitData::ONE_BLOCK_RK2bis)
        use_save = true;

      //Choice of Normals
      if (ioData.dgcl.normals == DGCLData::AUTO) {
        if (ioData.ts.expl.type == ExplicitData::RUNGE_KUTTA_2 ||
            ioData.ts.expl.type == ExplicitData::ONE_BLOCK_RK2 ||
            ioData.ts.expl.type == ExplicitData::ONE_BLOCK_RK2bis)
          typeNormals = DGCLData::EXPLICIT_RK2;
        else
          typeNormals = DGCLData::IMPLICIT_LATEST_CFG;
      }
      else
        typeNormals = ioData.dgcl.normals;

      //Choice of Velocities
      if (ioData.dgcl.velocities == DGCLData::AUTO_VEL) {
        if (ioData.ts.expl.type == ExplicitData::RUNGE_KUTTA_2 ||
            ioData.ts.expl.type == ExplicitData::ONE_BLOCK_RK2 ||
            ioData.ts.expl.type == ExplicitData::ONE_BLOCK_RK2bis)
          typeVelocities = DGCLData::EXPLICIT_RK2_VEL;
        else
          typeVelocities = DGCLData::IMPLICIT_BACKWARD_EULER_VEL;
      }
      else
        typeVelocities = ioData.dgcl.velocities;
    }

  }
  if (ioData.problem.type[ProblemData::LINEARIZED])  {
    use_n = true;
    use_nm1 = true;
  }
}

GeoData::GeoData(const GeoData& d) :
    typeNormals(d.typeNormals), typeVelocities(d.typeVelocities), config(d.config),
    configSA(d.configSA), use_n(d.use_n), use_nm1(d.use_nm1), use_nm2(d.use_nm2), use_save(d.use_save)
{}

//------------------------------------------------------------------------------
