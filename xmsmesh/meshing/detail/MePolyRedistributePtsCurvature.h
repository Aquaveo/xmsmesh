//------------------------------------------------------------------------------
/// \file
/// \ingroup meshing
/// \copyright (C) Copyright Aquaveo 2018. Distributed under FreeBSD License
/// (See accompanying file LICENSE or https://aqaveo.com/bsd/license.txt)
//------------------------------------------------------------------------------
#pragma once

//----- Included files ---------------------------------------------------------
#include <vector>
#include <xmscore/stl/vector.h>
#include <xmscore/misc/base_macros.h>
#include <xmscore/misc/boost_defines.h>

//----- Forward declarations ---------------------------------------------------

//----- Namespace declaration --------------------------------------------------
namespace xms
{
//----- Constants / Enumerations -----------------------------------------------

//----- Structs / Classes ------------------------------------------------------
/// \brief Redistributes the point locations on a polyline or polygon based on curvature
/// \see MePolyRedistributePtsCurvatureImpl
class MePolyRedistributePtsCurvature
{
public:
  static BSHP<MePolyRedistributePtsCurvature> New();
  virtual ~MePolyRedistributePtsCurvature();

  /// \cond
  virtual VecPt3d Redistribute(const VecPt3d& points,
                               double a_featureSize,
                               double a_meanSpacing,
                               double a_minimumCurvature = 0.001,
                               bool a_smooth = false) = 0;

private:
  XM_DISALLOW_COPY_AND_ASSIGN(MePolyRedistributePtsCurvature);
  /// \endcond

protected:
  MePolyRedistributePtsCurvature();
};
} // namespace xms
