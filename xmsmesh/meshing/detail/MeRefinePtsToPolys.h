//------------------------------------------------------------------------------
/// \file
/// \brief
/// \ingroup meshing_detail
/// \copyright (C) Copyright Aquaveo 2018. Distributed under FreeBSD License
/// (See accompanying file LICENSE or https://aqaveo.com/bsd/license.txt)
//------------------------------------------------------------------------------
#pragma once

//----- Included files ---------------------------------------------------------
#include <vector>
#include <xmscore/points/ptsfwd.h>
#include <xmscore/misc/base_macros.h>
#include <xmscore/misc/boost_defines.h>

//----- Namespace declaration --------------------------------------------------
namespace xms
{
//----- Forward declarations ---------------------------------------------------
class MeRefinePoint;

//----- Constants / Enumerations -----------------------------------------------

//----- Structs / Classes ------------------------------------------------------
/// \brief Creates polygon from refine point information
/// \see MeRefinePtsToPolysImpl
class MeRefinePtsToPolys
{
public:
  static BSHP<MeRefinePtsToPolys> New();

  /// \cond
  virtual void SetRefinePoints(const std::vector<MeRefinePoint>& a_pts, double a_tol) = 0;
  virtual void RefPtsAsPolys(int a_polyId,
                             const std::vector<Pt3d>& a_outPoly,
                             const std::vector<std::vector<Pt3d>>& a_inPolys,
                             std::vector<std::vector<Pt3d>>& a_newInPolys,
                             std::vector<Pt3d>& a_refMeshPts,
                             std::vector<Pt3d>& a_refPtsTooClose) = 0;

private:
  XM_DISALLOW_COPY_AND_ASSIGN(MeRefinePtsToPolys);
  /// \endcond

protected:
  MeRefinePtsToPolys();
  virtual ~MeRefinePtsToPolys();
};
//----- Function prototypes ----------------------------------------------------

} // namespace xms
