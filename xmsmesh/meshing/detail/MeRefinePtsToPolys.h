//------------------------------------------------------------------------------
/// \file
/// \brief
/// \ingroup meshing_detail
/// \copyright (C) Copyright Aquaveo 2018. Distributed under the xmsng
///  Software License, Version 1.0. (See accompanying file
///  LICENSE_1_0.txt or copy at http://www.aquaveo.com/xmsng/LICENSE_1_0.txt)
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
  virtual void RefPtsAsPolys(const std::vector<Pt3d>& a_outPoly,
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
