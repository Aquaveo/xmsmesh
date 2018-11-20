//------------------------------------------------------------------------------
/// \file
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

//----- Forward declarations ---------------------------------------------------

//----- Namespace declaration --------------------------------------------------
namespace xms
{
class MePolyRedistributePts;

//----- Constants / Enumerations -----------------------------------------------

//----- Structs / Classes ------------------------------------------------------
/// \brief Generates mesh node locations by paving a polygon
/// \see MePolyPaverToMeshPtsImpl
class MePolyPaverToMeshPts
{
public:
  static BSHP<MePolyPaverToMeshPts> New();

  /// \cond
  virtual bool PolyToMeshPts(const std::vector<Pt3d>& a_outPoly,
                             const std::vector<std::vector<Pt3d>>& a_inPolys,
                             double a_bias,
                             double a_xyTol,
                             std::vector<Pt3d>& a_meshPts) = 0;

  virtual void SetRedistributor(BSHP<MePolyRedistributePts> a_) = 0;

private:
  XM_DISALLOW_COPY_AND_ASSIGN(MePolyPaverToMeshPts);
  /// \endcond

protected:
  MePolyPaverToMeshPts();
  virtual ~MePolyPaverToMeshPts();
};
//----- Function prototypes ----------------------------------------------------

} // namespace xms
