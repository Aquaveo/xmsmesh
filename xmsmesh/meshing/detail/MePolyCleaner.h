//------------------------------------------------------------------------------
/// \file
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

//----- Forward declarations ---------------------------------------------------

//----- Namespace declaration --------------------------------------------------
namespace xms
{
//----- Constants / Enumerations -----------------------------------------------

//----- Structs / Classes ------------------------------------------------------
class MePolyOffsetterOutput;

/// \brief cleans the output produced by MePolyOffsetter
/// \see MePolyCleanerImpl
class MePolyCleaner
{
public:
  static BSHP<MePolyCleaner> New();

  /// \cond
  virtual void CleanPolyOffset(const std::vector<Pt3d>& a_input,
                               int a_pType,
                               double a_tol,
                               MePolyOffsetterOutput& a_out) = 0;

  virtual void IntersectCleanInPolys(const std::vector<MePolyOffsetterOutput>& a_offsets,
                                     MePolyOffsetterOutput& a_out,
                                     double a_xyTol) = 0;

  virtual void IntersectCleanInOutPolys(const MePolyOffsetterOutput& a_offsets,
                                        MePolyOffsetterOutput& a_out,
                                        double a_xyTol) = 0;

protected:
  MePolyCleaner() {}
  virtual ~MePolyCleaner() {}

private:
  XM_DISALLOW_COPY_AND_ASSIGN(MePolyCleaner);
  /// \endcond
};
//----- Function prototypes ----------------------------------------------------

} // namespace xms
