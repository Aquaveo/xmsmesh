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
/// convenience class for holding output data from the MePolyOffsetter
class MePolyOffsetterOutput
{
public:
  std::vector<Pt3d> m_pts;                  ///< locations used by polygons
  std::vector<std::vector<size_t>> m_loops; ///< indexes of points that define loops
  std::vector<int> m_loopTypes;             ///< type of loop
};

/// \brief Does an internal offset from a polygon outer boundary (shrink) and
/// does an external offset from a polygon inner boundary (buffer)
/// \see MePolyOffsetterImpl
class MePolyOffsetter
{
public:
  static BSHP<MePolyOffsetter> New();

  /// enum to identify types of polygons created by this class
  enum polytype { OUTSIDE_POLY = 0, INSIDE_POLY, NEWOUT_POLY };

  /// \cond
  virtual bool Offset(const std::vector<Pt3d>& a_input,
                      polytype a_ptype,
                      MePolyOffsetterOutput& a_out,
                      double a_xyTol) = 0;

protected:
  MePolyOffsetter() {}
  virtual ~MePolyOffsetter() {}

private:
  XM_DISALLOW_COPY_AND_ASSIGN(MePolyOffsetter);
  /// \endcond
};
//----- Function prototypes ----------------------------------------------------

} // namespace xms
