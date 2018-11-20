//------------------------------------------------------------------------------
/// \file
/// \ingroup meshing_detail
/// \copyright (C) Copyright Aquaveo 2018. Distributed under FreeBSD License
/// (See accompanying file LICENSE or https://aqaveo.com/bsd/license.txt)
//------------------------------------------------------------------------------
#pragma once

//----- Included files ---------------------------------------------------------
#include <xmscore/misc/base_macros.h>
#include <xmscore/stl/vector.h>
#include <xmsmesh/meshing/detail/MePolyPts.h>

//----- Forward declarations ---------------------------------------------------

//----- Namespace declaration --------------------------------------------------
class MeIntersectPolysUnitTests;
namespace xms
{
class MePolyOffsetterOutput;
//----- Constants / Enumerations -----------------------------------------------

//----- Structs / Classes ------------------------------------------------------
/// \brief Intersect polygons that are a result of the paving process
/// \see MeIntersectPolysImpl
class MeIntersectPolys
{
  friend MeIntersectPolysUnitTests; ///< tests MeIntersectPolys
public:
  MeIntersectPolys();
  ~MeIntersectPolys();

  void SetupInIn(const std::vector<MePolyOffsetterOutput>& a_offsets, double a_xyTol);
  void SetupInOut(const MePolyOffsetterOutput& a_offsets, double a_xyTol);
  void CalcEnvelopes();
  void InInTrivialPolyCases();
  void InInDoIntersection();
  void InOutDoIntersection();
  void FillOutput(MePolyOffsetterOutput& a_out);
  void ClassifyPolys(const MePolyOffsetterOutput& a_input,
                     std::vector<std::vector<size_t>>& a_output);
  void DeleteBad_NEWOUT_POLY(MePolyOffsetterOutput& a_out, const VecPt3d& a_origOutsidePoly);

private:
  XM_DISALLOW_COPY_AND_ASSIGN(MeIntersectPolys); ///< prevent compiler generated copy/assign
  class impl;
  impl* m_p; ///< implementation class
};

//----- Function prototypes ----------------------------------------------------

} // namespace xms
