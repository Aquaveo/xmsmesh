#pragma once
//------------------------------------------------------------------------------
/// \file
/// \ingroup meshing
/// \copyright (C) Copyright Aquaveo 2018. Distributed under the xmsng
///  Software License, Version 1.0. (See accompanying file
///  LICENSE_1_0.txt or copy at http://www.aquaveo.com/xmsng/LICENSE_1_0.txt)
//------------------------------------------------------------------------------

//----- Included files ---------------------------------------------------------

// 3. Standard library headers
#include <vector>

// 4. External library headers

// 5. Shared code headers
#include <xmscore/misc/base_macros.h>
#include <xmscore/misc/boost_defines.h>
#include <xmscore/stl/vector.h>

//----- Forward declarations ---------------------------------------------------

//----- Namespace declaration --------------------------------------------------

namespace xms
{
//----- Constants / Enumerations -----------------------------------------------

//----- Forward declarations ---------------------------------------------------
class XmUGrid;

//----- Structs / Classes ------------------------------------------------------

class MeQuadBlossom
{
public:
  static BSHP<MeQuadBlossom> New(BSHP<XmUGrid> a_ugrid);
  MeQuadBlossom();
  virtual ~MeQuadBlossom();

  /// \cond
  virtual int PreMakeQuads() = 0;
  virtual BSHP<XmUGrid> MakeQuads(bool a_splitBoundaryPoints,
                                  bool a_useAngle) = 0;
  /// \endcond
  
  static double EstimatedRunTimeInMinutes(int a_numPoints);
  static BSHP<XmUGrid> SplitToQuads(BSHP<XmUGrid> a_ugrid);
  
private:
  /// \cond
  XM_DISALLOW_COPY_AND_ASSIGN(MeQuadBlossom);
  /// \endcond
}; // class MeQuadBlossom

//----- Function prototypes ----------------------------------------------------

} // namespace xms
