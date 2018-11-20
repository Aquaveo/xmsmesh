#pragma once
//------------------------------------------------------------------------------
/// \file
/// \ingroup meshing_detail
/// \copyright (C) Copyright Aquaveo 2018. Distributed under FreeBSD License
/// (See accompanying file LICENSE or https://aqaveo.com/bsd/license.txt)
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

class MeBadQuadRemover
{
public:
  static BSHP<MeBadQuadRemover> New(BSHP<XmUGrid> a_ugrid);
  MeBadQuadRemover();
  virtual ~MeBadQuadRemover();

  /// \cond
  virtual BSHP<XmUGrid> RemoveBadQuads(double a_maxAspect = 0.7) = 0;

private:
  XM_DISALLOW_COPY_AND_ASSIGN(MeBadQuadRemover);
  /// \endcond
}; // class MeBadQuadRemover

//----- Function prototypes ----------------------------------------------------

} // namespace xms
