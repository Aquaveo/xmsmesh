//------------------------------------------------------------------------------
/// \file
/// \ingroup meshing_detail
/// \copyright (C) Copyright Aquaveo 2018. Distributed under the xmsng
///  Software License, Version 1.0. (See accompanying file
///  LICENSE_1_0.txt or copy at http://www.aquaveo.com/xmsng/LICENSE_1_0.txt)
//------------------------------------------------------------------------------
#pragma once

//----- Included files ---------------------------------------------------------
//#include <xmscore/points/ptsfwd.h>
#include <xmscore/misc/base_macros.h>
#include <xmscore/misc/boost_defines.h>
#include <xmscore/stl/vector.h>

//----- Forward declarations ---------------------------------------------------

//----- Namespace declaration --------------------------------------------------
namespace xms
{
class Observer;

//----- Constants / Enumerations -----------------------------------------------

//----- Structs / Classes ------------------------------------------------------
/// \brief Generates a mesh with triangle/quad cells using the patch method
/// \see MePolyPatcherImpl
class MePolyPatcher
{
public:
  static BSHP<MePolyPatcher> New();

  /// \cond
  virtual bool MeshIt(const VecPt3d& a_outPoly,
                      const VecInt& a_polyCorners,
                      double a_xytol,
                      VecPt3d& a_points,
                      VecInt& a_cells) = 0;
  virtual void SetObserver(BSHP<Observer> a_) = 0;

protected:
  MePolyPatcher() {}
  virtual ~MePolyPatcher() {}

private:
  XM_DISALLOW_COPY_AND_ASSIGN(MePolyPatcher);
  /// \endcond
};
//----- Function prototypes ----------------------------------------------------

} // namespace xms
