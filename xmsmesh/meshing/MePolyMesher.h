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
#include <boost/shared_ptr.hpp>

// 5. Shared code headers

#include <xmscore/stl/vector.h>
#include <xmscore/misc/base_macros.h> // for XM_DISALLOW_COPY_AND_ASSIGN

//----- Forward declarations ---------------------------------------------------

//----- Namespace declaration --------------------------------------------------

namespace xms
{
class MeMultiPolyMesherIo;
class Observer;
//----- Constants / Enumerations -----------------------------------------------

//----- Structs / Classes ------------------------------------------------------

//----- Function prototypes ----------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// \brief Fills a polygon with a mesh (points and cells). Honors the polygon
/// boundary.
/// \see MePolyMesherImpl
class MePolyMesher
{
public:
  static boost::shared_ptr<MePolyMesher> New();

  MePolyMesher();
  virtual ~MePolyMesher();

  /// \cond
  virtual bool MeshIt(const MeMultiPolyMesherIo& a_input,
                      size_t a_polyIdx,
                      VecPt3d& a_points,
                      VecInt& a_triangles,
                      VecInt& a_cell) = 0;

  virtual void SetObserver(boost::shared_ptr<Observer> a) = 0;
  virtual void GetProcessedRefinePts(VecPt3d& a_pts) = 0;

private:
  XM_DISALLOW_COPY_AND_ASSIGN(MePolyMesher);
  /// \endcond
}; // MePolyMesher

} // namespace xms
