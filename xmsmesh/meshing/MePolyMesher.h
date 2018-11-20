#pragma once
//------------------------------------------------------------------------------
/// \file
/// \ingroup meshing
/// \copyright (C) Copyright Aquaveo 2018. Distributed under FreeBSD License
/// (See accompanying file LICENSE or https://aqaveo.com/bsd/license.txt)
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

  virtual void GetProcessedRefinePts(VecPt3d& a_pts) = 0;

private:
  XM_DISALLOW_COPY_AND_ASSIGN(MePolyMesher);
  /// \endcond
}; // MePolyMesher

} // namespace xms
