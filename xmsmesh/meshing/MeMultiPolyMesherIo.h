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

// 4. External library headers
#include <xmscore/stl/vector.h>
#include <xmscore/misc/boost_defines.h>

// 5. Shared code headers

//----- Forward declarations ---------------------------------------------------

//----- Namespace declaration --------------------------------------------------

namespace xms
{
class InterpBase;
//----- Constants / Enumerations -----------------------------------------------

//----- Structs / Classes ------------------------------------------------------

//----- Function prototypes ----------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// \class MePolyInput
/// \brief Meshing inputs for one polygon.
class MePolyInput
{
public:
  /// \brief Constructor
  /// \param[in] a_outPoly: see member variable
  /// \param[in] a_insidePolys: see member variable
  /// \param[in] a_bias: see member variable
  /// \param[in] a_sizeFunction: see member variable
  /// \param[in] a_polyCorners: see member variable
  /// \param[in] a_elevFunction: see member variable
  MePolyInput(const VecPt3d& a_outPoly = VecPt3d(),
              const VecPt3d2d& a_insidePolys = VecPt3d2d(),
              double a_bias = 1.0,
              const BSHP<InterpBase> a_sizeFunction = nullptr,
              const VecInt& a_polyCorners = VecInt(),
              const BSHP<InterpBase> a_elevFunction = nullptr)
  : m_outPoly(a_outPoly)
  , m_insidePolys(a_insidePolys)
  , m_bias(a_bias)
  , m_sizeFunction(a_sizeFunction)
  , m_constSizeFunction(-1.0)
  , m_constSizeBias(-1.0)
  , m_polyCorners(a_polyCorners)
  , m_elevFunction(a_elevFunction)
  , m_removeInternalFourTrianglePts(false)
  , m_polyId(-1)
  {
  }

  /// Required. Outer polygons. Clockwise. 1st pt != last.
  VecPt3d m_outPoly;

  /// Optional. Inner polygons (holes). Counter clockwise. 1st pt != last.
  VecPt3d2d m_insidePolys;

  /// Optional. Factor for transitioning between areas of high refinement to
  /// less refinement.
  double m_bias;

  /// Optional. Size function for scalar paving.
  BSHP<InterpBase> m_sizeFunction;
  /// Optional. Constant value size function.
  double m_constSizeFunction;
  /// Optional. Transition factor for constant size function.
  double m_constSizeBias;

  /// Optional. Corner nodes for creating meshes using the patch algorithm.
  /// 3 per outer poly (not 4 - outer poly index point [0] is assumed to be
  /// a corner).
  VecInt m_polyCorners;

  /// Optional. Elevation function for interpolating z coordinate of mesh points.
  BSHP<InterpBase> m_elevFunction;

  /// Optional. Outer boundary locations to remove after the paving process.
  /// Used by the ugAutoCorrectCells class
  VecPt3d m_boundPtsToRemove;

  /// Optional. Remove internal points that are only connected to 4 cells.
  /// Used by the ugAutoCorrectCells class
  bool m_removeInternalFourTrianglePts;

  /// Optional. Polygon id. Useful for reporting errors if calling software stores
  /// ids for polygons and would like to report those back to the user.
  int m_polyId;

  /// Optional. Seed points. If the user has some other methodology for creating
  /// point inside the polygon then those points can be specified here. If these
  /// points are specified then the paving is not performed. These points will
  /// not be used if the meshing option is patch.
  VecPt3d m_seedPoints;
}; // MePolyInput

////////////////////////////////////////////////////////////////////////////////
/// \class MeRefinePoint
/// \brief A refine point used in meshing.
class MeRefinePoint
{
public:
  /// \brief Constructor
  /// \param[in] a_pt: see member variable
  /// \param[in] a_size: see member variable
  /// \param[in] a_createMeshPoint: see member variable
  MeRefinePoint(const Pt3d& a_pt, double a_size, bool a_createMeshPoint)
  : m_pt(a_pt)
  , m_size(a_size)
  , m_createMeshPoint(a_createMeshPoint)
  {
  }

  /// Location of refine point or hard points. Hard points
  /// are points that must be included in the final mesh but have no user
  /// specified size associated with them
  Pt3d m_pt;

  /// Element size at the refine point. A negative value indicates a
  /// "hard point."
  double m_size;

  /// Should a mesh node/point be created at the refine point.
  bool m_createMeshPoint;
}; // MeRefinePoint

////////////////////////////////////////////////////////////////////////////////
/// \class MeMultiPolyMesherIo
/// \brief Provides the input to meshing multiple polygons and holds the output.
class MeMultiPolyMesherIo
{
public:
  /// \brief Constructor
  MeMultiPolyMesherIo()
  : m_polys()
  , m_refPts()
  , m_checkTopology(false)
  , m_returnCellPolygons(true)
  , m_cellPolygons()
  {
  }

  /// Required (but some data is optional). Inputs for each polygon.
  std::vector<MePolyInput> m_polys;

  /// Optional. Refine points.
  std::vector<MeRefinePoint> m_refPts;

  /// Optional. If true, checks polygon input topology for errors.
  bool m_checkTopology;

  /// If true, returns the polygon index of each cell.
  bool m_returnCellPolygons;

  // Output:
  VecPt3d m_points;      ///< The points of the resulting mesh.
  VecInt m_cells;        ///< The cells of the resulting mesh, as a stream.
  VecInt m_cellPolygons; ///< Polygon index of each cell.

}; // MeMultiPolyMesherIo

} // namespace xms
