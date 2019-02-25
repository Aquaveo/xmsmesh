//------------------------------------------------------------------------------
/// \file
/// \ingroup meshing
/// \copyright (C) Copyright Aquaveo 2018. Distributed under FreeBSD License
/// (See accompanying file LICENSE or https://aqaveo.com/bsd/license.txt)
//------------------------------------------------------------------------------

//----- Included files ---------------------------------------------------------

// 1. Precompiled header

// 2. My own header
#include <xmsmesh/meshing/MeMultiPolyMesher.h>

// 3. Standard library headers
#include <fstream>
#include <numeric>
#include <set>
#include <sstream>

// 4. External library headers
#include <boost/format.hpp>
#include <boost/unordered_map.hpp>
#include <xmscore/misc/StringUtil.h>
#include <xmscore/misc/Progress.h>
#include <xmscore/misc/XmError.h>
#include <xmscore/stl/set.h>
#include <xmscore/stl/vector.h>
#include <xmsinterp/geometry/GmPtSearch.h>
#include <xmsinterp/geometry/geoms.h>
#include <xmsinterp/interpolate/InterpIdw.h>
#include <xmsinterp/interpolate/InterpLinear.h>

// 5. Shared code headers
#include <xmsmesh/meshing/MeMultiPolyMesherIo.h>
#include <xmsmesh/meshing/MePolyMesher.h>

// 6. Non-shared code headers

//----- Forward declarations ---------------------------------------------------

//----- External globals -------------------------------------------------------

//----- Namespace declaration --------------------------------------------------

namespace xms
{
//----- Constants / Enumerations -----------------------------------------------

//----- Classes / Structs ------------------------------------------------------
class MeMultiPolyMesherImpl : public MeMultiPolyMesher
{
public:
  MeMultiPolyMesherImpl();
  ~MeMultiPolyMesherImpl() {}

  virtual bool MeshIt(MeMultiPolyMesherIo& a_io) override;

private:
  void AppendMesh(VecPt3d& a_points, const VecInt& a_triangles, VecInt& a_cells);
  void AddUniquePoints(const VecPt3d& a_points, VecInt& a_oldDups, VecInt& a_newDups);
  void RenumberNewMesh(int a_oldNumPts,
                       size_t a_numNewPts,
                       const VecInt& a_oldDups,
                       const VecInt& a_newDups,
                       VecInt& a_cells) const;
  void AppendNewCells(const VecInt& a_cells);
  void ReportUnusedRefinePts(const MeMultiPolyMesherIo& a_io, const VecPt3d& a_usedPts);
  void EnsureProperPolygonInputs(MeMultiPolyMesherIo& a_io);
  bool ValidateInput(const MeMultiPolyMesherIo& a_io);
  void CheckForIntersections(const MeMultiPolyMesherIo& a_io, std::string& a_errors) const;
  bool ExtentsOverlap(const Pt3d& oneMn,
                      const Pt3d& oneMx,
                      const Pt3d& two1,
                      const Pt3d& two2) const;

  BSHP<VecPt3d> m_pts; ///< Mesh points. BSHP because of PtSearch::VectorThatGrowsToSearch
  VecInt m_cells;      ///< Mesh cells as a stream
  int m_cellCount;     ///< Number of cells
  boost::unordered_map<std::pair<double, double>, int> m_ptHash; ///< hashes points
};                                                               // class MeMultiPolyMesherImpl

////////////////////////////////////////////////////////////////////////////////
/// \brief Struct defining the location of a polygon line segment in the list of
/// polygons in MeMultiPolyMesherIo.
struct SegmentLocation
{
public:
  /// \brief Constructor
  /// \param a_poly: Index of polygon in MeMultiPolyMesherIo::m_polys.
  /// \param a_polySeg: Index of segment in MePolyInput::m_outPoly
  /// \param a_inPoly: Index of inner polygon in MePolyInput::m_insidePolys
  /// \param a_inPolySeg: Index of segment on inner polygon
  SegmentLocation(int a_poly = 0, int a_polySeg = 0, int a_inPoly = 0, int a_inPolySeg = 0)
  : m_poly(a_poly)
  , m_polySeg(a_polySeg)
  , m_inPoly(a_inPoly)
  , m_inPolySeg(a_inPolySeg)
  {
  }

  int m_poly;      ///< Index of polygon in MeMultiPolyMesherIo::m_polys
  int m_polySeg;   ///< Index of segment in MePolyInput::m_outPoly
  int m_inPoly;    ///< Index of inner polygon in MePolyInput::m_insidePolys
  int m_inPolySeg; ///< Index of segment on inner polygon
};                 // struct SegmentLocation

//----- Internal functions -----------------------------------------------------
namespace
{
//------------------------------------------------------------------------------
/// \brief Writes the Interpolation data to a file if a debug file is present
/// \param a_os[in]: the file
/// \param a_interp[in]: the interpolation class
//------------------------------------------------------------------------------
void iWriteInterpDataToDebugFile(std::fstream& a_os, BSHP<InterpBase> a_interp)
{
  BSHP<VecPt3d> ptsPtr = a_interp->GetPts();
  BSHP<InterpIdw> idw = BDPC<InterpIdw>(a_interp);
  BSHP<InterpLinear> linear = BDPC<InterpLinear>(a_interp);
  XM_ENSURE_TRUE(ptsPtr && (idw || linear));
  if (idw)
    a_os << "IDW";
  else
    a_os << "LINEAR";

  VecPt3d& pts(*ptsPtr);
  a_os << pts.size() << "\n";
  for (const auto& p : pts)
    a_os << STRstd(p.x) << " " << STRstd(p.y) << "\n";
} // iWriteInterpDataToDebugFile
//------------------------------------------------------------------------------
/// \brief Writes the MeMultiPolyMesherIo to a file if a debug file is present
/// \param a_io[in] the input to the mesher
//------------------------------------------------------------------------------
void iWriteInputsToDebugFile(MeMultiPolyMesherIo& a_io)
{
  std::ifstream f("c:\\temp\\xmsmesh_write_inputs.dbg");
  if (!f.good())
    return;
  std::fstream os;
  os.open("C:\\temp\\xmsmesh_inputs.txt", std::fstream::out);
  if (os.bad())
    return;
  for (size_t i = 0; i < a_io.m_polys.size(); ++i)
  {
    MePolyInput& poly(a_io.m_polys[i]);
    os << "BEGIN_POLYGON\nOUTSIDE " << poly.m_outPoly.size() << "\n";
    for (auto& p : poly.m_outPoly)
      os << STRstd(p.x) << " " << STRstd(p.y) << "\n";
    for (auto& v : poly.m_insidePolys)
    {
      os << "INSIDE " << v.size() << "\n";
      for (auto& p : v)
        os << STRstd(p.x) << " " << STRstd(p.y) << "\n";
    }
    os << "BIAS " << STRstd(poly.m_bias) << "\n";
    if (poly.m_sizeFunction)
    {
      os << "SIZE_FUNCTION\n";
      iWriteInterpDataToDebugFile(os, poly.m_sizeFunction);
    }
    if (poly.m_elevFunction)
    {
      os << "ELEVATION_FUNCTION\n";
      iWriteInterpDataToDebugFile(os, poly.m_elevFunction);
    }
    if (poly.m_constSizeFunction != -1)
    {
      os << "CONST_SIZE_FUNCTION " << STRstd(poly.m_constSizeFunction) << "\n";
    }
    if (!poly.m_polyCorners.empty())
    {
      os << "PATCH_CORNERS ";
      for (size_t j = 0; j < poly.m_polyCorners.size(); ++j)
        os << poly.m_polyCorners[j] << " ";
      os << "\n";
    }
    if (poly.m_relaxationMethod != "")
    {
      os << "RELAXATION_METHOD " << poly.m_relaxationMethod << "\n";
    }
    os << "END_POLYGON\n";
  }
  if (a_io.m_checkTopology)
  {
    os << "CHECK_TOPOLOGY\n";
  }
  if (a_io.m_returnCellPolygons)
  {
    os << "RETURN_CELL_POLYGONS\n";
  }
  if (!a_io.m_refPts.empty())
  {
    os << "REFINE_POINTS " << a_io.m_refPts.size() << "\n";
    for (auto& p : a_io.m_refPts)
    {
      os << STRstd(p.m_pt.x) << " " << STRstd(p.m_pt.y) << " " << STRstd(p.m_size) << " "
         << p.m_createMeshPoint << "\n";
    }
  }

} // iWriteInputsToDebugFile

} // unnamed namespace
//----- Class / Function definitions -------------------------------------------

//////////////////////////////////////////////////////////////////////////////
/// \class MeMultiPolyMesherImpl
/// \brief Creates a mesh from multiple polygons that will honor polygon
/// boundaries.
//////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief Constructor.
//------------------------------------------------------------------------------
MeMultiPolyMesherImpl::MeMultiPolyMesherImpl()
: m_pts(new VecPt3d())
, m_cells()
, m_cellCount(0)
{
} // MeMultiPolyMesherImpl::MeMultiPolyMesherImpl
//------------------------------------------------------------------------------
/// \brief Creates a triangle mesh from the input polygons. The polygons can
/// not overlap.
/// \param[in] a_io: Input/output of polygons and options for generating a mesh.
///   MeMultiPolyMesherIo::m_returnCellPolygons is true.
/// \return true if successful.
//------------------------------------------------------------------------------
bool MeMultiPolyMesherImpl::MeshIt(MeMultiPolyMesherIo& a_io)
{
  iWriteInputsToDebugFile(a_io);
  EnsureProperPolygonInputs(a_io);
  if (!ValidateInput(a_io))
  {
    return false;
  }

  // Mesh each polygon and merge the triangles together into one mesh
  BSHP<MePolyMesher> pm = MePolyMesher::New();
  std::stringstream ss;
  ss << "Meshing polygon 1 of " << a_io.m_polys.size();
  Progress prog(ss.str());

  VecPt3d pts, refinePts, tmpPts;
  VecInt tris, cells, cellPolygons;
  for (size_t i = 0; i < a_io.m_polys.size(); ++i)
  {
    pts.resize(0);
    tris.resize(0);
    if (pm->MeshIt(a_io, i, pts, tris, cells))
    {
      pm->GetProcessedRefinePts(tmpPts);
      refinePts.insert(refinePts.end(), tmpPts.begin(), tmpPts.end());
      AppendMesh(pts, tris, cells);

      // Assign cell polygons
      if (a_io.m_returnCellPolygons)
      {
        cellPolygons.resize(m_cellCount, (int)i);
      }
    }

    // Update progress
    ss.str("");
    ss << "Meshing polygon " << i + 2 << " of " << a_io.m_polys.size();
    prog.UpdateMessage(ss.str());
    prog.ProgressStatus((double)i / a_io.m_polys.size());
  }

  // Move memory and cleanup
  a_io.m_points.swap(*m_pts);
  a_io.m_cells.swap(m_cells);
  a_io.m_cellPolygons.swap(cellPolygons);
  m_ptHash.clear();
  m_cellCount = 0;

  // report unused refine points
  ReportUnusedRefinePts(a_io, refinePts);
  return true;
} // MeMultiPolyMesherImpl::MeshIt
//------------------------------------------------------------------------------
/// \brief Remove last point of polygon if it is the same as the first point and
/// make sure the polygon points are ordered correctly.
/// \param a_io: The input/output parameters.
//------------------------------------------------------------------------------
void MeMultiPolyMesherImpl::EnsureProperPolygonInputs(MeMultiPolyMesherIo& a_io)
{
  double tol(1e-9);
  // remove repeated first and last points
  for (size_t i = 0; i < a_io.m_polys.size(); ++i)
  {
    MePolyInput& polyInput = a_io.m_polys[i];
    auto &first = polyInput.m_outPoly.front();
    auto &last = polyInput.m_outPoly.back();
    if (gmEqualPointsXY(first, last, tol))
    {
      polyInput.m_outPoly.pop_back();
    }

    for (size_t j = 0; j < polyInput.m_insidePolys.size(); ++j)
    {
      auto &first = polyInput.m_insidePolys[j].front();
      auto &last = polyInput.m_insidePolys[j].back();
      if (gmEqualPointsXY(first, last, tol))
      {
        polyInput.m_insidePolys[j].pop_back();
      }
    }
  }
  // check point ordering. We want the outside polygon should be CW and the
  // inside polygons should be CCW
  for (size_t i = 0; i < a_io.m_polys.size(); ++i)
  {
    MePolyInput& polyInput = a_io.m_polys[i];
    double area = gmPolygonArea(&polyInput.m_outPoly[0], polyInput.m_outPoly.size());
    if (area > 0)
    {
      std::reverse(polyInput.m_outPoly.begin(), polyInput.m_outPoly.end());
      std::cout << "Reverse outer polygon\n";
    }

    for (size_t j = 0; j < polyInput.m_insidePolys.size(); ++j)
    {
      area = gmPolygonArea(&polyInput.m_insidePolys[j][0], polyInput.m_insidePolys[j].size());
      if (area < 0)
      {
        std::reverse(polyInput.m_insidePolys[j].begin(), polyInput.m_insidePolys[j].end());
        std::cout << "Reverse inner polygon\n";
      }
    }
  }
} // MeMultiPolyMesherImpl::EnsureProperPolygonInputs
//------------------------------------------------------------------------------
/// \brief Make sure the input makes sense.
/// \param a_io: The input/output parameters.
/// \return true if OK, false if there are errors in the input.
//------------------------------------------------------------------------------
bool MeMultiPolyMesherImpl::ValidateInput(const MeMultiPolyMesherIo& a_io)
{
  std::string errors;
  if (a_io.m_polys.empty())
  {
    errors += "Error: Per polygon input is empty. No polygons to mesh.\n";
  }
  else
  {
    for (size_t i = 0; i < a_io.m_polys.size(); ++i)
    {
      const MePolyInput& polyInput = a_io.m_polys[i];
      std::string id;
      {
        std::stringstream ss;
        if (polyInput.m_polyId > -1)
          ss << "id " << polyInput.m_polyId;
        else
          ss << "index " << i;
        id = ss.str();
      }

      // Check outer polygon size
      if (polyInput.m_outPoly.empty())
      {
        std::stringstream ss;
        ss << "Error: Outer polygon " << id << " is empty.\n";
        errors += ss.str();
      }

      // Check inner polygons (if any)
      for (size_t j = 0; j < polyInput.m_insidePolys.size(); ++j)
      {
        if (polyInput.m_insidePolys[j].empty())
        {
          std::stringstream ss;
          ss << "Error: Inner polygon " << j << " of outer polygon " << id << " is empty.\n";
          errors += ss.str();
        }
      }

      // Check patch corners size
      if (polyInput.m_polyCorners.size() != 0 && polyInput.m_polyCorners.size() != 3)
      {
        std::stringstream ss;
        ss << "Error: Polygon patch corners for polygon " << id << " is size "
           << polyInput.m_polyCorners.size() << ". It must be size 0 or 3.\n";
        errors += ss.str();
      }
    }

    if (a_io.m_checkTopology)
    {
      CheckForIntersections(a_io, errors);
    }
  }

  if (!errors.empty())
  {
    XM_ASSERT(false);
    XM_LOG(xmlog::warning, errors);
    return false;
  }
  return true;
} // MeMultiPolyMesherImpl::ValidateInput
//------------------------------------------------------------------------------
/// \brief Checks input to see if anything intersects. Will not catch polys
///        entirely inside or outside of where they are supposed to be.
///
/// We check by intersecting every line segment on every poly with every other
/// line segment on every other poly, if the extents overlap. This can be slow.
/// \param a_io: Mesher input/output
/// \param a_errors: Error string that may get appended to.
//------------------------------------------------------------------------------
void MeMultiPolyMesherImpl::CheckForIntersections(const MeMultiPolyMesherIo& a_io,
                                                  std::string& a_errors) const
{
  // Create vectors defining all poly segments and where they're located

  std::vector<const Pt3d*> segments;
  std::vector<SegmentLocation> segmentLocs;
  for (int i = 0; i < (int)a_io.m_polys.size(); ++i)
  {
    if (i > 0)
    {
      // Add a null to separate polygons
      segments.push_back(nullptr);
      segmentLocs.push_back(SegmentLocation());
    }

    const MePolyInput& inputPoly = a_io.m_polys[i];
    const VecPt3d& outerPoly = inputPoly.m_outPoly;

    // Do the outer poly
    for (int j = 0; j < (int)outerPoly.size(); ++j)
    {
      segments.push_back(&outerPoly[j]);
      segmentLocs.push_back(SegmentLocation(i, j, -1, -1));
    }
    // Do the last segment because polygon isn't closed
    segments.push_back(&outerPoly[0]);
    segmentLocs.push_back(SegmentLocation(i, 0, -1, -1));

    // Do the inner polys

    if (!inputPoly.m_insidePolys.empty())
    {
      // Add a null to separate polygons
      segments.push_back(nullptr);
      segmentLocs.push_back(SegmentLocation());
    }

    for (int k = 0; k < (int)inputPoly.m_insidePolys.size(); ++k)
    {
      if (k > 0)
      {
        // Add a null to separate polygons
        segments.push_back(nullptr);
        segmentLocs.push_back(SegmentLocation());
      }
      const VecPt3d& innerPoly = inputPoly.m_insidePolys[k];

      // Do the inner poly
      for (int m = 0; m < (int)innerPoly.size(); ++m)
      {
        segments.push_back(&innerPoly[m]);
        segmentLocs.push_back(SegmentLocation(i, -1, k, m));
      }
      // Do the last segment because polygon isn't closed
      segments.push_back(&innerPoly[0]);
      segmentLocs.push_back(SegmentLocation(i, -1, k, 0));
    }
  }

  // Intersect all segments with all other segments

  double xi, yi, zi1, zi2, tol(gmXyTol());
  size_t size = segments.size();

  // Loop through every segment
  for (size_t i = 1; i < segments.size(); ++i)
  {
    // Skip nulls which indicate the end of a poly
    if (segments[i] == nullptr)
    {
      ++i;
      continue;
    }
    const Pt3d* one1 = segments[i - 1];
    const Pt3d* one2 = segments[i];
    Pt3d oneMn, oneMx;
    gmAddToExtents(*one1, oneMn, oneMx);
    gmAddToExtents(*one2, oneMn, oneMx);
    oneMn -= tol;
    oneMx += tol;

    // Starting on the next segment from where we are, go to end of segments
    for (size_t j = i + 1; j < size; ++j)
    {
      // Skip nulls which indicate the end of a poly
      if (segments[j] == nullptr)
      {
        ++j;
        continue;
      }
      const Pt3d* two1 = segments[j - 1];
      const Pt3d* two2 = segments[j];
      if (ExtentsOverlap(oneMn, oneMx, *two1, *two2))
      {
        if (gmIntersectLineSegmentsWithTol(*one1, *one2, *two1, *two2, &xi, &yi, &zi1, &zi2, tol))
        {
          // See if we didn't just intersect on the ends
          if (!gmEqualPointsXY(one1->x, one1->y, xi, yi, tol) &&
              !gmEqualPointsXY(one2->x, one2->y, xi, yi, tol))
          {
            // Report the intersection
            std::stringstream ss;
            SegmentLocation& a0 = segmentLocs[i - 1];
            SegmentLocation& a1 = segmentLocs[i];
            SegmentLocation& b0 = segmentLocs[j - 1];
            SegmentLocation& b1 = segmentLocs[j];
            std::string s;
            if (a0.m_inPoly == -1)
            { // First segment on an outer poly
              if (b0.m_inPoly == -1)
              { // Second segment on an outer poly
                s = (boost::format("Error: Input polygon segments intersect."
                                   " The segment defined by points %d and %d"
                                   " of outer polygon %d"
                                   " intersects with the segment defined by points %d and %d"
                                   " of outer polygon %d.\n") %
                     a0.m_polySeg % a1.m_polySeg % a0.m_poly % b0.m_polySeg % b1.m_polySeg %
                     b0.m_poly)
                      .str();
              }
              else
              { // Second segment on an inner poly
                s = (boost::format("Error: Input polygon segments intersect."
                                   " The segment defined by points %d and %d"
                                   " of outer polygon %d"
                                   " intersects with the segment defined by points %d and %d"
                                   " of inner polygon %d of outer polygon %d.\n") %
                     a0.m_polySeg % a1.m_polySeg % a0.m_poly % b0.m_inPolySeg % b1.m_inPolySeg %
                     b0.m_inPoly % b0.m_poly)
                      .str();
              }
            }
            else
            { // First segment on an inner poly
              if (b0.m_inPoly == -1)
              { // Second segment on an outer poly
                s = (boost::format("Error: Input polygon segments intersect."
                                   " The segment defined by points %d and %d"
                                   " of inner polygon %d of outer polygon %d"
                                   " intersects with the segment defined by points %d and %d"
                                   " of outer polygon %d.\n") %
                     a0.m_inPolySeg % a1.m_inPolySeg % a0.m_inPoly % a0.m_poly % b0.m_polySeg %
                     b1.m_polySeg % b0.m_poly)
                      .str();
              }
              else
              { // Second segment on an inner poly
                s = (boost::format("Error: Input polygon segments intersect."
                                   " The segment defined by points %d and %d"
                                   " of inner polygon %d of outer polygon %d"
                                   " intersects with the segment defined by points %d and %d"
                                   " of inner polygon %d of outer polygon %d.\n") %
                     a0.m_inPolySeg % a1.m_inPolySeg % a0.m_inPoly % a0.m_poly % b0.m_inPolySeg %
                     b1.m_inPolySeg % b0.m_inPoly % b0.m_poly)
                      .str();
              }
            }
            a_errors += s;
          }
        }
      }
    }
  }
} // MeMultiPolyMesherImpl::CheckForIntersections
//------------------------------------------------------------------------------
/// \brief Given extents of one segment, and another segment, return true if
///        the extents of the two segments overlap.
/// \param a_oneMn: Min extents of first segment.
/// \param a_oneMx: Max extents of first segment.
/// \param a_two1: First point defining second segment.
/// \param a_two2: Second point defining second segment.
/// \return true if overlap, else false.
//------------------------------------------------------------------------------
bool MeMultiPolyMesherImpl::ExtentsOverlap(const Pt3d& a_oneMn,
                                           const Pt3d& a_oneMx,
                                           const Pt3d& a_two1,
                                           const Pt3d& a_two2) const
{
  Pt3d twoMn, twoMx;
  gmAddToExtents(a_two1, twoMn, twoMx);
  gmAddToExtents(a_two2, twoMn, twoMx);
  if (twoMn.x > a_oneMx.x)
  {
    return false;
  }
  if (twoMx.x < a_oneMn.x)
  {
    return false;
  }
  if (twoMn.y > a_oneMx.y)
  {
    return false;
  }
  if (twoMx.y < a_oneMn.y)
  {
    return false;
  }
  return true;
} // MeMultiPolyMesherImpl::ExtentsOverlap
//------------------------------------------------------------------------------
/// \brief Adds new points and triangles to existing mesh, hashing points and
///        renumbering.
/// \param a_points: New mesh points.
/// \param a_triangles: New mesh triangle cells.
/// \param a_cells: New mesh cells.
//------------------------------------------------------------------------------
void MeMultiPolyMesherImpl::AppendMesh(VecPt3d& a_points,
                                       const VecInt& a_triangles,
                                       VecInt& a_cells)
{
  // These correspond to defines in vtkCellType. There are classes downstream from
  // meshing that convert a cell stream into a vtkUnstructured grid that use these
  // magic numbers.
  const int VTK_TRI(5);

  VecInt cells;
  cells.swap(a_cells);
  int cellCount = 0;
  if (!a_triangles.empty())
  {
    cellCount = (int)a_triangles.size() / 3;
    cells.reserve(cellCount * 5);
    for (size_t i = 0; i < a_triangles.size(); i += 3)
    {
      cells.push_back(VTK_TRI);
      cells.push_back(3);
      cells.push_back(a_triangles[i + 0]);
      cells.push_back(a_triangles[i + 1]);
      cells.push_back(a_triangles[i + 2]);
    }
  }
  else
  {
    size_t i = 0;
    while (i < cells.size())
    {
      cellCount++;
      i++; // celltype = m_cells[i++];
      int npts = cells[i++];
      i += npts;
    }
  }
  m_cellCount += cellCount;

  if ((*m_pts).empty())
  {
    std::pair<std::pair<double, double>, int> pd;
    for (size_t i = 0; i < a_points.size(); ++i)
    {
      pd.first.first = a_points[i].x;
      pd.first.second = a_points[i].y;
      pd.second = (int)i;
      m_ptHash.insert(pd);
    }
    // First time, just add points to m_pts and return
    m_pts->swap(a_points);
    m_cells.swap(cells);
    return;
  }

  VecInt oldDups, newDups;
  int oldNumPts = (int)(*m_pts).size();
  AddUniquePoints(a_points, oldDups, newDups);

  RenumberNewMesh(oldNumPts, a_points.size(), oldDups, newDups, cells);
  AppendNewCells(cells);
} // MeMultiPolyMesherImpl::AppendMesh
//------------------------------------------------------------------------------
/// \brief Adds new mesh points that aren't in old mesh, to old mesh.
///
/// Also gives back lists of duplicate points found (if any) in a_oldDups and
/// a_newDups which are the same size and correspond with each other
/// (a_oldDups[3] => a_newDups[3])
/// \param a_points:  New mesh points.
/// \param a_oldDups: Indices in m_pts (old mesh) of duplicate points.
/// \param a_newDups: Indices in a_points (new mesh) of duplicate points.
//------------------------------------------------------------------------------
void MeMultiPolyMesherImpl::AddUniquePoints(const VecPt3d& a_points,
                                            VecInt& a_oldDups,
                                            VecInt& a_newDups)
{
  a_oldDups.clear();
  a_newDups.clear();
  XM_ENSURE_TRUE_VOID_NO_ASSERT(!a_points.empty());

  // const double tol = 1e-9;
  // int idxOld;
  std::pair<std::pair<double, double>, int> pd;
  auto itEnd = m_ptHash.end();
  for (size_t i = 0; i < a_points.size(); ++i)
  {
    pd.first.first = a_points[i].x;
    pd.first.second = a_points[i].y;
    pd.second = (int)m_pts->size();
    auto it = m_ptHash.find(pd.first);
    if (it == itEnd)
    {
      m_ptHash.insert(pd);
      m_pts->push_back(a_points[i]);
    }
    else
    {
      a_oldDups.push_back(it->second);
      a_newDups.push_back((int)i);
    }
  }
} // MeMultiPolyMesherImpl::AddUniquePoints
//------------------------------------------------------------------------------
/// \brief Renumber points referred to in the new mesh to remove duplicate
///        points and to start new mesh numbering at size of existing mesh.
/// \param a_oldNumPts: Number of points we had before adding this mesh.
/// \param a_numNewPts: Number of new, unique points
/// \param a_oldDups: Indices in m_pts of duplicate points.
/// \param a_newDups: Indices in a_points of duplicate points.
/// \param a_cells: Cells of mesh (same format as vtk cell stream)
//------------------------------------------------------------------------------
void MeMultiPolyMesherImpl::RenumberNewMesh(int a_oldNumPts,
                                            size_t a_numNewPts,
                                            const VecInt& a_oldDups,
                                            const VecInt& a_newDups,
                                            VecInt& a_cells) const
{
  XM_ENSURE_TRUE_VOID_NO_ASSERT(a_oldDups.size() == a_newDups.size());

  // Renumber the new points starting at the old number of points
  VecInt newPtsIdx(a_numNewPts, -1);
  if (a_newDups.empty())
  {
    std::iota(newPtsIdx.begin(), newPtsIdx.end(), a_oldNumPts);
  }
  else
  {
    // Set the new points that are duplicates equal to the old point indices
    for (size_t i = 0; i < a_newDups.size(); ++i)
    {
      newPtsIdx[a_newDups[i]] = a_oldDups[i];
    }

    // Renumber, skipping duplicate points we just assigned
    int cnt(a_oldNumPts);
    for (size_t i = 0; i < newPtsIdx.size(); ++i)
    {
      if (newPtsIdx[i] < 0)
      {
        newPtsIdx[i] = cnt;
        ++cnt;
      }
    }
  }

  // Fix the cell numbers to correspond to the new point numbers
  for (size_t i = 1; i < a_cells.size(); ++i)
  {
    int nPts = a_cells[i];
    for (int j = 1; j <= nPts; ++j)
    {
      a_cells[i + j] = newPtsIdx[a_cells[i + j]];
    }
    i += (1 + nPts);
  }
} // MeMultiPolyMesherImpl::RenumberNewMesh
//------------------------------------------------------------------------------
/// \brief Append the new cells to the existing cells
/// \param a_cells: New cells to append to existing mesh.
//------------------------------------------------------------------------------
void MeMultiPolyMesherImpl::AppendNewCells(const VecInt& a_cells)
{
  m_cells.insert(m_cells.end(), a_cells.begin(), a_cells.end());
} // MeMultiPolyMesherImpl::AppendNewCells
//------------------------------------------------------------------------------
/// \brief Reports refine points that were not used.
/// \param a_io: MeMultiPolyMesherIo class that was provided to generate
/// the mesh.
/// \param a_usedPts: The locations of points that were processed.
//------------------------------------------------------------------------------
void MeMultiPolyMesherImpl::ReportUnusedRefinePts(const MeMultiPolyMesherIo& a_io,
                                                  const VecPt3d& a_usedPts)
{
  SetPt3d setPts(a_usedPts.begin(), a_usedPts.end());
  SetPt3d::iterator itEnd = setPts.end();
  VecStr strLoc;
  for (size_t i = 0; i < a_io.m_refPts.size(); ++i)
  {
    if (setPts.find(a_io.m_refPts[i].m_pt) == itEnd)
    {
      std::stringstream ss;
      ss << "(" << a_io.m_refPts[i].m_pt.x << ", " << a_io.m_refPts[i].m_pt.y << ")";
      strLoc.push_back(ss.str());
    }
  }
  if (!strLoc.empty())
  {
    std::string msg =
      "The following refine points were not included by the meshing process "
      "because the points are located outside of all polygons.";
    for (size_t i = 0; i < strLoc.size(); ++i)
    {
      msg += "\n" + strLoc[i];
    }
    XM_LOG(xmlog::warning, msg);
  }
} // ReportUnusedRefinePts

//////////////////////////////////////////////////////////////////////////////
/// \class MeMultiPolyMesher
/// \see MeMultiPolyMesherImpl
//////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief Creates a class
/// \return MeMultiPolyMesher.
//------------------------------------------------------------------------------
BSHP<MeMultiPolyMesher> MeMultiPolyMesher::New()
{
  BSHP<MeMultiPolyMesher> ret(new MeMultiPolyMesherImpl);
  return ret;
} // MePolysToUGrid::New

} // namespace xms

#if CXX_TEST
////////////////////////////////////////////////////////////////////////////////
// UNIT TESTS
////////////////////////////////////////////////////////////////////////////////

#include <xmsmesh/meshing/MeMultiPolyMesher.t.h>

#include <xmscore/testing/TestTools.h>

//----- Namespace declaration --------------------------------------------------

// namespace xms {
using namespace xms;

////////////////////////////////////////////////////////////////////////////////
/// \class MeMultiPolyMesherUnitTests
/// \brief Tests for MeMultiPolyMesher.
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief tests creating the class
//------------------------------------------------------------------------------
void MeMultiPolyMesherUnitTests::testCreateClass()
{
  BSHP<MeMultiPolyMesher> m = MeMultiPolyMesher::New();
  TS_ASSERT(m);
} // MeMultiPolyMesherUnitTests::testCreateClass
//------------------------------------------------------------------------------
/// \brief Tests checking for bad input: self-intersecting outer poly.
/// \verbatim
///              10              *
///                            / |
///               0    *----/----*
///                    | /
///             -10    *
///                    0--------100
/// \endverbatim
//------------------------------------------------------------------------------
void MeMultiPolyMesherUnitTests::testCheckForIntersections1()
{
  MeMultiPolyMesherIo input;
  input.m_checkTopology = true;
  {
    MePolyInput poly;
    poly.m_outPoly = {{0, 0, 0}, {100, 0, 0}, {100, 10, 0}, {0, -10, 0}};
    input.m_polys.push_back(poly);
  }

  // Mesh them
  BSHP<MeMultiPolyMesher> mesher = MeMultiPolyMesher::New();
  bool asserting = xmAsserting();
  xmAsserting() = false;
  bool rv = mesher->MeshIt(input);
  xmAsserting() = asserting;
  TS_ASSERT_EQUALS(rv, false);
  std::string errors = XmLog::Instance().GetAndClearStackStr();
  std::string expected =
    "---Error: Input polygon segments intersect. The segment defined by points 0 and 1 of outer "
    "polygon 0 intersects with the segment defined by points 2 and 3 of outer polygon 0.\n"
    "\n\n";

  TS_ASSERT_EQUALS(expected, errors);
} // MeMultiPolyMesherUnitTests::testCheckForIntersections1
//------------------------------------------------------------------------------
/// \brief Tests checking for bad input: self-intersecting inner poly.
/// \verbatim
///             100    *-------------*
///                    |           * |
///                    |         / | |
///                    | *----/----* |
///                    | | /         |
///                    | *           |
///               0    *-------------*
///                    0------------100
/// \endverbatim
//------------------------------------------------------------------------------
void MeMultiPolyMesherUnitTests::testCheckForIntersections2()
{
  MeMultiPolyMesherIo input;
  input.m_checkTopology = true;
  {
    MePolyInput poly;
    poly.m_outPoly = {{0, 0, 0}, {100, 0, 0}, {100, 100, 0}, {0, 100, 0}};
    poly.m_insidePolys.push_back({{10, 50}, {90, 50}, {90, 90}, {10, 10}});
    input.m_polys.push_back(poly);
  }

  // Mesh them
  BSHP<MeMultiPolyMesher> mesher = MeMultiPolyMesher::New();
  bool asserting = xmAsserting();
  xmAsserting() = false;
  bool rv = mesher->MeshIt(input);
  xmAsserting() = asserting;
  TS_ASSERT_EQUALS(rv, false);
  std::string errors = XmLog::Instance().GetAndClearStackStr();
  std::string expected =
    "---Error: Input polygon segments intersect. The segment defined by points 0 and 1 of inner "
    "polygon 0 of outer polygon 0 intersects with the segment defined by points 2 and 3 of inner "
    "polygon 0 of outer polygon 0.\n"
    "\n\n";

  TS_ASSERT_EQUALS(expected, errors);
} // MeMultiPolyMesherUnitTests::testCheckForIntersections2
//------------------------------------------------------------------------------
/// \brief Tests checking for bad input: inner poly intersects outer poly.
/// \verbatim
///             100    *-------------*
///                    |             |
///                    |          *--|--*
///                    |          |  |  |
///                    |          *--|--*
///               0    *-------------*
///                    0------------100
/// \endverbatim
//------------------------------------------------------------------------------
void MeMultiPolyMesherUnitTests::testCheckForIntersections3()
{
  MeMultiPolyMesherIo input;
  input.m_checkTopology = true;
  {
    MePolyInput poly;
    poly.m_outPoly = {{0, 0, 0}, {100, 0, 0}, {100, 100, 0}, {0, 100, 0}};
    poly.m_insidePolys.push_back({{90, 10}, {110, 10}, {110, 20}, {90, 20}});
    input.m_polys.push_back(poly);
  }

  // Mesh them
  BSHP<MeMultiPolyMesher> mesher = MeMultiPolyMesher::New();
  bool asserting = xmAsserting();
  xmAsserting() = false;
  bool rv = mesher->MeshIt(input);
  xmAsserting() = asserting;
  TS_ASSERT_EQUALS(rv, false);
  std::string errors = XmLog::Instance().GetAndClearStackStr();
  std::string expected =
    "---Error: Input polygon segments intersect. The segment defined by points 1 and 2 of outer "
    "polygon 0 intersects with the segment defined by points 0 and 1 of inner polygon 0 of outer "
    "polygon 0.\n"
    "Error: Input polygon segments intersect. The segment defined by points 1 and 2 of outer "
    "polygon 0 intersects with the segment defined by points 2 and 3 of inner polygon 0 of outer "
    "polygon 0.\n"
    "\n\n";

  TS_ASSERT_EQUALS(expected, errors);
} // MeMultiPolyMesherUnitTests::testCheckForIntersections3
//------------------------------------------------------------------------------
/// \brief Tests checking for bad input: 2 outer polys overlap
/// \verbatim
///                      *-------------*
///             100    *-|-----------* |
///                    | |           | |
///                    | |           | |
///                    | |           | |
///                    | |           | |
///                    | *-------------*
///               0    *-------------*
///                    0------------100
/// \endverbatim
//------------------------------------------------------------------------------
void MeMultiPolyMesherUnitTests::testCheckForIntersections4()
{
  MeMultiPolyMesherIo input;
  input.m_checkTopology = true;
  {
    MePolyInput poly;
    poly.m_outPoly = {{0, 0, 0}, {100, 0, 0}, {100, 100, 0}, {0, 100, 0}};
    input.m_polys.push_back(poly);
    poly.m_outPoly = {{10, 10, 0}, {110, 10, 0}, {110, 110, 0}, {10, 110, 0}};
    input.m_polys.push_back(poly);
  }

  // Mesh them
  BSHP<MeMultiPolyMesher> mesher = MeMultiPolyMesher::New();
  bool asserting = xmAsserting();
  xmAsserting() = false;
  bool rv = mesher->MeshIt(input);
  xmAsserting() = asserting;
  TS_ASSERT_EQUALS(rv, false);
  std::string errors = XmLog::Instance().GetAndClearStackStr();
  std::string expected =
    "---Error: Input polygon segments intersect. The segment defined by points 1 and 2 of outer "
    "polygon 0 intersects with the segment defined by points 0 and 1 of outer polygon 1.\n"
    "Error: Input polygon segments intersect. The segment defined by points 2 and 3 of outer "
    "polygon 0 intersects with the segment defined by points 3 and 0 of outer polygon 1.\n"
    "\n\n";

  TS_ASSERT_EQUALS(expected, errors);
} // MeMultiPolyMesherUnitTests::testCheckForIntersections4
//------------------------------------------------------------------------------
/// \brief Tests checking for bad input: 2 inner polys overlap
/// \verbatim
///             100-   *--------------*
///                    |    *-----*   |
///                    |  *-|---* |   |
///                    |  | |   | |   |
///                    |  | *-----*   |
///                    |  *-----*     |
///               0-   *--------------*
///
///                    0------------100
/// \endverbatim
//------------------------------------------------------------------------------
void MeMultiPolyMesherUnitTests::testCheckForIntersections5()
{
  MeMultiPolyMesherIo input;
  input.m_checkTopology = true;
  {
    MePolyInput poly;
    poly.m_outPoly = {{0, 0, 0}, {100, 0, 0}, {100, 100, 0}, {0, 100, 0}};
    poly.m_insidePolys.push_back({{10, 10}, {60, 10}, {60, 60}, {10, 60}});
    poly.m_insidePolys.push_back({{40, 40}, {90, 40}, {90, 90}, {40, 90}});
    input.m_polys.push_back(poly);
  }

  // Mesh them
  BSHP<MeMultiPolyMesher> mesher = MeMultiPolyMesher::New();
  bool asserting = xmAsserting();
  xmAsserting() = false;
  bool rv = mesher->MeshIt(input);
  xmAsserting() = asserting;
  TS_ASSERT_EQUALS(rv, false);
  std::string errors = XmLog::Instance().GetAndClearStackStr();
  std::string expected =
    "---Error: Input polygon segments intersect. The segment defined by points 1 and 2 of inner "
    "polygon 0 of outer polygon 0 intersects with the segment defined by points 0 and 1 of inner "
    "polygon 1 of outer polygon 0.\n"
    "Error: Input polygon segments intersect. The segment defined by points 2 and 3 of inner "
    "polygon 0 of outer polygon 0 intersects with the segment defined by points 3 and 0 of inner "
    "polygon 1 of outer polygon 0.\n"
    "\n\n";

  TS_ASSERT_EQUALS(expected, errors);
} // MeMultiPolyMesherUnitTests::testCheckForIntersections5

//} // namespace xms

#endif // CXX_TEST
