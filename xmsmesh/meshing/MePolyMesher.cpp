//------------------------------------------------------------------------------
/// \file
/// \ingroup meshing
/// \copyright (C) Copyright Aquaveo 2018. Distributed under the xmsng
///  Software License, Version 1.0. (See accompanying file
///  LICENSE_1_0.txt or copy at http://www.aquaveo.com/xmsng/LICENSE_1_0.txt)
//------------------------------------------------------------------------------

//----- Included files ---------------------------------------------------------

// 1. Precompiled header

// 2. My own header
#include <xmsmesh/meshing/MePolyMesher.h>

// 3. Standard library headers
#include <fstream>

// 4. External library headers
#pragma warning(push)
#pragma warning(disable : 4512) // boost code: no assignment operator
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#pragma warning(pop)

// 5. Shared code headers
#include <xmscore/math/math.h>
#include <xmscore/misc/XmConst.h>
#include <xmscore/misc/XmError.h>
#include <xmscore/misc/boost_defines.h> // BSHP
#include <xmscore/misc/Observer.h>
#include <xmscore/misc/XmLog.h>
#include <xmscore/stl/set.h>
#include <xmsinterp/geometry/geoms.h>
#include <xmsinterp/interpolate/InterpBase.h>
#include <xmsinterp/triangulate/detail/TrAutoFixFourTrianglePts.h>
#include <xmsinterp/triangulate/detail/TrOuterTriangleDeleter.h>
#include <xmsinterp/triangulate/TrTriangulatorPoints.h>
#include <xmsinterp/triangulate/TrBreaklineAdder.h>
#include <xmsinterp/triangulate/TrTin.h>
#include <xmsmesh/meshing/detail/MePolyPaverToMeshPts.h>
#include <xmsmesh/meshing/detail/MePolyPatcher.h>
#include <xmsmesh/meshing/detail/MeRefinePtsToPolys.h>
#include <xmsmesh/meshing/detail/MeRelaxer.h>
#include <xmsmesh/meshing/MeMeshUtils.h>
#include <xmsmesh/meshing/MeMultiPolyMesherIo.h>
#include <xmsmesh/meshing/MePolyRedistributePts.h>

// 6. Non-shared code headers

//----- Forward declarations ---------------------------------------------------

//----- External globals -------------------------------------------------------

//----- Namespace declaration --------------------------------------------------

namespace xms
{
//----- Constants / Enumerations -----------------------------------------------

//----- Classes / Structs ------------------------------------------------------
typedef boost::unordered_map<std::pair<double, double>, int> PtHash; ///< typdef for a short name

////////////////////////////////////////////////////////////////////////////////
class MePolyMesherImpl : public MePolyMesher
{
public:
  /// enumeration for boundary processing
  enum BoundaryEnum {
    BE_UNVISITED = 0,
    BE_VISITED = 1,
    BE_OUTSIDE = 2,
    BE_INSIDE = 4,
    BE_ONBOUNDARY = 8
  };

  MePolyMesherImpl();
  virtual ~MePolyMesherImpl();

public:
  virtual bool MeshIt(const MeMultiPolyMesherIo& a_input,
                      size_t a_polyIdx,
                      VecPt3d& a_points,
                      VecInt& a_triangles,
                      VecInt& a_cells) override;
  virtual bool MeshIt(const VecPt3d& a_outPoly,
                      const VecPt3d2d& a_inPolys,
                      double a_bias,
                      VecPt3d& a_points,
                      VecInt& a_triangles);

  //------------------------------------------------------------------------------
  /// \brief sets the observer class to get feedback on the meshing process
  /// \param a_: The observer.
  //------------------------------------------------------------------------------
  virtual void SetObserver(BSHP<Observer> a_) override { m_observer = a_; }
  virtual void GetProcessedRefinePts(std::vector<Pt3d>& a_pts) override;

  void TestWithPoints(const VecInt& a_outPoly,
                      const VecInt2d& a_inPolys,
                      const VecPt3d& a_points,
                      VecInt& a_triangles);

private:
  bool MeshFromInputs(std::vector<Pt3d>& a_points, VecInt& a_triangles, VecInt& a_cells);
  void SortPoly(VecPt3d& a_outPoly);
  bool ComputeExtents(Pt3d& a_mn, Pt3d& a_mx);
  void ComputeTolerance();
  void GenerateMeshPts();
  void ProcessBoundaryPtsFlaggedToRemove();
  void Triangulate();
  void FindAllPolyPointIdxs();
  void FindPolyPointIdxs(const VecPt3d& a_poly, VecInt& a_polyPtIdxs);
  void AddBreaklines();
  void DeleteTrianglesOutsidePolys();
  void AutoFixFourTrianglePts();
  void Relax();
  void ExportTinForDebug();

private:
  VecPt3d m_outPoly;      ///< outer boundary of polygon being meshed
  VecPt3d2d m_inPolys;    ///< inside boundaries or holes in polygon
  VecPt3d2d m_refPtPolys; ///< refine points
  VecPt3d m_refMeshPts;   ///< refine point that have been made into mesh nodes
  VecPt3d
    m_refPtsTooClose;     ///< refine points that can not be honored because of distance constraints
  BSHP<VecPt3d> m_points; ///< resulting mesh nodes
  BSHP<TrTin> m_tin;      ///< triangles that become mesh elements
  BSHP<MePolyRedistributePts>
    m_redist; ///< class for performing redistribution of points along a polyline
  BSHP<MePolyPaverToMeshPts> m_polyPaver; ///< class for paving from a polygon definition
  BSHP<MeRefinePtsToPolys>
    m_refineToPolys;         ///< class for creating polygon from refine point information
  BSHP<MeRelaxer> m_relaxer; ///< class for relaxing the location of mesh nodes
  VecInt m_cells;            ///< Cells generated by the Patcher
  double m_bias; ///< factor that affects how quickly the size of elements transitions in the mesh
  VecInt m_outPolyPtIdxs;    ///< indices to the mesh points that match the outer polygon boundary
  VecInt2d m_inPolyPtIdxs;   ///< indices to the mesh points that match the inner polygon boundaries
  VecInt m_refPtIdxs;        ///< indices to points mesh point that match refine points
  VecInt m_polyCorners;      ///< corner indexes used with mesh patch generation
  double m_xyTol;            ///< xy tolerance used in geometric comparisons
  bool m_testing;            ///< flag to indicate if we are testing the class
  BSHP<Observer> m_observer; ///< observer to send feedback on the meshing progress
  Pt3d m_min;                ///< min xy bound
  Pt3d m_max;                ///< max xy bound
  PtHash m_ptHash;           ///< hash for point locations
  BSHP<InterpBase> m_elev;   ///< interpolator to assign elevations to mesh points
  int m_polyId;              ///< id of the polygon
  VecPt3d m_boundPtsToRemove; ///< boundary points to remove after the paving process is complete
  bool m_removeInternalFourTrianglePts =
    false; ///< flag to indicate the removal of internal pts connected to 4 triangles will occur
};         // class MePolyMesherImpl

//----- Internal functions -----------------------------------------------------

//----- Class / Function definitions -------------------------------------------

//------------------------------------------------------------------------------
/// \brief Creates a polymesher class
/// \return MePolyMesher.
//------------------------------------------------------------------------------
BSHP<MePolyMesher> MePolyMesher::New()
{
  BSHP<MePolyMesher> polyMesher(new MePolyMesherImpl);
  return polyMesher;
} // MePolyMesher::New
//------------------------------------------------------------------------------
/// \brief
//------------------------------------------------------------------------------
MePolyMesher::MePolyMesher()
{
} // MePolyMesher::MePolyMesher
//------------------------------------------------------------------------------
/// \brief
//------------------------------------------------------------------------------
MePolyMesher::~MePolyMesher()
{
} // MePolyMesher::~MePolyMesher

////////////////////////////////////////////////////////////////////////////////
/// \class MePolyMesherImpl
/// \brief Creates a mesh inside a polygon.
///
/// The polygon may contain holes. The mesh is returned as a vector of points
/// and a vector of triangles.
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief
//------------------------------------------------------------------------------
MePolyMesherImpl::MePolyMesherImpl()
: MePolyMesher()
, m_outPoly()
, m_inPolys()
, m_points(new VecPt3d())
, m_tin(TrTin::New())
, m_redist(MePolyRedistributePts::New())
, m_polyPaver(MePolyPaverToMeshPts::New())
, m_refineToPolys(MeRefinePtsToPolys::New())
, m_relaxer(MeRelaxer::New())
, m_bias(1.0)
, m_outPolyPtIdxs()
, m_inPolyPtIdxs()
, m_refPtIdxs()
, m_xyTol(1e-9)
, m_testing(false)
, m_observer()
, m_polyId(-1)
{
} // MePolyMesherImpl::MePolyMesherImpl
//------------------------------------------------------------------------------
/// \brief
//------------------------------------------------------------------------------
MePolyMesherImpl::~MePolyMesherImpl()
{
} // MePolyMesherImpl::~MePolyMesherImpl
//------------------------------------------------------------------------------
/// \brief Perform MESH_PAVE, MESH_SPAVE, MESH_PATCH meshing on a polygon.
/// \param[in]  a_input:     Meshing input: polygons and optional inputs
/// \param[in]  a_polyIdx:   Index to the polygon in a_input to mesh
/// \param[out] a_points:    Computed mesh points.
/// \param[out] a_triangles: Computed mesh triangles from paving.
/// \param[out] a_cells:     Computed mesh cells from patch.
/// \return true if no errors encountered.
//------------------------------------------------------------------------------
bool MePolyMesherImpl::MeshIt(const MeMultiPolyMesherIo& a_input,
                              size_t a_polyIdx,
                              VecPt3d& a_points,
                              VecInt& a_triangles,
                              VecInt& a_cells)
{
  // reinitialize some internal classes
  m_tin = TrTin::New();
  m_redist = MePolyRedistributePts::New();
  m_polyPaver = MePolyPaverToMeshPts::New();
  m_refineToPolys = MeRefinePtsToPolys::New();
  m_relaxer = MeRelaxer::New();

  // outer polygons
  const MePolyInput& polyInput = a_input.m_polys[a_polyIdx];

  m_polyId = polyInput.m_polyId;
  m_outPoly = polyInput.m_outPoly;
  ComputeTolerance();

  // holes inside the outer polygons
  m_inPolys = polyInput.m_insidePolys;
  // bias term
  m_bias = polyInput.m_bias;
  // refine pts
  m_refineToPolys->SetRefinePoints(a_input.m_refPts, m_xyTol);
  m_refineToPolys->RefPtsAsPolys(polyInput.m_polyId, m_outPoly, m_inPolys, m_refPtPolys, m_refMeshPts,
                                 m_refPtsTooClose);
  // size function
  if (!polyInput.m_sizeFunction)
  {
    VecPt3d2d inPolys(m_inPolys);
    inPolys.insert(inPolys.end(), m_refPtPolys.begin(), m_refPtPolys.end());
    m_redist->SetSizeFuncFromPoly(m_outPoly, inPolys, m_bias);
  }
  else
  {
    m_redist->SetSizeFunc(polyInput.m_sizeFunction);
  }
  if (polyInput.m_constSizeFunction != -1.0)
  {
    m_redist->SetConstantSizeFunc(polyInput.m_constSizeFunction);
    if (polyInput.m_constSizeBias != -1.0)
      m_redist->SetConstantSizeBias(polyInput.m_constSizeBias);
  }
  m_polyPaver->SetRedistributor(m_redist);

  // patch
  if (!polyInput.m_polyCorners.empty())
    m_polyCorners = polyInput.m_polyCorners;
  else
    SortPoly(m_outPoly);

  m_elev = polyInput.m_elevFunction;
  m_boundPtsToRemove = polyInput.m_boundPtsToRemove;
  m_removeInternalFourTrianglePts = polyInput.m_removeInternalFourTrianglePts;
  return MeshFromInputs(a_points, a_triangles, a_cells);
} // MePolyMesherImpl::MeshIt
//------------------------------------------------------------------------------
/// \brief Perform MESH_PAVE, MESH_SPAVE meshing on a polygon. Return the
///        mesh points through a_points and triangles through a_triangles.
/// \param[in]  a_outPoly: Outer polygon. Clockwise. 1st pt != last
/// \param[in]  a_inPolys: Inner polygons (holes). Counter clockwise. 1st pt != last
/// \param[in]  a_bias:
/// \param[out] a_points: Computed mesh points.
/// \param[out] a_triangles: Computed mesh triangles.
/// \return true on success.
//------------------------------------------------------------------------------
bool MePolyMesherImpl::MeshIt(const VecPt3d& a_outPoly,
                              const VecPt3d2d& a_inPolys,
                              double a_bias,
                              VecPt3d& a_points,
                              VecInt& a_triangles)
{
  m_outPoly = a_outPoly;
  ComputeTolerance();
  SortPoly(m_outPoly);
  m_inPolys = a_inPolys;
  m_bias = a_bias;
  VecInt cells;
  return MeshFromInputs(a_points, a_triangles, cells);
} // MePolyMesherImp::MeshIt
//------------------------------------------------------------------------------
/// \brief Start polygon with lower left most point so results are consistent.
/// \param[in,out] a_outPoly: The polygon.
//------------------------------------------------------------------------------
void MePolyMesherImpl::SortPoly(VecPt3d& a_outPoly)
{
  // Find the lower left most point
  double minX(XM_DBL_HIGHEST), minY(XM_DBL_HIGHEST);
  size_t minIdx = 0;
  for (size_t i = 0; i < a_outPoly.size(); ++i)
  {
    if (a_outPoly[i].x < minX)
    {
      minX = a_outPoly[i].x;
      minY = a_outPoly[i].y;
      minIdx = i;
    }
    else if (EQ_TOL(a_outPoly[i].x, minX, m_xyTol))
    {
      if (a_outPoly[i].y < minY)
      { // X is same as min. Check y
        minX = a_outPoly[i].x;
        minY = a_outPoly[i].y;
        minIdx = i;
      }
    }
  }
  std::rotate(a_outPoly.begin(), a_outPoly.begin() + minIdx, a_outPoly.end());
} // MePolyMesherImpl::SortPoly
//------------------------------------------------------------------------------
/// \brief Gets the refine points that were inside the polygon, both points
/// that are included in the meshing process and those that were not.
/// \param a_pts Locations of refine points used inside of this polygon.
//------------------------------------------------------------------------------
void MePolyMesherImpl::GetProcessedRefinePts(std::vector<Pt3d>& a_pts)
{
  a_pts.resize(0);
  a_pts = m_refMeshPts;
  a_pts.insert(a_pts.end(), m_refPtsTooClose.begin(), m_refPtsTooClose.end());
} // MePolyMesherImpl::GetProcessedRefinePts
//------------------------------------------------------------------------------
/// \brief Creates the mesh from inputs that have set member variables in the
/// class.
/// \param[out] a_points:    Points filled by meshing.
/// \param[out] a_triangles: Triangles created by paving meshing.
/// \param[out] a_cells:     Cells created by patch meshing.
/// \return true if no errors encountered.
//------------------------------------------------------------------------------
bool MePolyMesherImpl::MeshFromInputs(VecPt3d& a_points, VecInt& a_triangles, VecInt& a_cells)
{
  try
  {
    GenerateMeshPts();
    if (m_cells.empty())
    { // paving
      ProcessBoundaryPtsFlaggedToRemove();
      Triangulate();
      FindAllPolyPointIdxs();
      AutoFixFourTrianglePts();
      AddBreaklines();
      DeleteTrianglesOutsidePolys();
      Relax();
      // Re-add breaklines and delete outer polys because relaxing can swap edges
      AddBreaklines();
      DeleteTrianglesOutsidePolys();
      a_triangles.swap(m_tin->Triangles());
      a_cells.resize(0);
    }
    else
    { // patch
      a_cells.swap(m_cells);
      a_triangles.resize(0);
    }
    a_points.swap(*m_points);
    m_tin->Clear();
    m_polyCorners.clear();
    if (m_elev)
      for (Pt3d& p : a_points)
        p.z = (double)m_elev->InterpToPt(p);
  }
  catch (std::exception& e)
  {
    std::string msg = e.what();
    meModifyMessageWithPolygonId(m_polyId, msg);
    XM_LOG(xmlog::error, msg);
    return false;
  }
  return true;
} // MePolyMesherImpl::MeshFromInputs
//------------------------------------------------------------------------------
/// \brief Used only for testing. Test the class by supplying the polygons
///        and mesh points.
/// \param[in] a_outPoly: Outer polygon.
/// \param[in] a_inPolys: Inner polygons.
/// \param[in] a_points: The points.
/// \param[out] a_triangles: The created triangles.
//------------------------------------------------------------------------------
void MePolyMesherImpl::TestWithPoints(const VecInt& a_outPoly,
                                      const VecInt2d& a_inPolys,
                                      const VecPt3d& a_points,
                                      VecInt& a_triangles)
{
  m_testing = true;

  // Convert polys from idxs to Pt3d
  VecPt3d outPoly(a_outPoly.size());
  for (size_t i = 0; i < a_outPoly.size(); ++i)
  {
    outPoly[i] = a_points[a_outPoly[i]];
  }
  VecPt3d2d inPolys(a_inPolys.size());
  for (size_t j = 0; j < a_inPolys.size(); ++j)
  {
    inPolys[j].resize(a_inPolys[j].size());
    for (size_t i = 0; i < a_inPolys[j].size(); ++i)
    {
      inPolys[j][i] = a_points[a_inPolys[j][i]];
    }
  }

  *m_points = a_points;
  VecPt3d points;
  MeshIt(outPoly, inPolys, 1.0, points, a_triangles);
} // MePolyMesherImpl::TestWithPoints
//------------------------------------------------------------------------------
/// \brief Computes the extents (min, max) of the polygon.
/// \param[out] a_mn The minimum coordinates of the polygon extents.
/// \param[out] a_mx The maximum coordinates of the polygon extents.
/// \return true if there are any points and the extents were computed.
//------------------------------------------------------------------------------
bool MePolyMesherImpl::ComputeExtents(Pt3d& a_mn, Pt3d& a_mx)
{
  if (!m_outPoly.empty())
  {
    a_mn = XM_DBL_HIGHEST;
    a_mx = XM_DBL_LOWEST;
    for (size_t i = 0; i < m_outPoly.size(); ++i)
    {
      gmAddToExtents(m_outPoly[i], a_mn, a_mx);
    }
    return true;
  }
  return false;
} // MePolyMesherImpl::ComputeExtents
//------------------------------------------------------------------------------
/// \brief Computes a tolerance to use based on point extents.
//------------------------------------------------------------------------------
void MePolyMesherImpl::ComputeTolerance()
{
  if (ComputeExtents(m_min, m_max))
  {
    m_xyTol = gmComputeXyTol(m_min, m_max);
  }
} // MePolyMesherImpl::ComputeTolerance
//------------------------------------------------------------------------------
/// \brief Creates the points in interior of the input polygon that are used
/// to create cells.
//------------------------------------------------------------------------------
void MePolyMesherImpl::GenerateMeshPts()
{
  if (m_testing)
    return;
  if (!m_polyCorners.empty())
  {
    BSHP<MePolyPatcher> patcher(MePolyPatcher::New());
    patcher->MeshIt(m_polyId, m_outPoly, m_polyCorners, m_xyTol, *m_points, m_cells);
  }
  else
  {
    m_polyPaver->SetObserver(m_observer);
    // add the refine point polygons
    VecPt3d2d inPolys(m_inPolys);
    inPolys.insert(inPolys.end(), m_refPtPolys.begin(), m_refPtPolys.end());
    // generate mesh points
    m_polyPaver->PolyToMeshPts(m_outPoly, inPolys, m_bias, m_xyTol, *m_points);

    // add the mesh refine points
    m_points->insert(m_points->end(), m_refMeshPts.begin(), m_refMeshPts.end());
  }
} // MePolyMesherImpl::GenerateMeshPts
//------------------------------------------------------------------------------
/// \brief Remove boundary points
//------------------------------------------------------------------------------
void MePolyMesherImpl::ProcessBoundaryPtsFlaggedToRemove()
{
  if (m_boundPtsToRemove.empty())
    return;
  boost::unordered_set<std::pair<double, double>> ptSet;
  for (auto& p : m_boundPtsToRemove)
  {
    std::pair<double, double> pp(p.x, p.y);
    ptSet.insert(pp);
  }
  auto itEnd = ptSet.end();
  for (size_t i = m_outPoly.size(); i > 0; --i)
  {
    size_t ix = i - 1;
    std::pair<double, double> pp(m_outPoly[ix].x, m_outPoly[ix].y);
    auto it = ptSet.find(pp);
    if (it != itEnd)
    {
      m_outPoly.erase(m_outPoly.begin() + ix);
    }
  }
  VecPt3d& pts(*m_points);
  for (size_t i = pts.size(); i > 0; --i)
  {
    size_t ix = i - 1;
    std::pair<double, double> pp(pts[ix].x, pts[ix].y);
    auto it = ptSet.find(pp);
    if (it != itEnd)
    {
      pts.erase(pts.begin() + ix);
    }
  }
} // MePolyMesherImpl::ProcessBoundaryPtsFlaggedToRemove
//------------------------------------------------------------------------------
/// \brief Triangulate the mesh points.
//------------------------------------------------------------------------------
void MePolyMesherImpl::Triangulate()
{
  m_tin->SetPoints(m_points);
  TrTriangulatorPoints client(*m_points.get(), m_tin->Triangles(), &m_tin->TrisAdjToPts());
  client.SetObserver(m_observer);
  bool rv = client.Triangulate();
  if (!rv)
  {
    throw std::runtime_error("Error triangulating mesh points.");
  }
} // MePolyMesherImpl::Triangulate
//------------------------------------------------------------------------------
/// \brief Exports the tin for debugging purposes.
//------------------------------------------------------------------------------
void MePolyMesherImpl::ExportTinForDebug()
{
#if _DEBUG
  std::ofstream ofs;
  ofs.open("C:\\temp\\TrBreaklineAdder.tin");
  m_tin->ExportTinFile(ofs);
#endif
} // MePolyMesherImpl::ExportTinForDebug
//------------------------------------------------------------------------------
/// \brief Add the polys as breaklines in the tin.
//------------------------------------------------------------------------------
void MePolyMesherImpl::AddBreaklines()
{
  BSHP<TrBreaklineAdder> adder = TrBreaklineAdder::New();
  adder->SetObserver(m_observer);
  adder->SetTin(m_tin, m_xyTol);

  // Close the polys
  VecInt outPolyPtIdxs = m_outPolyPtIdxs;
  outPolyPtIdxs.push_back(m_outPolyPtIdxs.front());
  VecInt2d inPolyPtIdxs = m_inPolyPtIdxs;
  for (size_t i = 0; i < m_inPolyPtIdxs.size(); ++i)
  {
    inPolyPtIdxs[i].push_back(m_inPolyPtIdxs[i].front());
  }

  adder->AddBreakline(outPolyPtIdxs);
  adder->AddBreaklines(inPolyPtIdxs);
} // MePolyMesherImpl::AddBreaklines
//------------------------------------------------------------------------------
/// \brief Delete triangles outside the polygon boundary or in polygon holes.
///
/// The polygon may be concave or have holes so we may need to delete some
/// triangles. This method finds triangles along the polygon and classifies
/// them as in or out and deletes all that are out.
//------------------------------------------------------------------------------
void MePolyMesherImpl::DeleteTrianglesOutsidePolys()
{
  if (m_observer)
  {
    m_observer->BeginOperationString("Deleting triangles outside of polygon.");
  }

  // Close the polygons by making the last point the same as the first and
  // create one 2D array for outer and all inner
  VecInt2d allPolys;
  allPolys.push_back(m_outPolyPtIdxs);
  allPolys.back().push_back(m_outPolyPtIdxs.front());
  for (size_t i = 0; i < m_inPolyPtIdxs.size(); ++i)
  {
    allPolys.push_back(m_inPolyPtIdxs[i]);
    allPolys.back().push_back(m_inPolyPtIdxs[i].front());
  }

  BSHP<TrOuterTriangleDeleter> deleter = TrOuterTriangleDeleter::New();
  deleter->Delete(allPolys, m_tin);

  if (m_observer)
  {
    m_observer->EndOperation();
  }
} // MePolyMesherImpl::DeleteTrianglesOutsidePolys
//------------------------------------------------------------------------------
/// \brief Delete internal points that are only connected to 4 triangles and
/// retriangulates.
//------------------------------------------------------------------------------
void MePolyMesherImpl::AutoFixFourTrianglePts()
{
  if (!m_removeInternalFourTrianglePts)
    return;
  BSHP<TrAutoFixFourTrianglePts> fixer = TrAutoFixFourTrianglePts::New();
  fixer->SetObserver(m_observer);

  // can't delete boundary pts or refine points
  VecInt noDeletePts(m_outPolyPtIdxs);
  for (size_t i = 0; i < m_inPolys.size(); ++i)
  {
    for (auto& ix : m_inPolyPtIdxs[i])
      noDeletePts.push_back(ix);
  }
  for (auto& ix : m_refPtIdxs)
    noDeletePts.push_back(ix);
  fixer->SetUndeleteablePtIdxs(noDeletePts);

  size_t ptsBefore(m_tin->Points().size());
  fixer->Fix(m_tin);
  if (ptsBefore != m_tin->Points().size())
    FindAllPolyPointIdxs();
} // MePolyMesherImpl::AutoFixFourTrianglePts
//------------------------------------------------------------------------------
/// \brief Relaxes the mesh points for better element quality
//------------------------------------------------------------------------------
void MePolyMesherImpl::Relax()
{
  // Put polygons into one vector
  VecInt fixedPoints(m_outPolyPtIdxs.begin(), m_outPolyPtIdxs.end());
  for (size_t i = 0; i < m_inPolyPtIdxs.size(); ++i)
  {
    fixedPoints.insert(fixedPoints.end(), m_inPolyPtIdxs[i].begin(), m_inPolyPtIdxs[i].end());
  }
  // add refine points to the fixed
  fixedPoints.insert(fixedPoints.end(), m_refPtIdxs.begin(), m_refPtIdxs.end());
  m_relaxer->Relax(fixedPoints, m_tin);
} // MePolyMesherImpl::Relax
//------------------------------------------------------------------------------
/// \brief Find the indices of the poly points among m_points. They will most
///        likely be at the front of m_points so we don't use an rtree to
///        find them.
//------------------------------------------------------------------------------
void MePolyMesherImpl::FindAllPolyPointIdxs()
{
  m_ptHash.clear();
  std::pair<std::pair<double, double>, int> pd;
  size_t nPts(m_points->size());
  Pt3d* pts(nPts > 0 ? &(*m_points)[0] : NULL);
  for (size_t i = 0; i < nPts; ++i)
  {
    pd.first.first = pts[i].x;
    pd.first.second = pts[i].y;
    pd.second = (int)i;
    m_ptHash.insert(pd);
  }
  FindPolyPointIdxs(m_outPoly, m_outPolyPtIdxs);
  m_inPolyPtIdxs.resize(m_inPolys.size(), VecInt());
  for (size_t i = 0; i < m_inPolys.size(); ++i)
  {
    FindPolyPointIdxs(m_inPolys[i], m_inPolyPtIdxs[i]);
  }
  // get the refine point idxs
  m_refPtIdxs.resize(0);
  FindPolyPointIdxs(m_refMeshPts, m_refPtIdxs);
  for (size_t i = 0; i < m_refPtPolys.size(); ++i)
  {
    VecInt tmp;
    FindPolyPointIdxs(m_refPtPolys[i], tmp);
    m_refPtIdxs.insert(m_refPtIdxs.end(), tmp.begin(), tmp.end());
  }
} // MePolyMesherImpl::FindAllPolyPointIdxs
//------------------------------------------------------------------------------
/// \brief See FindAllPolyPointIdxs.
/// \param[in] a_poly: The polygon.
/// \param[out] a_polyPtIdxs: Polygon point indices.
//------------------------------------------------------------------------------
void MePolyMesherImpl::FindPolyPointIdxs(const VecPt3d& a_poly, VecInt& a_polyPtIdxs)
{
  a_polyPtIdxs.resize(0);
  size_t nPts(a_poly.size());
  auto itEnd = m_ptHash.end();
  auto it = m_ptHash.begin();
  std::pair<double, double> pd;
  for (size_t i = 0; i < nPts; ++i)
  {
    pd.first = a_poly[i].x;
    pd.second = a_poly[i].y;
    it = m_ptHash.find(pd);
    if (it == itEnd)
    {
      throw std::runtime_error("Polygon point not in mesh points.");
    }
    else
      a_polyPtIdxs.push_back(it->second);
  }
} // MePolyMesherImpl::FindPolyPointIdxs

} // namespace xms

#if CXX_TEST
////////////////////////////////////////////////////////////////////////////////
// UNIT TESTS
////////////////////////////////////////////////////////////////////////////////

#include <xmsmesh/meshing/MePolyMesher.t.h>

#include <xmscore/testing/TestTools.h>

//----- Namespace declaration --------------------------------------------------

// namespace xms {
using namespace xms;

namespace
{
//------------------------------------------------------------------------------
/// \brief
//------------------------------------------------------------------------------
static VecPt3d iArrayToVecPt3d(int* a_array, int a_size)
{
  VecPt3d v(a_size / 2);
  for (int i = 0; i < a_size; i += 2)
  {
    v[i / 2].x = a_array[i];
    v[i / 2].y = a_array[i + 1];
  }
  return v;
} // iArrayToVecPt3d

} // namespace unnamed

////////////////////////////////////////////////////////////////////////////////
/// \class MePolyMesherUnitTests
/// \brief Tests for MePolyMesher.
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief Tests the part of MePolyMesher that triangulates, adds breaklines,
///        and removes outer triangles.
/// \verbatim
//  Before
//
//  10-  17-----18------19------20------21
//    |  |\  8  /|\ 11  /|\ 22  / \ 26  /|
//    |  | \   / | \   / | \   /   \   //|
//    |  |  \ /  |  \ /30|23\ / 27  \ / /|
//   5-  |6 13 12|9 14   |  15------16  /|
//    |  |  / \  |  /|\  |  / \ 25  /|28/|
//    |  | /   \ | / | \ | /   \   / | / |
//    |  |/  7  \|/  |  \|/ 24  \ /29| / |
//   0-  9------10 10|3 11------12   | / |
//    |  |\ 13  / \  |  / \ 15  / \  |/  |
//    |  | \   /   \ | /   \   /   \ |/  |
//    |  |  \ /  5  \|/ 16  \ / 21  \|/  |
//  -5-  | 0 5-------6-------7-------8 18|
//    |  |  / \  4  / \ 14  / \ 19  / \  |
//    |  | /   \   /   \   /   \   /   \ |
//    |  |/  1  \ /  2  \ / 20  \ / 17  \|
// -10-  0-------1-------2-------3-------4
//
//       |-------|-------|-------|-------|
//       0      10      20      30      40
//
//  After swapping
//
//  10-  17-----18------19------20------21
//    |  |\  8  / \ 11  /|\ 22  / \ 26  /|
//    |  | \   /   \   / | \   /   \   //|
//    |  |  \ / 12  \ /30|23\ / 27  \ / /|
//   5-  |6 13------14   |  15------16  /|
//    |  |  /|\  9  /|\  |  /|\ 25  /|28/|
//    |  | / | \   / | \ | / | \   / | / |
//    |  |/  |  \ /  |  \|/  |24\ /29| / |
//   0-  9 13|7 10 10|3 11 15|  12   | / |
//    |  |\  |  / \  |  / \  |  / \  |/  |
//    |  | \ | /   \ | /   \ | /   \ |/  |
//    |  |  \|/  5  \|/ 16  \|/ 21  \|/  |
//  -5-  | 0 5-------6-------7-------8 18|
//    |  |  / \  4  / \ 14  / \ 19  / \  |
//    |  | /   \   /   \   /   \   /   \ |
//    |  |/  1  \ /  2  \ / 20  \ / 17  \|
// -10-  0-------1-------2-------3-------4
//
//       |-------|-------|-------|-------|
//       0      10      20      30      40
//
//  After removing triangles (*** means a hole)
//
//  10-  17-----18------19------20------21
//    |  |\  7  / \  9  /|\ 17  / \ 21  /
//    |  | \   /   \   / | \   /   \   /
//    |  |  \ / 10  \ /24|18\ / 22  \ /
//   5-  |6 13------14   |  15------16
//    |  |  /|******/|\  |  /|\ 20  /|
//    |  | / |*****/ | \ | /*| \   / |
//    |  |/  |****/  |  \|/**|19\ /23|
//   0-  9 11|**10  8|3 11***|  12   |
//    |  |\  |**/ \  |  /****|  / \  |
//    |  | \ |*/   \ | /*****| /   \ |
//    |  |  \|/  5  \|/******|/ 16  \|
//  -5-  | 0 5-------6-------7-------8
//    |  |  / \  4  / \ 12  / \ 14  / \
//    |  | /   \   /   \   /   \   /   \
//    |  |/  1  \ /  2  \ / 15  \ / 13  \
// -10-  0-------1-------2-------3-------4
//
//       |-------|-------|-------|-------|
//       0      10      20      30      40
///\endverbatim
//------------------------------------------------------------------------------
void MePolyMesherUnitTests::test1()
{
  // Set up tin
  int meshPtsA[] = {0,  -10, 10, -10, 20, -10, 30, -10, 40, -10, 5,  -5, 15, -5, 25,
                    -5, 35,  -5, 0,   0,  10,  0,  20,  0,  30,  0,  5,  5,  15, 5,
                    25, 5,   35, 5,   0,  10,  10, 10,  20, 10,  30, 10, 40, 10};
  VecPt3d points = iArrayToVecPt3d(meshPtsA, 44);

  // Set up polygons
  int outPolyA[] = {9, 17, 18, 19, 20, 21, 16, 8, 4, 3, 2, 1, 0};
  VecInt outPoly(&outPolyA[0], &outPolyA[13]);
  int inPoly1[] = {14, 13, 5, 10};
  int inPoly2[] = {7, 15, 11, 6};
  VecInt2d inPolys;
  inPolys.push_back(VecInt(&inPoly1[0], &inPoly1[4]));
  inPolys.push_back(VecInt(&inPoly2[0], &inPoly2[4]));

  // Create the mesher, set the test points, and mesh it

  BSHP<MePolyMesher> mesher = MePolyMesher::New();
  BSHP<MePolyMesherImpl> mesherImp = BDPC<MePolyMesherImpl>(mesher);
  VecInt triangles;
  mesherImp->TestWithPoints(outPoly, inPolys, points, triangles);

  // Verify triangles after
  int trisA[] = {0,  5,  9,  5,  0,  1,  1,  2,  6,  14, 6,  11, 1,  6,  5,  5,  6,  10, 9,
                 13, 17, 17, 13, 18, 14, 10, 6,  18, 14, 19, 14, 18, 13, 9,  5,  13, 6,  2,
                 7,  3,  4,  8,  3,  8,  7,  7,  2,  3,  12, 7,  8,  19, 15, 20, 15, 19, 11,
                 12, 15, 7,  15, 12, 16, 20, 16, 21, 16, 20, 15, 12, 8,  16, 11, 19, 14};
  VecInt trisAfter(&trisA[0], &trisA[75]);
  TS_ASSERT_EQUALS(trisAfter, triangles);

} // MePolyMesherUnitTests::test1

  //} // namespace xms

#endif // CXX_TEST
