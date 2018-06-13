//------------------------------------------------------------------------------
/// \file
/// \ingroup meshing_detail
/// \copyright (C) Copyright Aquaveo 2018. Distributed under the xmsng
///  Software License, Version 1.0. (See accompanying file
///  LICENSE_1_0.txt or copy at http://www.aquaveo.com/xmsng/LICENSE_1_0.txt)
//------------------------------------------------------------------------------

//----- Included files ---------------------------------------------------------

// 1. Precompiled header

// 2. My own header
#include <xmsmesh/meshing/detail/MeRefinePtsToPolys.h>

// 3. Standard library headers
#include <map>
#include <sstream>

// 4. External library headers

// 5. Shared code headers
#include <xmscore/math/math.h>
#include <xmscore/points/pt.h>
#include <xmsinterp/geometry/GmPolygon.h>
#include <xmsinterp/geometry/GmPtSearch.h>
#include <xmsmesh/meshing/MeMultiPolyMesherIo.h>
#include <xmscore/misc/XmLog.h>

// 6. Non-shared code headers

//----- Forward declarations ---------------------------------------------------

//----- External globals -------------------------------------------------------

//----- Namespace declaration --------------------------------------------------
namespace xms
{
//----- Constants / Enumerations -----------------------------------------------

//----- Classes / Structs ------------------------------------------------------

//----- Internal functions -----------------------------------------------------

//----- Class / Function definitions -------------------------------------------
class MeRefinePtsToPolysImpl : public MeRefinePtsToPolys
{
public:
  MeRefinePtsToPolysImpl()
  : m_TwoTimesSqrtThree(3.4641016151377545870548926830117)
  , m_SqrtThree(1.7320508075688772935274463415059)
  , m_SqrtThreeOverTwo(0.86602540378443864676372317075294)
  , m_xyTol(1e-9)
  {
  }

  //------------------------------------------------------------------------------
  /// \brief Sets the refine points.
  /// \param[in] a_pts: The refine points.
  /// \param[in] a_tol: The xy tolerance.
  //------------------------------------------------------------------------------
  virtual void SetRefinePoints(const std::vector<MeRefinePoint>& a_pts, double a_tol) override
  {
    m_pts = a_pts;
    m_xyTol = a_tol;
  }

  virtual void RefPtsAsPolys(const std::vector<Pt3d>& a_outPoly,
                             const std::vector<std::vector<Pt3d>>& a_inPolys,
                             std::vector<std::vector<Pt3d>>& a_newInPolys,
                             std::vector<Pt3d>& a_refMeshPts,
                             std::vector<Pt3d>& a_refPtsTooClose) override;

  void FindPtsInsidePolygon(const std::vector<Pt3d>& a_outPoly,
                            const std::vector<std::vector<Pt3d>>& a_inPolys,
                            std::vector<Pt3d>& a_refPtsTooClose);
  void CheckRefPtsTooCloseToOtherRefPts(std::multimap<double, size_t>& a_mapSizeIdx,
                                        std::vector<Pt3d>& a_refPtsTooClose);
  void CreateNewInsidePolygons(std::vector<std::vector<Pt3d>>& a_newInPolys,
                               std::vector<Pt3d>& a_refMeshPts,
                               std::vector<Pt3d>& a_refPtsProcessed);
  std::vector<Pt3d> TriPolyAtPoint(const Pt3d& a_pt, double a_size);
  std::vector<Pt3d> HexPolyAtPoint(const Pt3d& a_pt, double a_size);

  double m_TwoTimesSqrtThree;          ///< precalculated constant
  double m_SqrtThree;                  ///< precalculated constant
  double m_SqrtThreeOverTwo;           ///< precalculated constant
  double m_xyTol;                      ///< tolerance used for geometric comparison
  std::vector<MeRefinePoint> m_pts;    ///< the refine points
  std::vector<size_t> m_ptsInsidePoly; ///< indexes to the refine point that are inside of the
                                       ///< current polygon being considered
};
//------------------------------------------------------------------------------
/// \brief Creates a new instance of this class
/// \return MeRefinePtsToPolys.
//------------------------------------------------------------------------------
BSHP<MeRefinePtsToPolys> MeRefinePtsToPolys::New()
{
  BSHP<MeRefinePtsToPolys> ret(new MeRefinePtsToPolysImpl);
  return ret;
} // MeRefinePtsToPolys::New
//------------------------------------------------------------------------------
/// \brief
//------------------------------------------------------------------------------
MeRefinePtsToPolys::MeRefinePtsToPolys()
{
}
//------------------------------------------------------------------------------
/// \brief
//------------------------------------------------------------------------------
MeRefinePtsToPolys::~MeRefinePtsToPolys()
{
}

////////////////////////////////////////////////////////////////////////////////
/// \class MeRefinePtsToPolysImpl
/// \brief Creates polygons from refine point input.
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief Creates new inside polygons from refine points that are inside of
/// the polygon being considered
/// \param[in] a_outPoly The points making up the outside polygon
/// \param[in] a_inPolys Inner boundaries for the a_outPoly
/// \param[out] a_newInPolys New inside polygons created from refine points
/// \param[out] a_refMeshPts Locations of refine points that are also mesh node
/// locations.
/// \param[out] a_refPtsProcessed Locations of other refine points that are not
/// also mesh nodes.
//------------------------------------------------------------------------------
void MeRefinePtsToPolysImpl::RefPtsAsPolys(const std::vector<Pt3d>& a_outPoly,
                                           const std::vector<std::vector<Pt3d>>& a_inPolys,
                                           std::vector<std::vector<Pt3d>>& a_newInPolys,
                                           std::vector<Pt3d>& a_refMeshPts,
                                           std::vector<Pt3d>& a_refPtsProcessed)
{
  m_ptsInsidePoly.resize(0);
  a_newInPolys.resize(0);
  a_refMeshPts.resize(0);
  a_refPtsProcessed.resize(0);
  FindPtsInsidePolygon(a_outPoly, a_inPolys, a_refPtsProcessed);
  CreateNewInsidePolygons(a_newInPolys, a_refMeshPts, a_refPtsProcessed);
} // MeRefinePtsToPolysImpl::RefPtsAsPolys
//------------------------------------------------------------------------------
/// \brief Finds the refine points that are inside of the polygon.
/// \param a_outPoly The points defining the outer loop of the polygon.
/// Clockwise. 1st pt != last pt.
/// \param a_inPolys The points defining inner loops inside the outer loop of
/// the polygon. Can be empty. Counter clockwise. 1st pt != last pt.
/// \param a_refPtsProcessed Vector of point locations that have been
/// process by the meshing algorithm.
//------------------------------------------------------------------------------
void MeRefinePtsToPolysImpl::FindPtsInsidePolygon(const std::vector<Pt3d>& a_outPoly,
                                                  const std::vector<std::vector<Pt3d>>& a_inPolys,
                                                  std::vector<Pt3d>& a_refPtsProcessed)
{
  BSHP<GmPolygon> gmPoly = GmPolygon::New();
  gmPoly->Setup(a_outPoly, a_inPolys);
  std::multimap<double, size_t> mapSizeIdx;
  for (size_t i = 0; i < m_pts.size(); ++i)
  {
    if (gmPoly->Within(m_pts[i].m_pt))
    {
      // see if this point is far enough away from the boundary of the
      // polygon. If not then push a message that this point was too close
      // to the polygon boundary.
      double dist = gmPoly->MinDistanceToBoundary(m_pts[i].m_pt);
      if (dist < m_pts[i].m_size)
      {
        std::stringstream ss, loc;
        loc << "(" << m_pts[i].m_pt.x << ", " << m_pts[i].m_pt.y << ")";
        ss << "Refine point at location: " << loc.str()
           << " is too close to the polygon boundary with specified size: " << m_pts[i].m_size
           << ". The point was not inserted by the meshing process. Specify a "
              "size smaller than "
           << dist << " for the point to be included by the meshing process.";
        XM_LOG(xmlog::error, ss.str());
        a_refPtsProcessed.push_back(m_pts[i].m_pt);
      }
      else
        mapSizeIdx.insert(std::make_pair(m_pts[i].m_size, i));
    }
  }
  // now check points that are too close to other points
  CheckRefPtsTooCloseToOtherRefPts(mapSizeIdx, a_refPtsProcessed);
} // MeRefinePtsToPolysImpl::FindPtsInsidePolygon
//------------------------------------------------------------------------------
/// \brief Checks on refine points that are inside of the polygon to make sure
/// that they are not too close to one another. Points are sorted based on the
/// refinement size so that preference is given to the smallest refinement
/// size.
/// \param a_mapSizeIdx A map containing the refinement size and the index to
/// the point.
/// \param a_refPtsProcessed Vector of point locations that have been
/// process by the meshing algorithm.
//------------------------------------------------------------------------------
void MeRefinePtsToPolysImpl::CheckRefPtsTooCloseToOtherRefPts(
  std::multimap<double, size_t>& a_mapSizeIdx,
  std::vector<Pt3d>& a_refPtsProcessed)
{
  std::vector<int> removePt(a_mapSizeIdx.size(), 0);
  std::multimap<double, size_t>::iterator it(a_mapSizeIdx.begin()), itEnd(a_mapSizeIdx.end());
  for (size_t i = 0; it != itEnd; ++i, ++it)
  {
    if (removePt[i])
      continue;

    size_t iIdx(it->second);
    double iSize(m_pts[iIdx].m_size);
    const Pt3d& pi(m_pts[iIdx].m_pt);
    m_ptsInsidePoly.push_back(iIdx);

    std::multimap<double, size_t>::iterator it2(it);
    ++it2;
    for (size_t j = i; it2 != itEnd; ++j, ++it2)
    {
      size_t jIdx(it2->second);
      const Pt3d& pj(m_pts[jIdx].m_pt);
      double jSize(m_pts[jIdx].m_size);
      double dist(Mdist(pi.x, pi.y, pj.x, pj.y));
      double totalSize(iSize + jSize);
      if (dist < totalSize)
      {
        double target = dist - iSize;
        removePt[j] = true;
        std::stringstream ss, loc;
        loc << "(" << pj.x << ", " << pj.y << ")";
        ss << "Refine point at location: " << loc.str() << " with specified size: " << jSize
           << " is too close to another refine point. The point was not "
              "inserted by the meshing process. Specify a size smaller than "
           << target << " for the point to be included by the meshing process.";
        XM_LOG(xmlog::error, ss.str());
        a_refPtsProcessed.push_back(pj);
      }
    }
  }
} // MeRefinePtsToPolysImpl::CheckRefPtsTooCloseToOtherRefPts
//------------------------------------------------------------------------------
/// \brief Creates new inside polygons from the refine points and appends
/// these polygons to the current vector of inside polys.
/// \param a_newInPolys New inside polygons with refine points included.
/// \param a_refMeshPts Vector of refine point locations that are included
/// in the final mesh by the meshing algorithm.
/// \param a_refPtsProcessed Vector of point locations that have been
/// processed by the meshing algorithm but are not represented by mesh nodes
/// in the final mesh.
//------------------------------------------------------------------------------
void MeRefinePtsToPolysImpl::CreateNewInsidePolygons(std::vector<std::vector<Pt3d>>& a_newInPolys,
                                                     std::vector<Pt3d>& a_refMeshPts,
                                                     std::vector<Pt3d>& a_refPtsProcessed)
{
  for (size_t i = 0; i < m_ptsInsidePoly.size(); ++i)
  {
    size_t idx = m_ptsInsidePoly[i];
    std::vector<Pt3d> p;
    if (m_pts[idx].m_size > m_xyTol)
    {
      if (m_pts[idx].m_createMeshPoint)
      {
        p = HexPolyAtPoint(m_pts[idx].m_pt, m_pts[idx].m_size);
        a_refMeshPts.push_back(m_pts[idx].m_pt);
      }
      else
      {
        p = TriPolyAtPoint(m_pts[idx].m_pt, m_pts[idx].m_size);
        a_refPtsProcessed.push_back(m_pts[idx].m_pt);
      }
      a_newInPolys.push_back(p);
    }
    else if (m_pts[idx].m_createMeshPoint)
      a_refMeshPts.push_back(m_pts[idx].m_pt);
  }
} // MeRefinePtsToPolysImpl::CreateNewInsidePolygons
//------------------------------------------------------------------------------
/// \brief Creates a triangle polygon that surrounds a refine point. The
/// refine point itself will _NOT_ be included in the mesh.
/// \param a_pt The location of the refine point.
/// \param a_size The element edge size a the refine point.
/// \return A polygon for a_pt refine point.
//------------------------------------------------------------------------------
std::vector<Pt3d> MeRefinePtsToPolysImpl::TriPolyAtPoint(const Pt3d& a_pt, double a_size)
{
  // Equilateral triangle with edge length "l"
  //
  //                                      * (x, y + l / sqrt(3))
  //                                     / \
//                                    /   \
//                                   /     \
//                                  /       \
//                                 /         \
//                                /           \
//                               /             \
//                              /               \
//                             /                 \
//                            /         *         \
//                           /        (x,y)        \
//                          /                       \
//                         /                         \
//                        /                           \
//(x-l/2, y-l/2*sqrt(3)) *-----------------------------* (x+l/2, y-l/2*sqrt(3))
  //
  //                       |------------ l --------------|

  double sizeOver2(a_size / 2);
  double sizeOver2TimesSqrt3(a_size / m_TwoTimesSqrtThree);
  double sizeOverSqrt3(a_size / m_SqrtThree);
  std::vector<Pt3d> rval(3, a_pt);
  rval[0].x -= sizeOver2;
  rval[0].y -= sizeOver2TimesSqrt3;
  rval[1].x += sizeOver2;
  rval[1].y -= sizeOver2TimesSqrt3;
  rval[2].y += sizeOverSqrt3;
  return rval;
} // MeRefinePtsToPolysImpl::TriPolyAtPoint
//------------------------------------------------------------------------------
/// \brief Creates a hexagon polygon that surrounds a refine point. The
/// refine point itself _WILL_ be included in the mesh.
/// \param a_pt The location of the refine point.
/// \param a_size The element edge size a the refine point.
/// \return A polygon for a_pt refine point.
//------------------------------------------------------------------------------
std::vector<Pt3d> MeRefinePtsToPolysImpl::HexPolyAtPoint(const Pt3d& a_pt, double a_size)
{
  // Hexagon with element edge length equal to "l"
  // the center is at (x, y)
  //
  //
  //   (x-l/2,y+l*sqrt(3)/2)
  //                     *-----------*(x+l/2,y+l*sqrt(3)/2)
  //                    / \         / \
//                   /   \       /   \
//                  /     \     /     \
//                 /       \   /       \
//                /         \ /         \
//      (x-l/2,y)*-----------*-----------* (x+l/2,y)
  //                \         / \         /
  //                 \       /   \       /
  //                  \     /     \     /
  //                   \   /       \   /
  //                    \ /         \ /
  //                     *-----------*(x+l/2,y-l*sqrt(3)/2)
  //   (x-l/2,y-l*sqrt(3)/2)
  //
  // squence to define inner poly that will have the inside pt
  // 0,1,2,0,2,3,0,3,4,0,4,5,0,5,6,0,6,1
  //
  double sizeOver2(a_size / 2);
  double sizeTimesSqrtThreeOverTwo(a_size * m_SqrtThreeOverTwo);
  std::vector<Pt3d> pts(6, a_pt);
  pts[0].x += a_size;
  pts[1].x += sizeOver2;
  pts[1].y += sizeTimesSqrtThreeOverTwo;
  pts[2].x -= sizeOver2;
  pts[2].y += sizeTimesSqrtThreeOverTwo;
  pts[3].x -= a_size;
  pts[4].x -= sizeOver2;
  pts[4].y -= sizeTimesSqrtThreeOverTwo;
  pts[5].x += sizeOver2;
  pts[5].y -= sizeTimesSqrtThreeOverTwo;
  return pts;
  // int sequence[] = {0,1,2,0,2,3,0,3,4,0,4,5,0,5,6,0,6,1};
  // std::vector<Pt3d> rval(18, a_pt);
  // for (size_t i=0; i<rval.size(); ++i) rval[i] = pts[sequence[i]];
  // return rval;
} // MeRefinePtsToPolysImpl::HexPolyAtPoint

} // namespace xms

#ifdef CXX_TEST
////////////////////////////////////////////////////////////////////////////////

#include <xmsmesh/meshing/detail/MeRefinePtsToPolys.t.h>

#include <xmscore/testing/TestTools.h>
#include <xmsinterp/geometry/geoms.h>

// namespace xms {
using namespace xms;

////////////////////////////////////////////////////////////////////////////////
/// \class MeRefinePtsToPolysUnitTests
/// \brief tester for the MeRefinePtsToPolys class
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief tests creating the class
//------------------------------------------------------------------------------
void MeRefinePtsToPolysUnitTests::testCreateClass()
{
  BSHP<MeRefinePtsToPolys> b = MeRefinePtsToPolys::New();
  TS_ASSERT(b);
} // MeRefinePtsToPolysUnitTests::testCreateClass
//------------------------------------------------------------------------------
/// \brief tests creating a triangle polygon at the refined point
//------------------------------------------------------------------------------
void MeRefinePtsToPolysUnitTests::testTriPolyAtPoint()
{
  MeRefinePtsToPolysImpl p;
  Pt3d pt(1, 1, 0);
  std::vector<Pt3d> pts = p.TriPolyAtPoint(pt, 1);
  std::vector<Pt3d> basePts = {{.5, .7113, 0}, {1.5, .7113, 0}, {1, 1.5773, 0}};
  TS_ASSERT_DELTA_VECPT3D(basePts, pts, 1e-3);
} // MeRefinePtsToPolysUnitTests::testTriPolyAtPoint
//------------------------------------------------------------------------------
/// \brief test creating a hex polygon at the refine point
//------------------------------------------------------------------------------
void MeRefinePtsToPolysUnitTests::testHexPolyAtPoint()
{
  MeRefinePtsToPolysImpl p;
  Pt3d pt(1, 1, 0);
  std::vector<Pt3d> pts = p.HexPolyAtPoint(pt, 1);
  std::vector<Pt3d> hexPts = {{2, 1, 0}, {1.5, 1.866, 0}, {.5, 1.866, 0},
                              {0, 1, 0}, {.5, .1339, 0},  {1.5, .1339, 0}};
  TS_ASSERT_DELTA_VECPT3D(hexPts, pts, 1e-3);
  // int seq[] = {0,1,2,0,2,3,0,3,4,0,4,5,0,5,6,0,6,1};
  // std::vector<Pt3d> basePts(18);
  // for (size_t i=0; i<basePts.size(); ++i) basePts[i] = hexPts[seq[i]];
  // TS_ASSERT_DELTA_VECPT3D(basePts, pts, 1e-3);
} // MeRefinePtsToPolysUnitTests::testHexPolyAtPoint
//------------------------------------------------------------------------------
/// \brief tests multiple refine points to polygons
//------------------------------------------------------------------------------
void MeRefinePtsToPolysUnitTests::testRefPtsAsPolys()
{
  MeRefinePtsToPolysImpl p;
  double tol(1e-9);
  std::vector<MeRefinePoint> refPts = {{{4, 4, 0}, 1, 1}, {{8, 8, 0}, 1, 0}};
  p.SetRefinePoints(refPts, tol);
  std::vector<Pt3d> outPoly = {{0, 0, 0}, {0, 10, 0}, {10, 10, 0}, {10, 0, 0}};
  std::vector<Pt3d> refMeshPts, tooClose;
  std::vector<std::vector<Pt3d>> inPolys, newInPolys;
  p.RefPtsAsPolys(outPoly, inPolys, newInPolys, refMeshPts, tooClose);
  TS_ASSERT_EQUALS(1, refMeshPts.size());
  TS_ASSERT_EQUALS(2, newInPolys.size());
  if (2 != newInPolys.size())
    return;

  std::vector<Pt3d> hexPts = {{5, 4, 0}, {4.5, 4.866, 0},  {3.5, 4.866, 0},
                              {3, 4, 0}, {3.5, 3.1339, 0}, {4.5, 3.1339, 0}};
  TS_ASSERT_DELTA_VECPT3D(hexPts, newInPolys[0], 1e-3);
  // int seq[] = {0,1,2,0,2,3,0,3,4,0,4,5,0,5,6,0,6,1};
  std::vector<Pt3d> basePts(18);
  // for (size_t i=0; i<basePts.size(); ++i) basePts[i] = hexPts[seq[i]];
  // TS_ASSERT_DELTA_VECPT3D(basePts, newInPolys[0], 1e-3);

  basePts = {{7.5, 7.71132, 0}, {8.5, 7.71132, 0}, {8, 8.57735, 0}};
  TS_ASSERT_DELTA_VECPT3D(basePts, newInPolys[1], 1e-3);
} // MeRefinePtsToPolysUnitTests::testRefPtsAsPolys
//------------------------------------------------------------------------------
/// \brief test refine points that are too close to the polygon boundary
//------------------------------------------------------------------------------
void MeRefinePtsToPolysUnitTests::testRefinePtsTooCloseToBoundary()
{
  XmLog::Instance().GetAndClearStackStr();
  MeRefinePtsToPolysImpl p;
  double tol(1e-9);
  std::vector<MeRefinePoint> refPts = {{{.5, 4, 0}, 1, 1}, {{9.001, 8, 0}, 1, 0}};
  p.SetRefinePoints(refPts, tol);
  std::vector<Pt3d> outPoly = {{0, 0, 0}, {0, 10, 0}, {10, 10, 0}, {10, 0, 0}};
  std::vector<Pt3d> refMeshPts, tooClose;
  std::vector<std::vector<Pt3d>> inPolys, newInPolys;
  p.RefPtsAsPolys(outPoly, inPolys, newInPolys, refMeshPts, tooClose);
  TS_ASSERT_EQUALS(0, newInPolys.size());
  TS_ASSERT_EQUALS(2, XmLog::Instance().ErrCount());
  TS_ASSERT_STACKED_ERRORS(
    "---Refine point at location: (0.5, 4) is too close to the polygon boundary "
    "with specified size: 1. The point was not inserted by the meshing "
    "process. Specify a size smaller than 0.5 for the point to be included "
    "by the meshing process.\n\n"
    "---Refine point at location: (9.001, 8) is too close to the polygon "
    "boundary with specified size: 1. The point was not inserted by the "
    "meshing process. Specify a size smaller than 0.999 for the point to be "
    "included by the meshing process.\n\n");
} // MeRefinePtsToPolysUnitTests::testRefinePtsTooCloseToBoundary

//} // namespace xms
#endif // CXX_TEST