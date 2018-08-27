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
#include <xmsmesh/meshing/detail/MePolyPts.h>

// 3. Standard library headers

// 4. External library headers
#pragma warning(push)
#pragma warning(disable : 4512) // boost code: no assignment operator
#include <boost/geometry/geometry.hpp>
#pragma warning(pop)

// 5. Shared code headers
#include <xmscore/math/math.h>
#include <xmsinterp/geometry/GmBoostTypes.h> // GmBstPoly3d, XmBstRing
#include <xmsinterp/geometry/geoms.h>
#include <xmsinterp/geometry/GmPtSearch.h>
#include <xmsmesh/meshing/detail/MePolyOffsetter.h>
#include <xmscore/misc/boost_defines.h>
#include <xmscore/misc/XmError.h>

// 6. Non-shared code headers

//----- Forward declarations ---------------------------------------------------

//----- External globals -------------------------------------------------------

//----- Namespace declaration --------------------------------------------------
namespace xms
{
/// map for segment intersection indices
typedef std::map<size_t, std::multimap<double, size_t>> SegMap;
/// map used to see where segments cross each other
typedef std::map<const std::vector<size_t>*, SegMap> MapSegMap;

//----- Constants / Enumerations -----------------------------------------------

//----- Classes / Structs ------------------------------------------------------
namespace
{
namespace bg = boost::geometry;

#define T_TOL 1e-13 ///< tolerance used in multipoly intersector

} // unnamed namespace
/// \brief Implementation of MePolyPts
class MePolyPts::impl
{
public:
  impl()
  : m_xyTol(1e-9)
  , m_pts(new std::vector<Pt3d>())
  {
  }
  ~impl() {}

  double m_xyTol;                ///< tolerance for geometric comparisons
  BSHP<std::vector<Pt3d>> m_pts; ///< point locations
  MapSegMap m_segCross;          ///< map used to see where segment cross each other
  BSHP<GmPtSearch> m_ps;         ///< spatial index for searching points
};

//----- Internal functions -----------------------------------------------------

//----- Class / Function definitions -------------------------------------------

//------------------------------------------------------------------------------
/// \brief Get the next index for a segment. The segments are stored as a
/// vector and it forms a closed loop so the last item in the vector makes
/// a segment with the first item in the vector.
/// \param a_i Index into the vector a_seg.
/// \param a_seg Vector of indexes that make up a closed loop polyline.
/// \return index to the next segment index
//------------------------------------------------------------------------------
static size_t iNextSegIdx(size_t a_i, const std::vector<size_t>& a_seg)
{
  size_t next(a_i + 1);
  if (next >= a_seg.size())
    next = 0;
  return next;
} // iNextSegIdx
//------------------------------------------------------------------------------
/// \brief Adds an intersection to a map of intersection information.
/// \param a_i Index for a segment.
/// \param a_iSeg Vector of segments that a_i refers to.
/// \param a_j Index for a point location.
/// \param a_t Parametric value for the location of point a_j along segment a_i.
/// \param a_mapSegCross Map that holds information about intersections of line
/// segments.
//------------------------------------------------------------------------------
static void iAddIntersection(size_t a_i,
                             const std::vector<size_t>& a_iSeg,
                             size_t a_j,
                             double a_t,
                             MapSegMap& a_mapSegCross)
{
  if (a_t >= 0 && a_t <= 1)
  {
    a_mapSegCross[&a_iSeg][a_i].insert(std::make_pair(a_t, a_j));
  }
} // iAddIntersection
//------------------------------------------------------------------------------
/// \brief hash the points and fill in a vector with the new index assigned to
/// each point. If there are no duplicates the a_newIdx will be 0,1,2,3,4
/// If there are duplicates then a_newIdx will be 0,1,0,1,4,5,0...
/// \param a_newIdx Vector of new indexes for points.
/// \param a_pts Vector of point locations.
/// \param a_xyTol Tolerance used to determine if 2 locations should be
/// considered the same location.
//------------------------------------------------------------------------------
static void iHashPts(std::vector<size_t>& a_newIdx, BSHP<std::vector<Pt3d>> a_pts, double a_xyTol)
{
  a_newIdx.resize(0);
  if (a_pts->empty())
    return;

  BSHP<GmPtSearch> ps = GmPtSearch::New(true);
  ps->PtsToSearch(a_pts);

  a_newIdx.resize(a_pts->size(), 0);
  for (size_t i = 0; i < a_newIdx.size(); ++i)
    a_newIdx[i] = i;

  std::vector<int> n;
  for (size_t i = 0; i < a_pts->size(); ++i)
  {
    ps->PtsWithinDistanceToPtInRtree((int)i, (*a_pts)[i], a_xyTol, n);
    for (size_t j = 0; j < n.size(); ++j)
    {
      if ((size_t)n[j] > i)
      {
        a_newIdx[n[j]] = i;
      }
    }
  }
} // iHashPts
//------------------------------------------------------------------------------
/// \brief Returns the envelopes of the segments
/// \param a_segs Vector of indexes defining the segments. This is a closed
/// loop polyline the first point connects to the last point in the vector.
/// \param a_pts Vector of locations referenced by the indexes in a_segs.
/// \return vector of envelopes, one for each segment
//------------------------------------------------------------------------------
static std::vector<GmBstBox3d> iCalcSegEnvelopes(const std::vector<size_t>& a_segs,
                                                 const std::vector<Pt3d>& a_pts)
{
  GmBstBox3d b;
  std::vector<GmBstBox3d> envelopes;
  envelopes.reserve(a_segs.size());
  for (size_t i = 0; i < a_segs.size(); ++i)
  {
    size_t ix0(a_segs[i]);
    size_t next(iNextSegIdx(i, a_segs));
    size_t ix1(a_segs[next]);

    const Pt3d &p0(a_pts[ix0]), &p1(a_pts[ix1]);
    b.max_corner() = b.min_corner() = p0;
    b.max_corner().z = b.min_corner().z = 0;
    if (b.min_corner().x > p1.x)
      b.min_corner().x = p1.x;
    if (b.max_corner().x < p1.x)
      b.max_corner().x = p1.x;

    if (b.min_corner().y > p1.y)
      b.min_corner().y = p1.y;
    if (b.max_corner().y < p1.y)
      b.max_corner().y = p1.y;

    envelopes.push_back(b);
  }
  return envelopes;
} // iCalcSegEnvelopes
//------------------------------------------------------------------------------
/// \brief Returns true if the envelopes overlap or touch
/// \param a_i Index of an envelope of a segment.
/// \param a_j Index of an envelope of a segment.
/// \param a_iEnv Vector of envelopes referenced by a_i.
/// \param a_jEnv Vector of envelopes referenced by a_j.
/// \return true if the envelopes overlap or touch
//------------------------------------------------------------------------------
static bool iEnvelopesOverlapOrTouch(size_t a_i,
                                     size_t a_j,
                                     const std::vector<GmBstBox3d>& a_iEnv,
                                     const std::vector<GmBstBox3d>& a_jEnv)
{
  const Pt3d &iMin(a_iEnv[a_i].min_corner()), &iMax(a_iEnv[a_i].max_corner()),
    &jMin(a_jEnv[a_j].min_corner()), &jMax(a_jEnv[a_j].max_corner());
  if (jMin.x > iMax.x || jMin.y > iMax.y || iMin.x > jMax.x || iMin.y > jMax.y)
    return false;
  return true;
} // iEnvelopesOverlapOrTouch
//------------------------------------------------------------------------------
/// \brief Computes the T value for the point on the line a_p0, a_p1
/// \param a_p0 Location of the 1st point defining a line segment.
/// \param a_p1 Location of the 2nd point defining a line segment.
/// \param a_iPt Location of the intersection that lies on the line segment
/// defined by a_p0, a_p1.
/// \param a_xyTol Tolerance used in geometric calculations.
/// \return Parametric value for location of intersection.
//------------------------------------------------------------------------------
static double iCalcIntersectionT(const Pt3d& a_p0, const Pt3d& a_p1, Pt3d& a_iPt, double a_xyTol)
{
  double tval(-1);
  double denomX = a_p1.x - a_p0.x;
  double denomY = a_p1.y - a_p0.y;
  if (fabs(denomX) > fabs(denomY))
  {
    tval = (a_iPt.x - a_p0.x) / denomX;
  }
  else
  {
    tval = (a_iPt.y - a_p0.y) / denomY;
  }
  if (gmEqualPointsXY(a_p0, a_iPt, a_xyTol) || (tval < 0 && EQ_TOL(tval, 0, T_TOL)))
  {
    a_iPt = a_p0;
    tval = 0;
  }
  if (gmEqualPointsXY(a_p1, a_iPt, a_xyTol) || (tval > 1 && EQ_TOL(tval, 1, T_TOL)))
  {
    a_iPt = a_p1;
    tval = 1;
  }
  return tval;
} // iCalcIntersectionT
//------------------------------------------------------------------------------
/// \brief Check segments to see if they share a point index
/// \param a_i Index to a segment.
/// \param a_j Index to a segment.
/// \param a_iSeg Vector of segments referred to by a_i.
/// \param a_jSeg Vector of segments referred to by a_j.
/// \param a_pts Vector of point locations. The segments just store
/// indexes that refer to locations in this vector.
/// \param a_xyTol Tolerance used in geometric calculations.
/// \param a_mapSegCross Map of intersection information.
/// \return true if the segments share a point
//------------------------------------------------------------------------------
static bool iCheckSharedPtIdx(size_t a_i,
                              size_t a_j,
                              const std::vector<size_t>& a_iSeg,
                              const std::vector<size_t>& a_jSeg,
                              const std::vector<Pt3d>& a_pts,
                              double a_xyTol,
                              MapSegMap& a_mapSegCross)
{
  size_t nexti(iNextSegIdx(a_i, a_iSeg)), nextj(iNextSegIdx(a_j, a_jSeg));
  size_t p00(a_iSeg[a_i]), p01(a_iSeg[nexti]), p10(a_jSeg[a_j]), p11(a_jSeg[nextj]);
  if (p00 != p10 && p00 != p11 && p01 != p10 && p01 != p11)
    return false;
  // check if this is the same segment (can happen after hashing points)
  if ((p00 == p10 || p00 == p11) && (p01 == p10 || p01 == p11))
    return true;

  const Pt3d &s1p0(a_pts[p00]), &s1p1(a_pts[p01]);
  size_t idx(0);
  if (p00 == p10 || p01 == p10)
    idx = p11;
  else
    idx = p10;

  Pt3d p = a_pts[idx];
  if (gmOnLineAndBetweenEndpointsWithTol(s1p0, s1p1, p.x, p.y, a_xyTol))
  {
    double t = iCalcIntersectionT(s1p0, s1p1, p, a_xyTol);
    iAddIntersection(a_i, a_iSeg, idx, t, a_mapSegCross);
  }
  return true;
} // iCheckSharedPtIdx
//------------------------------------------------------------------------------
/// \brief Intersects 2 segments. Normally there is only 1 intersection between
/// 2 segments but if the lines are colinear and they overlap then we get the
/// endpoint of the overlap as the intersection points.
/// \param a_s1p0 Location of 1st point on segment 1.
/// \param a_s1p1 Location of 2nd point on segment 1.
/// \param a_s2p0 Location of 1st point on segment 2.
/// \param a_s2p1 Location of 2nd point on segment 2.
/// \param a_xyTol Tolerance used in geometric calculations.
/// \param a_pts Locations of intersections between 2 segments.
/// \return true if the segments intersect
//------------------------------------------------------------------------------
static bool iIntersectTwoSegments(const Pt3d& a_s1p0,
                                  const Pt3d& a_s1p1,
                                  const Pt3d& a_s2p0,
                                  const Pt3d& a_s2p1,
                                  double a_xyTol,
                                  std::vector<Pt3d>& a_pts)
{
  bool rval(false);
  a_pts.resize(0);

  // check parallel
  double dx1 = a_s1p1.x - a_s1p0.x;
  double dy1 = a_s1p1.y - a_s1p0.y;
  double dx2 = a_s2p1.x - a_s2p0.x;
  double dy2 = a_s2p1.y - a_s2p0.y;
  double cross = (dx1 * dy2) - (dy1 * dx2);
  if (EQ_TOL(0.0, cross, T_TOL))
  { // if we are parallel then we just want to
    // check if p1 and seg2 is on seg1
    if (gmOnLineAndBetweenEndpointsWithTol(a_s1p0, a_s1p1, a_s2p1.x, a_s2p1.y, a_xyTol))
    {
      a_pts.push_back(a_s2p1);
      rval = true;
    }
    if (gmOnLineAndBetweenEndpointsWithTol(a_s2p0, a_s2p1, a_s1p1.x, a_s1p1.y, a_xyTol))
    {
      a_pts.push_back(a_s1p1);
      rval = true;
    }
  }
  else
  {
    double z2;
    Pt3d pt;
    if (gmIntersectLineSegmentsWithTol(a_s1p0, a_s1p1, a_s2p0, a_s2p1, &pt.x, &pt.y, &pt.z, &z2,
                                       a_xyTol))
    {
      a_pts.push_back(pt);
      rval = true;
    }
  }
  return rval;
} // iIntersectTwoSegments
//------------------------------------------------------------------------------
/// \brief Sees if 2 segments intersects and if they do then information is
/// added to the a_mapSegCross.
/// \param a_i Index to a segment.
/// \param a_j Index to a segment.
/// \param a_iSeg Vector of segments referred to by a_i.
/// \param a_jSeg Vector of segments referred to by a_j.
/// \param a_shptrPts Vector of point locations. The segments just store
/// indexes that refer to locations in this vector.
/// \param a_xyTol Tolerance used in geometric calculations.
/// \param a_mapSegCross Map of intersection information.
/// \param a_ps A PtSearch class. Will hash newly calculated locations with
/// existing point locations.
//------------------------------------------------------------------------------
static void iCheckIntersection(size_t a_i,
                               size_t a_j,
                               const std::vector<size_t>& a_iSeg,
                               const std::vector<size_t>& a_jSeg,
                               BSHP<std::vector<Pt3d>>& a_shptrPts,
                               double a_xyTol,
                               MapSegMap& a_mapSegCross,
                               BSHP<GmPtSearch> a_ps)
{
  std::vector<Pt3d>& spts(*a_shptrPts);
  if (iCheckSharedPtIdx(a_i, a_j, a_iSeg, a_jSeg, spts, a_xyTol, a_mapSegCross))
    return;

  // first segment
  size_t nexti = iNextSegIdx(a_i, a_iSeg);
  Pt3d &s1p0(spts[a_iSeg[a_i]]), &s1p1(spts[a_iSeg[nexti]]);
  // second segment
  size_t nextj = iNextSegIdx(a_j, a_jSeg);
  Pt3d &s2p0(spts[a_jSeg[a_j]]), &s2p1(spts[a_jSeg[nextj]]);

  std::vector<Pt3d> pts;
  bool intersect = iIntersectTwoSegments(s1p0, s1p1, s2p0, s2p1, a_xyTol, pts);
  if (intersect)
  {
    for (size_t i = 0; i < pts.size(); ++i)
    {
      // calculate the "t" value for the intersection on the first segment
      double t = iCalcIntersectionT(s1p0, s1p1, pts[i], a_xyTol);
      double t2 = iCalcIntersectionT(s2p0, s2p1, pts[i], a_xyTol);
      // if t = 1 then we want this intersection to go on the second segment
      // with the t value on the second segment so that intersections on that
      /// segment are sorted correctly
      if (1 == t)
      {
        t = iCalcIntersectionT(s2p0, s2p1, pts[i], a_xyTol);
        iAddIntersection(a_j, a_jSeg, a_iSeg[nexti], t, a_mapSegCross);
      }
      // > 0 is here because at some point in our looping we will consider
      // the same location for both t = 0 and 1 if an endpoint of a segment
      // touches another segment. To keep from adding it twice we only add
      // when t > 0. See testCase2 as an example.
      else if (t > 0)
      {
        if (1 == t2 || 0 == t2)
        {
          if (1 == t2)
            iAddIntersection(a_i, a_iSeg, a_jSeg[nextj], t, a_mapSegCross);
          // ignore 0 because we will hit this point when t2 is 1
        }
        else
        {
          size_t ptIdx;
          if (!a_ps)
          {
            ptIdx = (int)spts.size();
            spts.push_back(pts[i]);
          }
          else
          {
            int idx;
            a_ps->AddPtToVectorIfUnique(pts[i], a_xyTol, idx);
            ptIdx = (size_t)idx;
          }
          iAddIntersection(a_i, a_iSeg, ptIdx, t, a_mapSegCross);
          iAddIntersection(a_j, a_jSeg, ptIdx, t2, a_mapSegCross);
        }
      }
      if (t < 0 || t > 1)
      {
        XM_ASSERT(0);
      }
    }
  }
} // iCheckIntersection
//------------------------------------------------------------------------------
/// \brief Returns true if a polygon is inside of another polygon
/// \param a_poly A polygon used to test a_polyToTest
/// \param a_polyToTest Vector of indices defining a polygon. Want to determine
/// if it is inside of a_poly.
/// \param a_pts Vector of locations referenced by a_polyToTest.
/// \return true if a_polyToTest is inside of a_poly
//------------------------------------------------------------------------------
static bool iPolyInsideOfPoly(const GmBstPoly3d& a_poly,
                              const std::vector<size_t>& a_polyToTest,
                              const std::vector<Pt3d>& a_pts)
{
  bool out = false;
  for (size_t i = 0; !out && i < a_polyToTest.size(); ++i)
  {
    const Pt3d& p(a_pts[a_polyToTest[i]]);
    if (!bg::covered_by(p, a_poly))
      out = true;
  }
  if (out)
    return false;
  return true;
} // iPolyInsideOfPoly
//------------------------------------------------------------------------------
/// \brief Given a list of polygons, finds the polygons that are inside of other
/// polygons.
/// \param a_polyInsideOfPoly 2D vector of polygon indexes.
/// \param a_loops 2D vector of point indexes that define polygons.
/// \param a_polyPts Class with information on locations of points that define
/// the polygons.
//------------------------------------------------------------------------------
static void iFindPolysInsideOfOtherPolys(std::vector<std::vector<size_t>>& a_polyInsideOfPoly,
                                         std::list<std::vector<size_t>>& a_loops,
                                         MePolyPts& a_polyPts)
{
  std::vector<Pt3d>& pts(a_polyPts.Pts());
  // create set of point indices that make up polygons
  std::list<std::vector<size_t>>::iterator it, it2;

  // determine which polygons are inside of other polygons
  size_t cnt_it(0);
  a_polyInsideOfPoly.resize(a_loops.size());
  for (it = a_loops.begin(); it != a_loops.end(); ++it, ++cnt_it)
  {
    // create polygon for iterator it
    GmBstPoly3d poly;
    for (size_t i = 0; i < it->size(); ++i)
    {
      bg::exterior_ring(poly).push_back(pts[(*it)[i]]);
    }
    bg::exterior_ring(poly).push_back(pts[(*it)[0]]);

    size_t cnt_it2(0);
    for (it2 = a_loops.begin(); it2 != a_loops.end(); ++it2, ++cnt_it2)
    {
      if (it == it2)
        continue;
      if (iPolyInsideOfPoly(poly, *it2, pts))
      {
        a_polyInsideOfPoly[cnt_it2].push_back(cnt_it);
      }
    }
  }
} // iFindPolysInsideOfOtherPolys
//------------------------------------------------------------------------------
/// \brief Calculates polygon areas
/// \param a_loops 2D vector of point indexes that define polygons.
/// \param a_pts Locations referred to by indexes in a_loops.
/// \return vector of areas, one for each polygon
//------------------------------------------------------------------------------
static std::vector<double> iCalcPolyAreas(const std::list<std::vector<size_t>>& a_loops,
                                          const std::vector<Pt3d>& a_pts)
{
  std::vector<double> areas;
  areas.reserve(a_loops.size());
  std::list<std::vector<size_t>>::const_iterator it(a_loops.begin());
  for (; it != a_loops.end(); ++it)
  {
    std::vector<Pt3d> pts(it->size(), Pt3d());
    for (size_t i = 0; i < it->size(); ++i)
      pts[i] = a_pts[(*it)[i]];
    areas.push_back(gmPolygonArea(&pts[0], (int)pts.size()));
  }
  return areas;
} // iCalcPolyAreas
//------------------------------------------------------------------------------
/// \brief Removes a repeated segment (2 numbers) from a sequence
/// \param a_sequence List of indexes defining a closed loop polyline.
//------------------------------------------------------------------------------
static void iRemoveRepeatedSegmentFromSequence(std::list<size_t>& a_sequence)
{
  if (a_sequence.size() < 3)
    return;

  std::list<size_t>::iterator it(a_sequence.begin());
  std::list<size_t>::iterator back1(it);
  ++it;
  std::list<size_t>::iterator back2(back1);
  ++back1;
  ++it;
  std::list<size_t>::iterator itEnd(a_sequence.end());
  while (it != itEnd && back1 != itEnd && back2 != itEnd)
  {
    if (*it == *back2)
    { // remove it and back1
      ++it;
      a_sequence.erase(back1, it);
      return;
    }
    else
    {
      ++back2;
      ++back1;
      ++it;
    }
  }
} // iRemoveRepeatedSegmentFromSequence

////////////////////////////////////////////////////////////////////////////////
/// \class MePolyPts
/// \brief Helper class used by PolyCleaner
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief constructor
//------------------------------------------------------------------------------
MePolyPts::MePolyPts()
: m_p(new MePolyPts::impl())
{
} // MpPolyPts::MpPolyPts
//------------------------------------------------------------------------------
/// \brief Returns the tolerance
/// \return the xy tolerance.
//------------------------------------------------------------------------------
double& MePolyPts::XyTol()
{
  return m_p->m_xyTol;
} // MpPolyPts::XyTol
//------------------------------------------------------------------------------
/// \brief Returns the vector of points.
/// \return The vector of points.
//------------------------------------------------------------------------------
std::vector<Pt3d>& MePolyPts::Pts()
{
  return *(m_p->m_pts);
} // MpPolyPts::Pts
//------------------------------------------------------------------------------
/// \brief Returns the vector of points.
/// \return Shared pointer to vector of points.
//------------------------------------------------------------------------------
BSHP<std::vector<Pt3d>> MePolyPts::PtsSharedPointer()
{
  return m_p->m_pts;
} // MpPolyPts::PtsSharedPointer
//------------------------------------------------------------------------------
/// \brief Hashes the point locations.
/// \return Vector of new indexes for points.
//------------------------------------------------------------------------------
std::vector<size_t> MePolyPts::HashPts()
{
  std::vector<size_t> idx;
  iHashPts(idx, PtsSharedPointer(), XyTol());
  return idx;
} // MePolyPts::HashPts
//------------------------------------------------------------------------------
/// \brief Returns the segments that will be used in CleanPolyOffset.
/// \return Segments that will be used in CleanPolyOffset.
//------------------------------------------------------------------------------
std::vector<size_t> MePolyPts::SegmentsForCleanPolyOffset()
{
  std::vector<size_t> tmpSegs(HashPts()), retSegs;
  // remove any consecutive indices in a_segs
  retSegs.reserve(Pts().size());
  if (!tmpSegs.empty())
    retSegs.push_back(tmpSegs[0]);
  for (size_t i = 1; i < tmpSegs.size(); ++i)
  {
    if (tmpSegs[i - 1] != tmpSegs[i])
      retSegs.push_back(tmpSegs[i]);
  }
  return retSegs;
} // MpPolyPts::SegmentsForCleanPolyOffset
//------------------------------------------------------------------------------
/// \brief Intersects the segments
/// \param a_segs Vector of indexes defining segments.
//------------------------------------------------------------------------------
void MePolyPts::IntersectSegs(const std::vector<size_t>& a_segs)
{
  std::vector<GmBstBox3d> envelopes(iCalcSegEnvelopes(a_segs, Pts()));
  for (size_t i = 0; i < a_segs.size(); ++i)
  {
    for (size_t j = i + 1; j < a_segs.size(); ++j)
    {
      // this wouldn't work with the tests with boost 1.55
      // if (bg::overlaps(m_envelopes[i], m_envelopes[j]))
      // if (bg::covered_by(m_envelopes[i], m_envelopes[j]))
      if (iEnvelopesOverlapOrTouch(i, j, envelopes, envelopes))
      { // see if the segments intersect
        CheckIntersectTwoSegs(i, j, a_segs, a_segs);
      }
    }
  }
} // MpPolyPts::IntersectSegs
//------------------------------------------------------------------------------
/// \brief Returns a new sequence that includes intersection indexes
/// \param a_segs List of indexes defining a closed loop polyline.
/// \return a new sequence
//------------------------------------------------------------------------------
std::list<size_t> MePolyPts::SequenceWithIntersects(std::vector<size_t>& a_segs)
{
  std::list<size_t> sequence;

  for (size_t i = 0; i < a_segs.size(); ++i)
  {
    sequence.push_back(a_segs[i]);

    SegMap::iterator it = m_p->m_segCross[&a_segs].find(i);
    if (it == m_p->m_segCross[&a_segs].end())
      continue;

    std::multimap<double, size_t>::iterator it1(it->second.begin());
    for (; it1 != it->second.end(); ++it1)
    {
      size_t& pIdx(it1->second);
      // don't add the same point twice (happens in testCase6c)
      if (sequence.back() != pIdx)
        sequence.push_back(pIdx);
    }
  }
  return sequence;
} // MpPolyPts::SequenceWithIntersects
//------------------------------------------------------------------------------
/// \brief Intersects 2 segments
/// \param a_i Index to a segment.
/// \param a_j Index to a segment.
/// \param a_iSeg Vector of indexes defining a closed loop polyline referred to
/// by a_i
/// \param a_jSeg Vector of indexes defining a closed loop polyline referred to
/// by a_j
//------------------------------------------------------------------------------
void MePolyPts::CheckIntersectTwoSegs(size_t a_i,
                                      size_t a_j,
                                      const std::vector<size_t>& a_iSeg,
                                      const std::vector<size_t>& a_jSeg)
{
  iCheckIntersection(a_i, a_j, a_iSeg, a_jSeg, m_p->m_pts, XyTol(), m_p->m_segCross, m_p->m_ps);
} // MpPolyPts::CheckIntersectTwoSegs
//------------------------------------------------------------------------------
/// \brief Calculates new loops from a sequence that has intersections included
/// \param a_sequence List of indexes defining a closed loop polyline. If the
/// polyline intersects with itself then multiple loops are extracted from
/// the polyline.
/// \param a_loops List of vector of indexes defining a closed loop polyline.
/// These are extracted from a_sequence
//------------------------------------------------------------------------------
void MePolyPts::CalcLoopsForCleanPolyOffset(std::list<size_t>& a_sequence,
                                            std::list<std::vector<size_t>>& a_loops)
{
  while (!a_sequence.empty())
  {
    // keep doing this until we don't find any more
    size_t sSize(0);
    while (sSize != a_sequence.size())
    {
      sSize = a_sequence.size();
      iRemoveRepeatedSegmentFromSequence(a_sequence);
    }

    std::map<size_t, std::list<size_t>::iterator> mapPidxItr;
    std::map<size_t, std::list<size_t>::iterator>::iterator it1;
    std::list<size_t>::iterator it(a_sequence.begin());
    bool loopfound = false;
    for (; !loopfound && it != a_sequence.end(); ++it)
    {
      size_t& ptIdx(*it);
      it1 = mapPidxItr.find(ptIdx);
      if (it1 == mapPidxItr.end())
      {
        mapPidxItr[ptIdx] = it;
      }
      else
      {
        a_loops.push_back(std::vector<size_t>(it1->second, it));
        a_sequence.erase(it1->second, it);
        loopfound = true;
      }
    }
    if (!loopfound)
    {
      a_loops.push_back(std::vector<size_t>(a_sequence.begin(), a_sequence.end()));
      a_sequence.clear();
    }
  }
} // MpPolyPts::CalcLoopsForCleanPolyOffset
//------------------------------------------------------------------------------
/// \brief Removes loops that are "backwards" for CleanPolyOffset. "backwards"
/// depends on the type of input (OUTSIDE_POLY or INSIDE_POLY)
/// \param a_loops List of vector of indexes defining a closed loop polyline.
/// \param a_pType Type of polygon that was offset either OUTSIDE_POLY or
/// INSIDE_POLY.
//------------------------------------------------------------------------------
void MePolyPts::RemoveBackwardLoopsForCleanPolyOffset(std::list<std::vector<size_t>>& a_loops,
                                                      int a_pType)
{
  std::vector<double> areas(iCalcPolyAreas(a_loops, Pts()));
  int cnt(0), out(MePolyOffsetter::OUTSIDE_POLY), in(MePolyOffsetter::INSIDE_POLY);
  std::list<std::vector<size_t>>::iterator it(a_loops.begin()), next;
  while (it != a_loops.end())
  {
    next = it;
    ++next;
    if ((a_pType == out && areas[cnt] > 0) || (a_pType == in && areas[cnt] < 0))
    {
      a_loops.erase(it);
    }
    it = next;
    cnt++;
  }
} // MpPolyPts::RemoveBackwardLoopsForCleanPolyOffset
//------------------------------------------------------------------------------
/// \brief Classifies loops extracted from a pave on an INSIDE_POLY and removes
/// invalid loops. There are more rules when dealing with paving outward from
/// polygons. Sometimes paving outward will create a new "outside" polygon that
/// we then pave inward on the next step in the algorithm. See testCase8. You
/// can also generate "inside" polygons that are inside of other "inside"
/// polygons. These are deleted. Again see testCase8.
/// \param a_loops List of vector of indexes defining a closed loop polyline.
/// \param a_loopType Type of polygon defined in a_loops.
//------------------------------------------------------------------------------
void MePolyPts::ClassifyLoopsFromInPolyAndRemoveInvalid(std::list<std::vector<size_t>>& a_loops,
                                                        std::vector<int>& a_loopType)
{
  std::vector<std::vector<size_t>> polyInsideOfPoly(a_loops.size());
  iFindPolysInsideOfOtherPolys(polyInsideOfPoly, a_loops, *this);

  std::vector<double> areas(iCalcPolyAreas(a_loops, Pts()));
  std::vector<int> validPoly(areas.size(), 1);
  for (size_t i = 0; i < areas.size(); ++i)
  {
    if (areas[i] < 0)
      validPoly[i] = false;
  }
  // if a polygon is not inside of another polygon then just check its area
  // to see if we should remove it
  size_t cnt_it(0);
  std::list<std::vector<size_t>>::iterator it(a_loops.begin()), next;
  for (; it != a_loops.end(); ++it, ++cnt_it)
  {
    if (polyInsideOfPoly[cnt_it].empty()) {}
    else if (!polyInsideOfPoly[cnt_it].empty())
    {
      if (validPoly[polyInsideOfPoly[cnt_it][0]])
      {
        validPoly[cnt_it] = validPoly[cnt_it] ? 0 : 1;
      }
    }
    // else
    //{
    //  XM_ASSERT(false);
    //  XM_LOG(xmlog::debug, "Implement section in "
    //    "PolyCleanerImp::ClassifyLoopsFromInsidePoly");
    //}
  }

  for (it = a_loops.begin(), cnt_it = 0; it != a_loops.end(); it = next, ++cnt_it)
  {
    next = it;
    ++next;
    if (!validPoly[cnt_it])
    {
      a_loops.erase(it);
    }
    else
    {
      int type(MePolyOffsetter::INSIDE_POLY); // INSIDE_POLY
      if (areas[cnt_it] < 0)
        type = MePolyOffsetter::NEWOUT_POLY; // OUTSIDE_POLY
      a_loopType.push_back(type);
    }
  }
} // MpPolyPts::ClassifyLoopsFromInPolyAndRemoveInvalid
//------------------------------------------------------------------------------
/// \brief Returns true if a_polyToTest is inside of a_poly
/// \param a_poly A polygon used to test a_polyToTest
/// \param a_polyToTest Vector of indices defining a polygon. Want to determine
/// if it is inside of a_poly.
/// \return true if a_polyToText is inside of a_poly
//------------------------------------------------------------------------------
bool MePolyPts::PolyInsideOfPoly(const std::vector<size_t>& a_poly,
                                 const std::vector<size_t>& a_polyToTest)
{
  std::vector<Pt3d>& pts(Pts());
  GmBstPoly3d poly;
  for (size_t i = 0; i < a_poly.size(); ++i)
  {
    bg::exterior_ring(poly).push_back(pts[a_poly[i]]);
  }
  bg::exterior_ring(poly).push_back(pts[a_poly[0]]);

  return iPolyInsideOfPoly(poly, a_polyToTest, pts);
} // MePolyPts::PolyInsideOfPoly
//------------------------------------------------------------------------------
/// \brief Returns true if a_polyToTest is inside of a_poly
/// \param a_poly Vector of points defining a polygon. 1st point repeated at
/// the end of the vector.
/// \param a_polyToTest Vector of indices defining a polygon. Want to determine
/// if it is inside of a_poly.
/// \return true if a_polyToText is inside of a_poly
//------------------------------------------------------------------------------
bool MePolyPts::PolyInsideOfPoly(const std::vector<Pt3d>& a_poly,
                                 const std::vector<size_t>& a_polyToTest)
{
  GmBstPoly3d poly;
  for (size_t i = 0; i < a_poly.size(); ++i)
  {
    bg::exterior_ring(poly).push_back(a_poly[i]);
  }
  std::vector<Pt3d>& pts(Pts());
  return iPolyInsideOfPoly(poly, a_polyToTest, pts);
} // MePolyPts::PolyInsideOfPoly
//------------------------------------------------------------------------------
/// \brief Returns the index of a location by hashing that location with the
/// current list of points. If the location does not exist yet then it is
/// added to the list of points.
/// \param a_pt A location.
/// \return index for a_pt.
//------------------------------------------------------------------------------
size_t MePolyPts::IdxFromPt3d(const Pt3d& a_pt)
{
  if (!m_p->m_ps)
  {
    m_p->m_ps = GmPtSearch::New(true);
    m_p->m_ps->VectorThatGrowsToSearch(m_p->m_pts);
  }
  int idx(0);
  m_p->m_ps->AddPtToVectorIfUnique(a_pt, m_p->m_xyTol, idx);
  return idx > -1 ? (size_t)idx : 0;
} // MePolyPts::IdxFromPt3d

} // namespace xms
