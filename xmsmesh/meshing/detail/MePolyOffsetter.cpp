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
#include <xmsmesh/meshing/detail/MePolyOffsetter.h>

// 3. Standard library headers
#include <sstream>

// 4. External library headers
#include <boost/unordered_set.hpp>
#include <xmscore/math/math.h>
#include <xmscore/points/pt.h>

// 5. Shared code headers
#include <xmsinterp/geometry/geoms.h>
#include <xmsmesh/meshing/detail/MePolyCleaner.h>
#include <xmscore/misc/XmError.h>

// 6. Non-shared code headers

//----- Forward declarations ---------------------------------------------------

//----- External globals -------------------------------------------------------

//----- Namespace declaration --------------------------------------------------
namespace xms
{
//----- Constants / Enumerations -----------------------------------------------

//----- Classes / Structs ------------------------------------------------------
class MePolyOffsetterImpl : public MePolyOffsetter
{
public:
  MePolyOffsetterImpl();

  virtual bool Offset(const std::vector<Pt3d>& a_input,
                      MePolyOffsetter::polytype a_pType,
                      MePolyOffsetterOutput& a_output,
                      double a_xyTol) override;

  bool Offset(const std::vector<Pt3d>& a_input,
              MePolyOffsetter::polytype a_pType,
              std::vector<std::vector<Pt3d>>& a_output,
              std::vector<MePolyOffsetter::polytype>& a_outPolyType);
  bool DoOffset(const std::vector<Pt3d>& a_input);
  void CheckToAddPoint(std::vector<Pt3d>& a_result, const Pt3d& a_pt);
  void ProcessAngleSegmentEnd(int npt_end,
                              double ang_end,
                              int in1,
                              int in2,
                              int in3,
                              double dx1,
                              double dy1,
                              const Pt3d* pts,
                              std::vector<Pt3d>& a_result);
  void SpecialRejection(const std::vector<Pt3d>& a_input, std::vector<Pt3d>& a_output);
  void SelfIntersection(std::vector<Pt3d>& a_pLine);
  void FindDuplicatesAndOrderLoops(const std::vector<Pt3d>& a_input);

  double m_xyTol;                    ///< tolerance for geometric comparisons
  bool m_setOffsetToZero;            ///< flag used in testing
  BSHP<MePolyCleaner> m_intersector; ///< class to clean the offset from the polygon
  MePolyOffsetter::polytype m_pType; ///< the type of polygon being offset
  MePolyOffsetterOutput m_output;    ///< the new polygons created by this class
};
//----- Internal functions -----------------------------------------------------

//----- Class / Function definitions -------------------------------------------

namespace
{
#define TWOPI 6.28318530717958647692528 ///< prior calculated constant
#define PIOVER6 0.523598775598299       ///< prior calculated constant
#define SEVENPIOVER3 7.330382858376180  ///< prior calculated constant
#define FOURPIOVER3 4.188790204786390   ///< prior calculated constant
#define SIN60 0.866025403784            ///< prior calculated constant
#define TWOPIOVER3 2.094395102393200    ///< prior calculated constant
#define PIOVER3 1.047197551196600       ///< prior calculated constant

} // unnamed namespace

//------------------------------------------------------------------------------
/// \brief Creates a BufferPoly class
/// \return MePolyOffsetter.
//------------------------------------------------------------------------------
BSHP<MePolyOffsetter> MePolyOffsetter::New()
{
  BSHP<MePolyOffsetter> ptr(new MePolyOffsetterImpl());
  return ptr;
} // PolyOffsetter::New

////////////////////////////////////////////////////////////////////////////////
/// \class MePolyOffsetterImpl
/// \brief Offsets a polyline (in or out). The polyline forms a closed loop.
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief constructor
//------------------------------------------------------------------------------
MePolyOffsetterImpl::MePolyOffsetterImpl()
: m_xyTol(1e-9)
, m_setOffsetToZero(false)
, m_intersector(MePolyCleaner::New())
, m_pType(OUTSIDE_POLY)
{
} // MePolyOffsetterImpl::MePolyOffsetterImpl
//------------------------------------------------------------------------------
/// \brief Takes an input polyline and offsets it based on the distance of the
/// segments that make up the polyline
/// \param a_input The input polyline.
/// \param a_pType The type of polygon: outside polygon or inside polygon
/// \param a_out The offset of the input polyline. This can be multiple
/// polylines
/// \param a_xyTol Tolerance for x,y plane for coincident points.
/// \return Returns false if the input polygon in invalid or the offset fails
//------------------------------------------------------------------------------
bool MePolyOffsetterImpl::Offset(const std::vector<Pt3d>& a_input,
                                 MePolyOffsetter::polytype a_pType,
                                 MePolyOffsetterOutput& a_out,
                                 double a_xyTol)
{
  m_xyTol = a_xyTol;
  m_pType = a_pType;
  std::vector<Pt3d> input(a_input);
  if (!DoOffset(input))
    return false;
  a_out = m_output;
  return true;
} // MePolyOffsetterImpl::Offset
//------------------------------------------------------------------------------
/// \brief Takes an input polyline and offsets it based on the distance of the
/// segments that make up the polyline
/// \param a_input The input polyline.
/// \param a_pType The type of polygon: outside polygon or inside polygon
/// \param a_output The offset of the input polyline. This can be multiple
/// polylines
/// \param a_outPolyType The type of polygon in a_output it is either:
/// OUTSIDE_POLY or INSIDE_POLY
/// \return Returns false if the input polygon in invalid or the offset fails
//------------------------------------------------------------------------------
bool MePolyOffsetterImpl::Offset(const std::vector<Pt3d>& a_input,
                                 MePolyOffsetter::polytype a_pType,
                                 std::vector<std::vector<Pt3d>>& a_output,
                                 std::vector<MePolyOffsetter::polytype>& a_outPolyType)
{
  m_pType = a_pType;
  // check that input polygon is valid
  // if (!InputValid(a_input)) return false;
  std::vector<Pt3d> input(a_input);

  // buffer the polygon
  if (!DoOffset(input))
    return false;

  a_outPolyType.resize(m_output.m_loopTypes.size());
  for (size_t i = 0; i < m_output.m_loopTypes.size(); ++i)
  {
    a_outPolyType[i] = (MePolyOffsetter::polytype)m_output.m_loopTypes[i];
  }
  a_output.assign(m_output.m_loops.size(), std::vector<Pt3d>());
  for (size_t i = 0; i < m_output.m_loops.size(); ++i)
  {
    a_output[i].reserve(m_output.m_loops[i].size());
    for (size_t j = 0; j < m_output.m_loops[i].size(); ++j)
    {
      a_output[i].push_back(m_output.m_pts[m_output.m_loops[i][j]]);
    }
  }

  return true;
} // MePolyOffsetterImpl::Buffer
//------------------------------------------------------------------------------
/// \brief buffers the polygon
/// \param[in] a_input: ???
/// \return true on success.
//------------------------------------------------------------------------------
bool MePolyOffsetterImpl::DoOffset(const std::vector<Pt3d>& a_input)
{
  bool rval(true);
  std::vector<Pt3d> result;

  const Pt3d* pts(&a_input[0]);
  int numpts = (int)a_input.size();
  // don't repeat last point
  if (gmEqualPointsXY(a_input.front(), a_input.back(), m_xyTol))
    --numpts;
  // compute starting angle
  int in1 = numpts - 3;
  int in2 = numpts - 2;
  int in3 = numpts - 1;
  double dx1 = pts[in1].x - pts[in2].x;
  double dy1 = pts[in1].y - pts[in2].y;
  double dx2 = pts[in3].x - pts[in2].x;
  double dy2 = pts[in3].y - pts[in2].y;
  double ang_beg = gmAngleBetween2DVectors(dx1, dy1, dx2, dy2);
  if (EQ_TOL(ang_beg, 0.0, 0.00001))
    ang_beg = TWOPI;
  // compute the number of points to add at the corner based on the angle
  //      0-  <90 degrees - -2 or -3 = block equilateral triangle creation
  //     90- <150 degrees - -1 = usually block equilateral triangle
  //    150- <210 degrees - 0 points around corner
  //    210- <270 degrees - 1 points - bisect corner
  //    270- <330 degrees - 2 points
  //    330-  360 degrees - 3 points
  // computed as the closest integer to (angle-120 degrees)/(60 degrees)
  int npt_beg = int((ang_beg + PIOVER6) / SEVENPIOVER3 * 7.0) - 3;
  in1 = in2;
  in2 = in3;
  // for each segments in the loop
  int num;
  for (in3 = 0, num = 0; in3 < numpts; in3++)
  {
    // find the angle at the end of the segment
    dx1 = pts[in1].x - pts[in2].x;
    dy1 = pts[in1].y - pts[in2].y;
    dx2 = pts[in3].x - pts[in2].x;
    dy2 = pts[in3].y - pts[in2].y;
    double ang_end = gmAngleBetween2DVectors(dx1, dy1, dx2, dy2);
    if (EQ_TOL(ang_end, 0.0, 0.00001))
      ang_end = TWOPI;
    int npt_end = int((ang_end + PIOVER6) / SEVENPIOVER3 * 7.0) - 3;
    // process the segment - the equilateral point
    //    if (npt_beg > -2 && npt_end > -2 && (ang_beg+ang_end > THREEPIOVER2)) {
    if (npt_beg > -2 && npt_end > -2 && (ang_beg + ang_end > FOURPIOVER3))
    {
      Pt3d newpt;
      newpt.x = (pts[in2].x + pts[in1].x) / 2.0 + (pts[in2].y - pts[in1].y) * SIN60;
      newpt.y = (pts[in2].y + pts[in1].y) / 2.0 + (pts[in1].x - pts[in2].x) * SIN60;
      newpt.z = (pts[in2].z + pts[in1].z) / 2.0;
      CheckToAddPoint(result, newpt);
    }
    // process the angle at the end of the segment
    ProcessAngleSegmentEnd(npt_end, ang_end, in1, in2, in3, dx1, dy1, pts, result);

    // set up for next segment
    in1 = in2;
    in2 = in3;
    ang_beg = ang_end;
    npt_beg = npt_end;
  }

  SpecialRejection(a_input, result);
  if (result.size() < 3)
    return false;
  SelfIntersection(result);
  FindDuplicatesAndOrderLoops(result);
  return rval;
} // MePolyOffsetterImpl::DoOffset
//------------------------------------------------------------------------------
/// \brief checks to see if a point can be added to the resulting line
/// \param a_result: ???
/// \param a_pt: The point.
//------------------------------------------------------------------------------
void MePolyOffsetterImpl::CheckToAddPoint(std::vector<Pt3d>& a_result, const Pt3d& a_pt)
{
  a_result.push_back(a_pt);
} // MePolyOffsetterImpl::CheckToAddPoint
//------------------------------------------------------------------------------
/// \brief handles inserting points around the end of a segment.
/// \param[in] npt_end: ???
/// \param[in] ang_end: ???
/// \param[in] in1: ???
/// \param[in] in2: ???
/// \param[in] in3: ???
/// \param[in] dx1: ???
/// \param[in] dy1: ???
/// \param[in] pts: ???
/// \param[in] a_result: ???
//------------------------------------------------------------------------------
void MePolyOffsetterImpl::ProcessAngleSegmentEnd(int npt_end,
                                                 double ang_end,
                                                 int in1,
                                                 int in2,
                                                 int in3,
                                                 double dx1,
                                                 double dy1,
                                                 const Pt3d* pts,
                                                 std::vector<Pt3d>& a_result)
{
  Pt3d newpt;
  double l1, l2, dl, offset, ptx, pty, ctheta, stheta, alpha;

  switch (npt_end)
  {
  case 1:
    l1 = Mdist(pts[in1].x, pts[in1].y, pts[in2].x, pts[in2].y);
    l2 = Mdist(pts[in3].x, pts[in3].y, pts[in2].x, pts[in2].y);
    offset = (l1 + l2) * 0.4330127019 / l1;
    ptx = pts[in2].x + offset * dx1;
    pty = pts[in2].y + offset * dy1;
    ctheta = cos(ang_end / 2);
    stheta = sin(ang_end / 2);
    newpt.x = ptx * ctheta - pty * stheta + (1.0 - ctheta) * pts[in2].x + stheta * pts[in2].y;
    newpt.y = ptx * stheta + pty * ctheta + (1.0 - ctheta) * pts[in2].y - stheta * pts[in2].x;
    newpt.z = (pts[in2].z + pts[in1].z) / 2.0;
    CheckToAddPoint(a_result, newpt);
    break;
  case 2:
    l1 = Mdist(pts[in1].x, pts[in1].y, pts[in2].x, pts[in2].y);
    l2 = Mdist(pts[in3].x, pts[in3].y, pts[in2].x, pts[in2].y);
    dl = l2 - l1;
    // compute angle between points
    alpha = (ang_end - TWOPIOVER3) / 3.0;
    // place first point
    offset = SIN60 * (l1 + dl * ((alpha - PIOVER6) / (3 * alpha - PIOVER3))) / l1;
    ptx = pts[in2].x + offset * dx1;
    pty = pts[in2].y + offset * dy1;
    ctheta = cos(alpha + PIOVER3);
    stheta = sin(alpha + PIOVER3);
    newpt.x = ptx * ctheta - pty * stheta + (1.0 - ctheta) * pts[in2].x + stheta * pts[in2].y;
    newpt.y = ptx * stheta + pty * ctheta + (1.0 - ctheta) * pts[in2].y - stheta * pts[in2].x;
    newpt.z = (pts[in2].z + pts[in1].z) / 2.0;
    CheckToAddPoint(a_result, newpt);
    // place second point
    offset = SIN60 * (l1 + dl * ((2 * alpha - PIOVER6) / (3 * alpha - PIOVER3))) / l1;
    ptx = pts[in2].x + offset * dx1;
    pty = pts[in2].y + offset * dy1;
    ctheta = cos(2 * alpha + PIOVER3);
    stheta = sin(2 * alpha + PIOVER3);
    newpt.x = ptx * ctheta - pty * stheta + (1.0 - ctheta) * pts[in2].x + stheta * pts[in2].y;
    newpt.y = ptx * stheta + pty * ctheta + (1.0 - ctheta) * pts[in2].y - stheta * pts[in2].x;
    newpt.z = (pts[in2].z + pts[in1].z) / 2.0;
    CheckToAddPoint(a_result, newpt);
    break;
  case 3:
    l1 = Mdist(pts[in1].x, pts[in1].y, pts[in2].x, pts[in2].y);
    l2 = Mdist(pts[in3].x, pts[in3].y, pts[in2].x, pts[in2].y);
    dl = l2 - l1;
    // compute angle between points
    alpha = (ang_end - TWOPIOVER3) / 4.0;
    // place first point
    offset = SIN60 * (l1 + dl * ((alpha - PIOVER6) / (4 * alpha - PIOVER3))) / l1;
    ptx = pts[in2].x + offset * dx1;
    pty = pts[in2].y + offset * dy1;
    ctheta = cos(alpha + PIOVER3);
    stheta = sin(alpha + PIOVER3);
    newpt.x = ptx * ctheta - pty * stheta + (1.0 - ctheta) * pts[in2].x + stheta * pts[in2].y;
    newpt.y = ptx * stheta + pty * ctheta + (1.0 - ctheta) * pts[in2].y - stheta * pts[in2].x;
    newpt.z = (pts[in2].z + pts[in1].z) / 2.0;
    CheckToAddPoint(a_result, newpt);
    // place second point
    offset = SIN60 * (l1 + dl * ((2 * alpha - PIOVER6) / (4 * alpha - PIOVER3))) / l1;
    ptx = pts[in2].x + offset * dx1;
    pty = pts[in2].y + offset * dy1;
    ctheta = cos(2 * alpha + PIOVER3);
    stheta = sin(2 * alpha + PIOVER3);
    newpt.x = ptx * ctheta - pty * stheta + (1.0 - ctheta) * pts[in2].x + stheta * pts[in2].y;
    newpt.y = ptx * stheta + pty * ctheta + (1.0 - ctheta) * pts[in2].y - stheta * pts[in2].x;
    newpt.z = (pts[in2].z + pts[in1].z) / 2.0;
    CheckToAddPoint(a_result, newpt);
    // place third point
    offset = SIN60 * (l1 + dl * ((3 * alpha - PIOVER6) / (4 * alpha - PIOVER3))) / l1;
    ptx = pts[in2].x + offset * dx1;
    pty = pts[in2].y + offset * dy1;
    ctheta = cos(3 * alpha + PIOVER3);
    stheta = sin(3 * alpha + PIOVER3);
    newpt.x = ptx * ctheta - pty * stheta + (1.0 - ctheta) * pts[in2].x + stheta * pts[in2].y;
    newpt.y = ptx * stheta + pty * ctheta + (1.0 - ctheta) * pts[in2].y - stheta * pts[in2].x;
    newpt.z = (pts[in2].z + pts[in1].z) / 2.0;
    CheckToAddPoint(a_result, newpt);
    break;
  }
} // MePolyOffsetterImpl::ProcessAngleSegmentEnd
//------------------------------------------------------------------------------
/// \brief a special case for removing points.
/// \param[in] a_input: ???
/// \param a_out: ???
//------------------------------------------------------------------------------
void MePolyOffsetterImpl::SpecialRejection(const std::vector<Pt3d>& a_input,
                                           std::vector<Pt3d>& a_out)
{
  // See if the beginning node needs to be deleted due to being a sharp angle
  /****************************************************
   *  AKZ-FIX - look for special rejection case
   *  - a thin triangle can offset to a triangle almost
   *    as big as the original with two or three edges
   *    intersecting of the original.
   *  We may want to handle all cases of the offset
   *    intersecting the original, but for now just
   *    delete this special case
   ****************************************************/
  const Pt3d* pts(&a_input[0]);
  int in(0);
  if (a_input.size() == 3 && a_out.size() == 3)
  {
    in = (int)gmPointInTriangleWithTol(&pts[0], &pts[1], &pts[2], a_out[0].x, a_out[0].y, m_xyTol);
    in += (int)gmPointInTriangleWithTol(&pts[0], &pts[1], &pts[2], a_out[1].x, a_out[1].y, m_xyTol);
    in += (int)gmPointInTriangleWithTol(&pts[0], &pts[1], &pts[2], a_out[2].x, a_out[2].y, m_xyTol);
    if (in < 2)
      a_out.resize(0);
  }
} // MePolyOffsetterImpl::SpecialRejection
//------------------------------------------------------------------------------
/// \brief Self intersection of a polyline.
/// \param a_pLine: The polyline.
//------------------------------------------------------------------------------
void MePolyOffsetterImpl::SelfIntersection(std::vector<Pt3d>& a_pLine)
{
  m_intersector->CleanPolyOffset(a_pLine, m_pType, m_xyTol, m_output);
} // MePolyOffsetterImpl::SelfIntersection
//------------------------------------------------------------------------------
/// \brief Finds duplicate points and orders the loops
/// \param[in] a_input Vector of point locations
//------------------------------------------------------------------------------
void MePolyOffsetterImpl::FindDuplicatesAndOrderLoops(const std::vector<Pt3d>& a_input)
{
  // figure out which points are repeated in the loops (these will be nodes
  // in a coverage)
  std::vector<int> vCnt(m_output.m_pts.size(), 0);
  for (size_t i = 0; i < m_output.m_loops.size(); ++i)
  {
    for (size_t j = 0; j < m_output.m_loops[i].size(); ++j)
    {
      vCnt[m_output.m_loops[i][j]]++;
    }
  }
  boost::unordered_set<size_t> setIdx;
  for (size_t i = 0; i < vCnt.size(); ++i)
  {
    if (vCnt[i] > 1)
    {
      setIdx.insert(i);
    }
  }
  // add all of the intersection points
  for (size_t i = a_input.size(); i < m_output.m_pts.size(); ++i)
  {
    if (setIdx.find(i) == setIdx.end() && vCnt[i] > 0)
    {
      setIdx.insert(i);
    }
  }

  for (size_t i = 0; i < m_output.m_loops.size(); ++i)
  {
    bool found = false;
    for (size_t j = 0; !found && j < m_output.m_loops[i].size(); ++j)
    {
      if (setIdx.find(m_output.m_loops[i][j]) != setIdx.end())
        found = true;
    }
    if (!found)
    {
      setIdx.insert(m_output.m_loops[i][0]);
    }
  }
  // order the points in the loops so that they start with one of the points
  // in setIdx
  for (size_t i = 0; i < m_output.m_loops.size(); ++i)
  {
    if (setIdx.find(m_output.m_loops[i][0]) == setIdx.end())
    {
      size_t start = 1;
      bool found = false;
      for (size_t j = 1; !found && j < m_output.m_loops[i].size(); ++j)
      {
        if (setIdx.find(m_output.m_loops[i][j]) != setIdx.end())
        {
          found = true;
          start = j;
        }
      }
      std::vector<size_t> l = m_output.m_loops[i];
      size_t cnt(0);
      while (cnt < l.size())
      {
        l[cnt] = m_output.m_loops[i][start];
        cnt++;
        start++;
        if (start >= l.size())
          start = 0;
      }
      m_output.m_loops[i] = l;
    }
  }
} // MePolyOffsetterImpl::FindDuplicatesAndOrderLoops

} // namespace xms

#ifdef CXX_TEST
////////////////////////////////////////////////////////////////////////////////

#include <xmsmesh/meshing/detail/MePolyOffsetter.t.h>

//#include <boost/assign.hpp>

#include <xmscore/testing/TestTools.h>

// namespace xms {
using namespace xms;

////////////////////////////////////////////////////////////////////////////////
/// \class MePolyOffsetterUnitTests
/// \brief tester for the MePolyOffsetter class
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief tests creating the class
//------------------------------------------------------------------------------
void MePolyOffsetterUnitTests::testCreateClass()
{
  BSHP<MePolyOffsetter> b = MePolyOffsetter::New();
  TS_ASSERT(b);
} // PolyOffsetterTest::testCreateClass
//------------------------------------------------------------------------------
/// \brief test offsetting from a box polygon
//------------------------------------------------------------------------------
void MePolyOffsetterUnitTests::testBox0()
{
  // clang-format off
  //   This is the input
  //
  // x=        0    4    8    12
  //
  // y=12     2--- 3--- 4--- 5
  //          |              |
  //          |              |
  // y-8      1              6
  //          |              |
  //          |              |
  // y=4      0              7
  //          |              |
  //          |              |
  // y=0     11---10--- 9--- 8
  //
  // clang-format on
  MePolyOffsetter::polytype ptype = MePolyOffsetter::OUTSIDE_POLY;
  MePolyOffsetterImpl pl;
  pl.m_setOffsetToZero = true;
  std::vector<Pt3d> input = {{0, 4, 0},  {0, 8, 0},   {0, 12, 0}, {4, 12, 0},
                             {8, 12, 0}, {12, 12, 0}, {12, 8, 0}, {12, 4, 0},
                             {12, 0, 0}, {8, 0, 0},   {4, 0, 0},  {0, 0, 0}};
  std::vector<std::vector<Pt3d>> output;
  std::vector<MePolyOffsetter::polytype> outPolyTypes;
  pl.Offset(input, ptype, output, outPolyTypes);
  std::vector<Pt3d> base = {{3.4641, 3.4641, 0.0}, {3.4641, 6.0, 0.0},    {3.4641, 8.5359, 0.0},
                            {6.0, 8.5359, 0.0},    {8.5359, 8.5359, 0.0}, {8.5359, 6.0, 0.0},
                            {8.5359, 3.4641, 0.0}, {6.0, 3.4641, 0.0}};
  TS_ASSERT_DELTA_VECPT3D(base, output[0], 1e-4);
} // PolyOffsetterTest::testBox0
//------------------------------------------------------------------------------
/// \brief tests offseting from the shape shown below
//------------------------------------------------------------------------------
void MePolyOffsetterUnitTests::testBox1()
{
  // clang-format off
  //   This is the input
  //
  // x=        0    4    8    12
  //
  // y=12     2--- 3--- 4--- 5
  //          |              |
  //          |              |
  // y-8      1              6
  //          |           /
  //          |         7
  //          |           \
  //y=4       0              8
  //          |              |
  //          |              |
  // y=0     12---11---10--- 9
  //
  // clang-format on
  MePolyOffsetter::polytype ptype = MePolyOffsetter::OUTSIDE_POLY;
  MePolyOffsetterImpl pl;
  pl.m_setOffsetToZero = true;
  std::vector<Pt3d> input = {{0, 4, 0},   {0, 8, 0},  {0, 12, 0}, {4, 12, 0}, {8, 12, 0},
                             {12, 12, 0}, {12, 8, 0}, {8, 6, 0},  {12, 4, 0}, {12, 0, 0},
                             {8, 0, 0},   {4, 0, 0},  {0, 0, 0}};
  std::vector<std::vector<Pt3d>> output;
  std::vector<MePolyOffsetter::polytype> outPolyTypes;
  pl.Offset(input, ptype, output, outPolyTypes);
  std::vector<Pt3d> base = {{3.4641, 3.4641, 0.0}, {3.4641, 6.0, 0.0},    {3.4641, 8.5359, 0.0},
                            {5.4609, 8.5359, 0.0}, {4.6853, 8.0031, 0.0}, {4.6853, 3.9969, 0.0},
                            {5.4609, 3.4641, 0.0}};
  TS_ASSERT_DELTA_VECPT3D(base, output[0], 1e-4);
} // PolyOffsetterTest::testBox1
//------------------------------------------------------------------------------
/// \brief tests offseting from the shape shown below
//------------------------------------------------------------------------------
void MePolyOffsetterUnitTests::testBox1a()
{
  // clang-format off
  // x=        0    4    8    12
  //
  // y=12    11---10--- 9--- 8
  //          |              |
  //          |              |
  // y-8      12              7
  //          |           /
  //          |         6
  //          |           \
  //y=4       0              5
  //          |              |
  //          |              |
  // y=0      1--- 2--- 3--- 4
  //
  // clang-format on
  MePolyOffsetter::polytype ptype = MePolyOffsetter::INSIDE_POLY;
  MePolyOffsetterImpl pl;
  pl.m_setOffsetToZero = true;
  std::vector<Pt3d> input = {{0, 4, 0},  {0, 0, 0},  {4, 0, 0},  {8, 0, 0},   {12, 0, 0},
                             {12, 4, 0}, {8, 6, 0},  {12, 8, 0}, {12, 12, 0}, {8, 12, 0},
                             {4, 12, 0}, {0, 12, 0}, {0, 8, 0}};
  std::vector<std::vector<Pt3d>> output;
  std::vector<MePolyOffsetter::polytype> outPolyTypes;
  pl.Offset(input, ptype, output, outPolyTypes);
  std::vector<Pt3d> base = {{-3.4641, 10, 0},      {-3.4641, 6, 0},        {-3.4641, 2, 0},
                            {-3.2551, -1.1847, 0}, {-1.1847, -3.2551, 0},  {2, -3.4641, 0},
                            {6, -3.4641, 0},       {10, -3.4641, 0},       {13.1847, -3.2551, 0},
                            {15.2551, -1.1847, 0}, {15.4641, 2, 0},        {15.1206, 5.9286, 0},
                            {15.1206, 6.0713, 0},  {15.4641, 10, 0},       {15.2551, 13.1847, 0},
                            {13.1847, 15.2551, 0}, {10, 15.4641, 0},       {6, 15.4641, 0},
                            {2, 15.4641, 0},       {-1.18479, 15.2551, 0}, {-3.2551, 13.1847, 0}};
  TS_ASSERT_DELTA_VECPT3D(base, output[0], 1e-4);
} // PolyOffsetterTest::testBox1
//------------------------------------------------------------------------------
/// \brief tests offseting from the shape shown below
//------------------------------------------------------------------------------
void MePolyOffsetterUnitTests::testCase1()
{
  // clang-format off
  // x=        0      4      8      12
  //
  // y=12    18-----17-----16----- 15
  //          |                    |
  //          |                    |
  //          |                    |
  //          |                    |
  // y-8     19         11--12--13-14
  //          |          |
  //          |         10
  //          |          |
  //          |          9
  //          |          |
  // y=4      0          8---7--6--5
  //          |                    |
  //          |                    |
  //          |                    |
  //          |                    |
  // y=0      1----- 2----- 3----- 4
  //
  // clang-format on
  MePolyOffsetter::polytype ptype = MePolyOffsetter::INSIDE_POLY;
  MePolyOffsetterImpl pl;
  pl.m_setOffsetToZero = true;
  std::vector<Pt3d> input{{0, 4, 0},    {0, 0, 0},    {4, 0, 0},  // 2
                          {8, 0, 0},    {12, 0, 0},   {12, 4, 0}, // 5
                          {10, 4, 0},   {8, 4, 0},    {6, 4, 0},  // 8
                          {6, 5.33, 0}, {6, 6.66, 0}, {6, 8, 0},  // 11
                          {8, 8, 0},    {10, 8, 0},   {12, 8, 0}, // 14
                          {12, 12, 0},  {8, 12, 0},   {4, 12, 0}, // 17
                          {0, 12, 0},   {0, 8, 0}};
  std::vector<std::vector<Pt3d>> output;
  std::vector<MePolyOffsetter::polytype> outPolyTypes;
  pl.Offset(input, ptype, output, outPolyTypes);
  std::vector<Pt3d> base{{7.1518, 5.7321, 0.0},   {7.1518, 5.9950, 0.0},   {7.1536, 6.2679, 0.0},
                         {9.0, 6.2679, 0.0},      {11.0, 6.2679, 0.0},     {12.7240, 6.0107, 0.0},
                         {14.8935, 6.9469, 0.0},  {15.4641, 10.0, 0.0},    {15.2552, 13.1848, 0.0},
                         {13.1848, 15.2552, 0.0}, {10.0, 15.4641, 0.0},    {6.0, 15.4641, 0.0},
                         {2.0, 15.4641, 0.0},     {-1.1848, 15.2552, 0.0}, {-3.2552, 13.1848, 0.0},
                         {-3.4641, 10.0, 0.0},    {-3.4641, 6.0, 0.0},     {-3.4641, 2.0, 0.0},
                         {-3.2552, -1.1848, 0.0}, {-1.1848, -3.2552, 0.0}, {2.0, -3.4641, 0.0},
                         {6.0, -3.4641, 0.0},     {10.0, -3.4641, 0.0},    {13.1848, -3.2552, 0.0},
                         {15.2552, -1.1848, 0.0}, {15.4641, 2.0, 0.0},     {14.8935, 5.0531, 0.0},
                         {12.7240, 5.9893, 0.0},  {11.0000, 5.7321, 0.0},  {9.0000, 5.7321, 0.0}};
  TS_ASSERT_DELTA_VECPT3D(base, output[0], 1e-4);
} // PolyOffsetterTest::testCase1
//------------------------------------------------------------------------------
/// \brief tests offseting from the shape shown below
//------------------------------------------------------------------------------
void MePolyOffsetterUnitTests::testCase1a()
{
  // clang-format off
  // x=        0      4      8      12
  //
  // y=12    18-----17-----16----- 15
  //          |                    |
  //          |                    |
  //          |                    |
  //          |                    |
  // y-8     19         11--12--13 |
  //          |          |        \14
  //          |         10
  //          |          |
  //          |          9
  //          |          |         5
  // y=4      0          8---7--6/ |
  //          |                    |
  //          |                    |
  //          |                    |
  //          |                    |
  // y=0      1----- 2----- 3----- 4
  //
  // clang-format on
  MePolyOffsetter::polytype ptype = MePolyOffsetter::INSIDE_POLY;
  MePolyOffsetterImpl pl;
  pl.m_setOffsetToZero = true;
  std::vector<Pt3d> input{{0, 4, 0},    {0, 0, 0},    {4, 0, 0},  // 2
                          {8, 0, 0},    {12, 0, 0},   {12, 5, 0}, // 5
                          {10, 4, 0},   {8, 4, 0},    {6, 4, 0},  // 8
                          {6, 5.33, 0}, {6, 6.66, 0}, {6, 8, 0},  // 11
                          {8, 8, 0},    {10, 8, 0},   {12, 7, 0}, // 14
                          {12, 12, 0},  {8, 12, 0},   {4, 12, 0}, // 17
                          {0, 12, 0},   {0, 8, 0}};
  std::vector<std::vector<Pt3d>> output;
  std::vector<MePolyOffsetter::polytype> outPolyTypes;
  pl.Offset(input, ptype, output, outPolyTypes);
  TS_ASSERT_EQUALS(2, output.size());
  if (2 != output.size())
    return;
  std::vector<Pt3d> base{{9.6076, 6.0, 0.0},   {9, 5.7320, 0.0},      {7.1518, 5.7320, 0.0},
                         {7.1518, 5.995, 0.0}, {7.1535, 6.2679, 0.0}, {9, 6.2679, 0.0}};
  TS_ASSERT_DELTA_VECPT3D(base, output[0], 1e-4);
  base = {{15.4709, 6.0, 0.0},     {16.3301, 9.5, 0.0},     {15.8881, 13.4151, 0.0},
          {13.2506, 15.4360, 0.0}, {10.0, 15.4641, 0.0},    {6.0, 15.4641, 0.0},
          {2.0, 15.4641, 0.0},     {-1.1848, 15.2551, 0.0}, {-3.2552, 13.1847, 0.0},
          {-3.4641, 10.0, 0.0},    {-3.4641, 6.0, 0.0},     {-3.4641, 2.0, 0.0},
          {-3.2552, -1.1847, 0.0}, {-1.1848, -3.2551, 0.0}, {2.0, -3.4641, 0.0},
          {6.0, -3.4641, 0.0},     {10.0, -3.4641, 0.0},    {13.2506, -3.4360, 0.0},
          {15.8881, -1.4151, 0.0}, {16.3301, 2.5, 0.0}};
  TS_ASSERT_DELTA_VECPT3D(base, output[1], 1e-4);
  std::vector<MePolyOffsetter::polytype> basePolyTypes = {MePolyOffsetter::NEWOUT_POLY,
                                                          MePolyOffsetter::INSIDE_POLY};
  TS_ASSERT_EQUALS_VEC(basePolyTypes, outPolyTypes);
} // PolyOffsetterTest::testCase1a
//------------------------------------------------------------------------------
/// \brief tests offseting from the shape shown below
//------------------------------------------------------------------------------
void MePolyOffsetterUnitTests::testCase1b()
{
  // clang-format off
  // x=        0      4      8      12
  //
  // y=12    18-----17-----16----- 15
  //          |                    |
  //          |                    |
  //          |                    |
  //          |                    |
  // y-8     19         11--12--13 |
  //          |          |        \14
  //          |         10
  //          |          |
  //          |          9
  //          |          |         5
  // y=4      0          8---7--6/ |
  //          |                    |
  //          |                    |
  //          |                    |
  //          |                    |
  // y=0      1----- 2----- 3----- 4
  //
  // clang-format on
  MePolyOffsetter::polytype ptype = MePolyOffsetter::INSIDE_POLY;
  MePolyOffsetterImpl pl;
  pl.m_setOffsetToZero = true;
  std::vector<Pt3d> input{{0, 4, 0},    {0, 0, 0},    {4, 0, 0},  // 2
                          {8, 0, 0},    {12, 0, 0},   {12, 5, 0}, // 5
                          {10, 4, 0},   {8, 4, 0},    {6, 4, 0},  // 8
                          {6, 5.33, 0}, {6, 6.66, 0}, {6, 8, 0},  // 11
                          {8, 8, 0},    {10, 8, 0},   {12, 7, 0}, // 14
                          {12, 12, 0},  {8, 12, 0},   {4, 12, 0}, // 17
                          {0, 12, 0},   {0, 8, 0}};
  MePolyOffsetterOutput out;
  pl.Offset(input, ptype, out, 1e-9);
  TS_ASSERT_EQUALS(2, out.m_loops.size());
  if (2 != out.m_loops.size())
    return;
  std::vector<Pt3d> base{{-3.4641, 10, 0.0},      {-3.4641, 6, 0.0},       {-3.4641, 2, 0.0},
                         {-3.2552, -1.1847, 0.0}, {-1.1848, -3.2551, 0.0}, {2.0, -3.4641, 0.0},
                         {6.0, -3.4641, 0.0},     {10.0, -3.4641, 0.0},    {13.2506, -3.4360, 0.0},
                         {15.8881, -1.4151, 0.0}, {16.3301, 2.5, 0.0},     {15.2735, 6.8037, 0.0},
                         {12.1011, 7.5270, 0.0},  {10.1340, 6.2320, 0.0},  {9.0, 5.7320, 0.0},
                         {7.0, 5.7320, 0.0},      {7.1518, 4.665, 0.0},    {7.1518, 5.995, 0.0},
                         {7.1605, 7.33, 0.0},     {7.0, 6.2679, 0.0},      {9.0, 6.2679, 0.0},
                         {10.1340, 5.7679, 0.0},  {12.1011, 4.4730, 0.0},  {15.2735, 5.1963, 0.0},
                         {16.3301, 9.5, 0.0},     {15.8881, 13.4152, 0.0}, {13.2506, 15.4360, 0.0},
                         {10.0, 15.4641, 0.0},    {6.0, 15.4641, 0.0},     {2.0, 15.4641, 0.0},
                         {-1.1847, 15.2552, 0.0}, {-3.2551, 13.1848, 0.0}, {15.4709, 6.0, 0.0},
                         {9.6077, 6.0, 0.0},      {7.1518, 5.7321, 0.0},   {7.1536, 6.2679, 0.0}};
  TS_ASSERT_DELTA_VECPT3D(base, out.m_pts, 1e-4);

  std::vector<MePolyOffsetter::polytype> basePolyTypes = {MePolyOffsetter::NEWOUT_POLY,
                                                          MePolyOffsetter::INSIDE_POLY};
  TS_ASSERT_EQUALS_VEC(basePolyTypes, out.m_loopTypes);

  std::vector<size_t> baseIdx = {33, 14, 34, 17, 35, 20};
  TS_ASSERT_EQUALS_VEC(baseIdx, out.m_loops[0]);
  baseIdx = {32, 24, 25, 26, 27, 28, 29, 30, 31, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  TS_ASSERT_EQUALS_VEC(baseIdx, out.m_loops[1]);
} // PolyOffsetterTest::testCase1b

//} // namespace xms
#endif
