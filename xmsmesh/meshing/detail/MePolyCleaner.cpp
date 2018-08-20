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
#include <xmsmesh/meshing/detail/MePolyCleaner.h>

// 3. Standard library headers
#include <iostream>

// 4. External library headers

// 5. Shared code headers
#include <xmsmesh/meshing/detail/MePolyOffsetter.h>
#include <xmsmesh/meshing/detail/MeIntersectPolys.h>
#include <xmsmesh/meshing/detail/MePolyPts.h>
#include <xmscore/misc/XmError.h>

// 6. Non-shared code headers

//----- Forward declarations ---------------------------------------------------
using namespace xms;

//----- External globals -------------------------------------------------------

//----- Namespace declaration --------------------------------------------------
namespace xms
{
//----- Constants / Enumerations -----------------------------------------------

//----- Classes / Structs ------------------------------------------------------
class MePolyCleanerImpl : public MePolyCleaner
{
public:
  MePolyCleanerImpl() {}

  virtual void CleanPolyOffset(const std::vector<Pt3d>& a_input,
                               int a_pType,
                               double a_tol,
                               MePolyOffsetterOutput& a_out) override;
  virtual void IntersectCleanInPolys(const std::vector<MePolyOffsetterOutput>& a_offsets,
                                     MePolyOffsetterOutput& a_out,
                                     double a_xyTol) override;
  virtual void IntersectCleanInOutPolys(const MePolyOffsetterOutput& a_offsets,
                                        MePolyOffsetterOutput& a_out,
                                        double a_xyTol) override;

  void FillOutputForCleanPolyOffset(MePolyOffsetterOutput& a_out,
                                    int a_pType,
                                    std::list<std::vector<size_t>>& a_loops,
                                    std::vector<int>& a_loopType,
                                    const std::vector<Pt3d>& a_pts);
};
//----- Internal functions -----------------------------------------------------
#if 0
//------------------------------------------------------------------------------
/// \brief Dumps the output loops and points to std::out for ease in debugging failed tests.
//------------------------------------------------------------------------------
static void iDumpOutput(const MePolyOffsetterOutput& a_output)
{
    std::cout << "o1.m_loops.size(): " << a_output.m_loops.size() << std::endl;
    int i  = -1;
    for (auto loop : a_output.m_loops)
    {
      ++i;
      std::cout << "loop type: " << a_output.m_loopTypes[i] << " loop number: " << i << std::endl;
      for (auto idx : loop)
      {
        std::cout << idx << ", ";
      }
      std::cout << std::endl << std::endl;

      for (auto idx : loop)
      {
        std::cout << "m_pts[" << idx << "]: " << "{" << a_output.m_pts[idx].x << ", "
          << a_output.m_pts[idx].y << "}" << std::endl;
      }
    }
}
#endif

//----- Class / Function definitions -------------------------------------------

//------------------------------------------------------------------------------
/// \brief Creates a new instance of this class
/// \return MePolyCleaner.
//------------------------------------------------------------------------------
BSHP<MePolyCleaner> MePolyCleaner::New()
{
  BSHP<MePolyCleaner> ptr(new MePolyCleanerImpl());
  return ptr;
} // PolyOffsetter::New

////////////////////////////////////////////////////////////////////////////////
/// \class MePolyCleanerImpl
/// \brief Intersects and cleans polygons generated from MePolyOffsetter
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief Takes an input polyline (forming a closed loop) and intersects it
/// with itself and fills the output with multiple closed loop polylines
/// \param a_input The input polyline forming a closed loop.
/// \param a_pType The type of polygon: OUTSIDE_POLY or INSIDE_POLY
/// \param a_tol Tolerance used for floating point geometry comparisons.
/// \param a_out Output class to hold new polygons
//------------------------------------------------------------------------------
void MePolyCleanerImpl::CleanPolyOffset(const std::vector<Pt3d>& a_input,
                                        int a_pType,
                                        double a_tol,
                                        MePolyOffsetterOutput& a_out)
{
  MePolyPts polyPts;
  polyPts.XyTol() = a_tol;
  polyPts.Pts() = a_input;
  // create a vector of segments - vector with the indexes to the points in
  // m_polyPts.m_pts
  std::vector<size_t> segs(polyPts.SegmentsForCleanPolyOffset());

  // Intersect the segments
  polyPts.IntersectSegs(segs);

  // Generate a list of segments that includes the intersections just
  // calculated. We use a list because we are going to iteratively remove
  // elements from the list and we didn't want to force a vector to keep
  // copying stuff around when we delete items in the list.
  std::list<size_t> sequence(polyPts.SequenceWithIntersects(segs));

  // create closed loops (polygons) using the list of segments. One polyline
  // came into this method "a_input". That line may have intersected itself
  // and is really multiple polygons. This next method extracts those multiple
  // polygons
  std::list<std::vector<size_t>> loops;
  polyPts.CalcLoopsForCleanPolyOffset(sequence, loops);

  // clean up the polygons that were extracted
  std::vector<int> loopType;
  if (MePolyOffsetter::OUTSIDE_POLY == a_pType || loops.size() < 2)
  { // when we pave inward or if we only extracted 1 polygon then
    // we delete any polygons with an area that has
    // the opposite sign from what we are expecting
    polyPts.RemoveBackwardLoopsForCleanPolyOffset(loops, a_pType);
  }
  else
  { // There are more rules when dealing with paving outward from polygons.
    // Sometimes paving outward will create a new "outside" polygon that we
    // then pave inward on the next step in the algorithm. See testCase8.
    // You can also generate "inside" polygons that are inside of other
    // "inside" polygons. These are deleted. Again see testCase8.
    polyPts.ClassifyLoopsFromInPolyAndRemoveInvalid(loops, loopType);
  }
  FillOutputForCleanPolyOffset(a_out, a_pType, loops, loopType, polyPts.Pts());
} // MePolyCleanerImpl::CleanPolyOffset
//------------------------------------------------------------------------------
/// \brief Takes a vector of PolyOffsetOutput and intersects and cleans all
/// INSIDE_POLY polygons .
/// \param a_offsets Results from the PolyOffsetter class
/// \param a_out Results from this operation. These will be the INSIDE_POLY
/// polygons considered in future operations done by the meshing class. There
/// can also be new OUTSIDE_POLY polygons created by this method. Any old
/// OUTSIDE_POLY polygons that are passed in will be passed back out.
/// \param a_xyTol Tolerance used for floating point geometry comparisons.
//------------------------------------------------------------------------------
void MePolyCleanerImpl::IntersectCleanInPolys(const std::vector<MePolyOffsetterOutput>& a_offsets,
                                              MePolyOffsetterOutput& a_out,
                                              double a_xyTol)
{
  MeIntersectPolys inPolys;
  inPolys.SetupInIn(a_offsets, a_xyTol);
  // compute envelopes for all INSIDE_POLY polygons
  inPolys.CalcEnvelopes();
  // do trivial cases of completely disjoint or completely inside
  inPolys.InInTrivialPolyCases();
  // intersect segments of potentially overlapping polys
  inPolys.InInDoIntersection();
  inPolys.FillOutput(a_out);
} // MePolyCleanerImpl::IntersectCleanInPolys
//------------------------------------------------------------------------------
/// \brief Takes a PolyOffsetOutput and intersects and cleans all
/// INSIDE_POLY polygons with OUTSIDE_POLY polygons. NEWOUT_POLY polygons
/// are not changed.
/// \param a_offsets Results from the IntersectCleanInPolys method
/// \param a_out Results from this operation.
/// \param a_xyTol Tolerance used for floating point geometry comparisons.
//------------------------------------------------------------------------------
void MePolyCleanerImpl::IntersectCleanInOutPolys(const MePolyOffsetterOutput& a_offsets,
                                                 MePolyOffsetterOutput& a_out,
                                                 double a_xyTol)
{
  MeIntersectPolys inPolys;
  inPolys.SetupInOut(a_offsets, a_xyTol);
  inPolys.CalcEnvelopes();
  inPolys.InOutDoIntersection();
  inPolys.FillOutput(a_out);
} // MePolyCleanerImpl::IntersectCleanInPolys
//------------------------------------------------------------------------------
/// \brief Fills the output variable for CleanPolyOffset method.
/// \param a_out:
/// \param a_pType:
/// \param a_loops:
/// \param a_loopType:
/// \param a_pts:
//------------------------------------------------------------------------------
void MePolyCleanerImpl::FillOutputForCleanPolyOffset(MePolyOffsetterOutput& a_out,
                                                     int a_pType,
                                                     std::list<std::vector<size_t>>& a_loops,
                                                     std::vector<int>& a_loopType,
                                                     const std::vector<Pt3d>& a_pts)
{
  a_out.m_pts = a_pts;
  a_out.m_loops.resize(0);
  a_out.m_loops.reserve(a_loops.size());
  std::list<std::vector<size_t>>::iterator it(a_loops.begin());
  for (; it != a_loops.end(); ++it)
  {
    a_out.m_loops.push_back(std::vector<size_t>(it->begin(), it->end()));
  }
  std::vector<int>& lt(a_out.m_loopTypes);
  if (!a_loopType.empty())
    lt.swap(a_loopType);
  else
    lt.assign(a_loops.size(), a_pType);
} // MePolyCleanerImpl::FillOutputForCleanPolyOffset

} // namespace xms

#ifdef CXX_TEST
////////////////////////////////////////////////////////////////////////////////

#include <xmsmesh/meshing/detail/MePolyCleaner.t.h>

#include <xmscore/testing/TestTools.h>

// namespace xms {
using namespace xms;

////////////////////////////////////////////////////////////////////////////////
/// \class MePolyCleanerUnitTests
/// \brief tester for the MePolyCleaner class
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief tests creating the class
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCreateClass()
{
  BSHP<MePolyCleaner> b = MePolyCleaner::New();
  TS_ASSERT(b);
} // MePolyCleanerUnitTests::testCreateClass
//------------------------------------------------------------------------------
/// \brief simple triangle SANITY CHECK
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCase0()
{
  // x =   0         1
  //
  // y=1              1
  //               /|
  //             /  |
  //           /    |
  //        /       |
  // y=0   0---------2
  MePolyCleanerImpl p;
  std::vector<Pt3d> input = {{0, 0, 0}, {1, 1, 0}, {1, 0, 0}};
  MePolyOffsetterOutput out;
  p.CleanPolyOffset(input, 0, 1e-9, out);
  TS_ASSERT_EQUALS(1, out.m_loops.size());
  if (out.m_loops.size() != 1)
    return;
  std::vector<size_t> baseLoop = {0, 1, 2};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[0]);
  TS_ASSERT_EQUALS(0, out.m_loopTypes[0]);
  TS_ASSERT_EQUALS_VEC(input, out.m_pts);
} // MePolyCleanerUnitTests::testCase0
//------------------------------------------------------------------------------
/// \brief simple square SANITY CHECK
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCase1()
{
  // x =   0         1
  //
  // y=1   1---------2
  //      |         |
  //      |         |
  //      |         |
  //      |         |
  // y=0   0---------3
  MePolyCleanerImpl p;
  std::vector<Pt3d> input = {{0, 0, 0}, {0, 1, 0}, {1, 1, 0}, {1, 0, 0}};
  MePolyOffsetterOutput out;
  p.CleanPolyOffset(input, 0, 1e-9, out);
  TS_ASSERT_EQUALS(1, out.m_loops.size());
  if (out.m_loops.size() != 1)
    return;
  std::vector<size_t> baseLoop = {0, 1, 2, 3};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[0]);
  TS_ASSERT_EQUALS(0, out.m_loopTypes[0]);
  TS_ASSERT_EQUALS_VEC(input, out.m_pts);
} // MePolyCleanerUnitTests::testCase1
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCase2()
{
  // x =   0         1
  //
  // y=2   2---------3
  //        \       |
  //          \     |
  //             \  |
  // y=1            \|1
  //               /|
  //             /  |
  //           /    |
  //        /       |
  // y=0   0---------4

  MePolyCleanerImpl p;
  std::vector<Pt3d> input = {{0, 0, 0}, {1, 1, 0}, {0, 2, 0}, {1, 2, 0}, {1, 0, 0}};
  MePolyOffsetterOutput out;
  p.CleanPolyOffset(input, 0, 1e-9, out);
  TS_ASSERT_EQUALS(2, out.m_loops.size());
  if (out.m_loops.size() != 2)
    return;
  std::vector<Pt3d> base = input;
  TS_ASSERT_EQUALS_VEC(base, out.m_pts);
  std::vector<size_t> baseLoop = {1, 2, 3};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[0]);
  baseLoop = {0, 1, 4};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[1]);
  std::vector<int> baseType = {0, 0};
  TS_ASSERT_EQUALS_VEC(baseType, out.m_loopTypes);
} // MePolyCleanerUnitTests::testCase2
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCase2a()
{
  // x =   0         1
  //
  // y=2   0---------1
  //        \       |
  //          \     |
  //             \  |
  // y=1            \|4
  //               /|
  //             /  |
  //           /    |
  //        /       |
  // y=0   3---------2

  MePolyCleanerImpl p;
  std::vector<Pt3d> input = {{0, 2, 0}, {1, 2, 0}, {1, 0, 0}, {0, 0, 0}, {1, 1, 0}};
  MePolyOffsetterOutput out;
  p.CleanPolyOffset(input, 0, 1e-9, out);
  TS_ASSERT_EQUALS(2, out.m_loops.size());
  if (out.m_loops.size() != 2)
    return;
  std::vector<Pt3d> base = input;
  TS_ASSERT_EQUALS_VEC(base, out.m_pts);
  std::vector<size_t> baseLoop = {4, 2, 3};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[0]);
  baseLoop = {0, 1, 4};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[1]);
  std::vector<int> baseType = {0, 0};
  TS_ASSERT_EQUALS_VEC(baseType, out.m_loopTypes);
} // MePolyCleanerUnitTests::testCase2a
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCase3()
{
  // x =   0         2
  //
  // y=6   2---------3
  //        \       |
  //          \     |
  //            \   |
  //              \ |
  //               \|
  //                |\
  //                |  \
  //y=3             |    1
  //                |  /
  //                |/
  //               /|
  //             /  |
  //           /    |
  //         /      |
  //       /        |
  // y=0   0---------4
  MePolyCleanerImpl p;
  std::vector<Pt3d> input = {{0, 0, 0}, {3, 3, 0}, {0, 6, 0}, {2, 6, 0}, {2, 0, 0}};
  MePolyOffsetterOutput out;
  p.CleanPolyOffset(input, 0, 1e-9, out);
  TS_ASSERT_EQUALS(2, out.m_loops.size());
  if (out.m_loops.size() != 2)
    return;
  std::vector<Pt3d> base = input;
  base.push_back(Pt3d(2, 2, 0));
  base.push_back(Pt3d(2, 4, 0));
  TS_ASSERT_EQUALS_VEC(base, out.m_pts);
  std::vector<size_t> baseLoop = {6, 2, 3};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[0]);
  baseLoop = {0, 5, 4};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[1]);
  std::vector<int> baseType = {0, 0};
  TS_ASSERT_EQUALS_VEC(baseType, out.m_loopTypes);
} // MePolyCleanerUnitTests::testCase3
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCase3a()
{
  // x =   0         2
  //
  // y=6   4---------0
  //        \       |
  //          \     |
  //            \   |
  //              \ |
  //               \|
  //                |\
  //                |  \
  //y=3             |    3
  //                |  /
  //                |/
  //               /|
  //             /  |
  //           /    |
  //         /      |
  //       /        |
  // y=0   2---------1
  MePolyCleanerImpl p;
  std::vector<Pt3d> input = {{2, 6, 0}, {2, 0, 0}, {0, 0, 0}, {3, 3, 0}, {0, 6, 0}};
  MePolyOffsetterOutput out;
  p.CleanPolyOffset(input, 0, 1e-9, out);
  TS_ASSERT_EQUALS(2, out.m_loops.size());
  if (out.m_loops.size() != 2)
    return;
  std::vector<Pt3d> base = input;
  base.push_back(Pt3d(2, 2, 0));
  base.push_back(Pt3d(2, 4, 0));
  TS_ASSERT_EQUALS_VEC(base, out.m_pts);
  std::vector<size_t> baseLoop = {5, 1, 2};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[0]);
  baseLoop = {0, 6, 4};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[1]);
  std::vector<int> baseType = {0, 0};
  TS_ASSERT_EQUALS_VEC(baseType, out.m_loopTypes);
} // MePolyCleanerUnitTests::testCase3a
//------------------------------------------------------------------------------
/// \brief triangle with duplicate point
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCase4()
{
  // x =   0         1
  //
  // y=1               3
  //                /|
  //              /  |
  //            /    |
  //         /       |
  // y=0 (2)1---------0
  MePolyCleanerImpl p;
  std::vector<Pt3d> input = {{1, 0, 0}, {0, 0, 0}, {0, 0, 0}, {1, 1, 0}};
  MePolyOffsetterOutput out;
  p.CleanPolyOffset(input, 0, 1e-9, out);
  TS_ASSERT_EQUALS(1, out.m_loops.size());
  if (out.m_loops.size() != 1)
    return;
  std::vector<Pt3d> base = input;
  TS_ASSERT_EQUALS_VEC(base, out.m_pts);
  std::vector<size_t> baseLoop = {0, 1, 3};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[0]);
  std::vector<int> baseType = {0};
  TS_ASSERT_EQUALS_VEC(baseType, out.m_loopTypes);
} // MePolyCleanerUnitTests::testCase4
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCase5()
{
  // x =   0                4
  //
  // y=6   5----------------0
  //        \             /
  //         \           /
  //          \         /
  //           \       /
  //            \     /
  //             \   /
  //              \ / 1
  // y=3         4 / \
  //             /   \
  //            /     \
  //           /       \
  //          /         \
  //         /           \
  //        /             \
  //       /               \
  //y=0   3-----------------2
  MePolyCleanerImpl p;
  std::vector<Pt3d> input = {{4, 6, 0}, {2, 3, 0}, {4, 0, 0}, {0, 0, 0}, {2, 3, 0}, {0, 6, 0}};
  MePolyOffsetterOutput out;
  p.CleanPolyOffset(input, 0, 1e-9, out);
  TS_ASSERT_EQUALS(2, out.m_loops.size());
  if (out.m_loops.size() != 2)
    return;
  std::vector<Pt3d> base = input;
  TS_ASSERT_EQUALS_VEC(base, out.m_pts);
  std::vector<size_t> baseLoop = {1, 2, 3};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[0]);
  baseLoop = {0, 1, 5};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[1]);
  std::vector<int> baseType = {0, 0};
  TS_ASSERT_EQUALS_VEC(baseType, out.m_loopTypes);
} // MePolyCleanerUnitTests::testCase5
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCase5a()
{
  // x =   0                4
  //
  // y=6   1----------------2
  //        \             /
  //         \           /
  //          \         /
  //           \       /
  //            \     /
  //             \   /
  //              \ / 3
  // y=3         0 / \
  //             /   \
  //            /     \
  //           /       \
  //          /         \
  //         /           \
  //        /             \
  //       /               \
  //y=0   5-----------------4
  MePolyCleanerImpl p;
  std::vector<Pt3d> input = {{4, 2, 0}, {0, 6, 0}, {4, 6, 0}, {4, 2, 0}, {4, 0, 0}, {0, 0, 0}};
  MePolyOffsetterOutput out;
  p.CleanPolyOffset(input, 0, 1e-9, out);
  TS_ASSERT_EQUALS(2, out.m_loops.size());
  if (out.m_loops.size() != 2)
    return;
  std::vector<Pt3d> base = input;
  TS_ASSERT_EQUALS_VEC(base, out.m_pts);
  std::vector<size_t> baseLoop = {0, 1, 2};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[0]);
  baseLoop = {0, 4, 5};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[1]);
  std::vector<int> baseType = {0, 0};
  TS_ASSERT_EQUALS_VEC(baseType, out.m_loopTypes);
} // MePolyCleanerUnitTests::testCase5a
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCase6()
{
  // x =   0                4
  //
  // y=6   1----------------2
  //      |                |
  //      |                |
  //      |                |
  //      |                |
  // y=4   4----------------3
  //      |
  //      |
  //      |
  //      |
  //      |
  // y=2   5-----------------6
  //      |                 |
  //      |                 |
  //      |                 |
  //      |                 |
  // y=0   0-----------------7
  MePolyCleanerImpl p;
  std::vector<Pt3d> input = {{0, 0, 0}, {0, 6, 0}, {4, 6, 0}, {4, 4, 0},
                             {0, 4, 0}, {0, 2, 0}, {4, 2, 0}, {4, 0, 0}};
  MePolyOffsetterOutput out;
  p.CleanPolyOffset(input, 0, 1e-9, out);
  TS_ASSERT_EQUALS(2, out.m_loops.size());
  if (out.m_loops.size() != 2)
    return;
  std::vector<Pt3d> base = input;
  TS_ASSERT_EQUALS_VEC(base, out.m_pts);
  std::vector<size_t> baseLoop = {4, 1, 2, 3};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[0]);
  baseLoop = {0, 5, 6, 7};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[1]);
  std::vector<int> baseType = {0, 0};
  TS_ASSERT_EQUALS_VEC(baseType, out.m_loopTypes);
} // MePolyCleanerUnitTests::testCase6
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCase6a()
{
  // x =   0                4
  //
  // y=6   6----------------7
  //      |                |
  //      |                |
  //      |                |
  //      |                |
  // y=4   1----------------0
  //      |
  //      |
  //      |
  //      |
  //      |
  // y=2   2-----------------3
  //      |                 |
  //      |                 |
  //      |                 |
  //      |                 |
  // y=0   5-----------------4
  MePolyCleanerImpl p;
  std::vector<Pt3d> input = {{2, 4, 0}, {0, 4, 0}, {0, 2, 0}, {4, 2, 0},
                             {4, 0, 0}, {0, 0, 0}, {0, 6, 0}, {4, 6, 0}};
  MePolyOffsetterOutput out;
  p.CleanPolyOffset(input, 0, 1e-9, out);
  TS_ASSERT_EQUALS(2, out.m_loops.size());
  if (out.m_loops.size() != 2)
    return;
  std::vector<Pt3d> base = input;
  TS_ASSERT_EQUALS_VEC(base, out.m_pts);
  std::vector<size_t> baseLoop = {2, 3, 4, 5};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[0]);
  baseLoop = {0, 1, 6, 7};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[1]);
  std::vector<int> baseType = {0, 0};
  TS_ASSERT_EQUALS_VEC(baseType, out.m_loopTypes);
} // MePolyCleanerUnitTests::testCase6a
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCase6b()
{
  // x =   0                4
  //
  // y=6   7----------------8
  //      |                |
  //      |                |
  //      |                |
  //      |                |
  // y=4   1----------------0
  //      |
  //      |
  //      6
  //      |
  //      |
  // y=2   2-----------------3
  //      |                 |
  //      |                 |
  //      |                 |
  //      |                 |
  // y=0   5-----------------4
  MePolyCleanerImpl p;
  std::vector<Pt3d> input = {{4, 4, 0}, {0, 4, 0}, {0, 2, 0}, {4, 2, 0}, {4, 0, 0},
                             {0, 0, 0}, {0, 3, 0}, {0, 6, 0}, {4, 6, 0}};
  MePolyOffsetterOutput out;
  p.CleanPolyOffset(input, 0, 1e-9, out);
  TS_ASSERT_EQUALS(2, out.m_loops.size());
  if (out.m_loops.size() != 2)
    return;
  std::vector<Pt3d> base = input;
  TS_ASSERT_EQUALS_VEC(base, out.m_pts);
  std::vector<size_t> baseLoop = {2, 3, 4, 5};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[0]);
  baseLoop = {0, 1, 7, 8};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[1]);
  std::vector<int> baseType = {0, 0};
  TS_ASSERT_EQUALS_VEC(baseType, out.m_loopTypes);
} // MePolyCleanerUnitTests::testCase6b
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCase6c()
{
  // x =     0                4
  //
  // y=6     9----------------10
  //        |                |
  //        |                |
  //        |                |
  //        |                |
  // y=4  (8)1----------------0
  //        |
  //        |
  //        7
  //        |
  //        |
  // y=2  (6)2-----------------3
  //        |                 |
  //        |                 |
  //        |                 |
  //        |                 |
  // y=0     5-----------------4
  MePolyCleanerImpl p;
  std::vector<Pt3d> input = {{4, 4, 0}, {0, 4, 0}, {0, 2, 0}, {4, 2, 0}, {4, 0, 0}, {0, 0, 0},
                             {0, 2, 0}, {0, 3, 0}, {0, 4, 0}, {0, 6, 0}, {4, 6, 0}};
  MePolyOffsetterOutput out;
  p.CleanPolyOffset(input, 0, 1e-9, out);
  TS_ASSERT_EQUALS(2, out.m_loops.size());
  if (out.m_loops.size() != 2)
    return;
  std::vector<Pt3d> base = input;
  TS_ASSERT_EQUALS_VEC(base, out.m_pts);
  std::vector<size_t> baseLoop = {2, 3, 4, 5};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[0]);
  baseLoop = {0, 1, 9, 10};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[1]);
  std::vector<int> baseType = {0, 0};
  TS_ASSERT_EQUALS_VEC(baseType, out.m_loopTypes);
} // MePolyCleanerUnitTests::testCase6c
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCase7()
{
  // x =     0                4                8
  //
  // y=6    10---------------------------------11
  //        |                                 |
  //        |                                 |
  //        |                                 |
  //        |                1(8)             |
  // y=4     9---------------------------------0
  //                         |
  //                         |
  //                         |
  //                         |
  //                         |
  // y=2     6----------------------------------3
  //        |                2(7)              |
  //        |                                  |
  //        |                                  |
  //        |                                  |
  // y=0     5----------------------------------4
  MePolyCleanerImpl p;
  std::vector<Pt3d> input = {{8, 4, 0}, {4, 4, 0}, {4, 2, 0}, {8, 2, 0}, {8, 0, 0}, {0, 0, 0},
                             {0, 2, 0}, {4, 2, 0}, {4, 4, 0}, {0, 4, 0}, {0, 6, 0}, {8, 6, 0}};
  MePolyOffsetterOutput out;
  p.CleanPolyOffset(input, 0, 1e-9, out);
  TS_ASSERT_EQUALS(2, out.m_loops.size());
  if (out.m_loops.size() != 2)
    return;
  std::vector<Pt3d> base = input;
  TS_ASSERT_EQUALS_VEC(base, out.m_pts);
  std::vector<size_t> baseLoop = {2, 3, 4, 5, 6};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[0]);
  baseLoop = {0, 1, 9, 10, 11};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[1]);
  std::vector<int> baseType = {0, 0};
  TS_ASSERT_EQUALS_VEC(baseType, out.m_loopTypes);
} // MePolyCleanerUnitTests::testCase7
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCase7a()
{
  // x =     0                4                8
  //
  // y=6    10---------------------------------11
  //        |                                 |
  //        |                                 |
  //        9----------------8                |
  //                         |                |
  // y=4                      -----------------0
  //                         | 1
  //                         |
  //                         |
  //                         |
  //                         | 2
  // y=2                      ------------------3
  //                         |                 |
  //        6----------------7                 |
  //        |                                  |
  //        |                                  |
  // y=0     5----------------------------------4
  MePolyCleanerImpl p;
  std::vector<Pt3d> input = {{8, 4, 0}, {4, 4, 0}, {4, 2, 0}, {8, 2, 0}, {8, 0, 0}, {0, 0, 0},
                             {0, 1, 0}, {4, 1, 0}, {4, 5, 0}, {0, 5, 0}, {0, 6, 0}, {8, 6, 0}};
  MePolyOffsetterOutput out;
  p.CleanPolyOffset(input, 0, 1e-9, out);
  TS_ASSERT_EQUALS(2, out.m_loops.size());
  if (out.m_loops.size() != 2)
    return;
  std::vector<Pt3d> base = input;
  TS_ASSERT_EQUALS_VEC(base, out.m_pts);
  std::vector<size_t> baseLoop = {2, 3, 4, 5, 6, 7};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[0]);
  baseLoop = {0, 1, 8, 9, 10, 11};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[1]);
  std::vector<int> baseType = {0, 0};
  TS_ASSERT_EQUALS_VEC(baseType, out.m_loopTypes);
} // MePolyCleanerUnitTests::testCase7a
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCase7b()
{
  // x =     0                4                8
  //
  // y=6    10---------------------------------11
  //        |                                 |
  //        |                                 |
  //        9----------------8                |
  //                         |                |
  // y=4                      -----------------0
  //                         | 1
  //                         |
  //        6----------------| 7
  //        |                |
  //        |                | 2
  // y=2     |                ------------------3
  //        |                                  |
  //        |                                  |
  //        |                                  |
  //        |                                  |
  // y=0     5----------------------------------4
  MePolyCleanerImpl p;
  std::vector<Pt3d> input = {{8, 4, 0}, {4, 4, 0}, {4, 2, 0}, {8, 2, 0}, {8, 0, 0}, {0, 0, 0},
                             {0, 3, 0}, {4, 3, 0}, {4, 5, 0}, {0, 5, 0}, {0, 6, 0}, {8, 6, 0}};
  MePolyOffsetterOutput out;
  p.CleanPolyOffset(input, 0, 1e-9, out);
  TS_ASSERT_EQUALS(2, out.m_loops.size());
  if (out.m_loops.size() != 2)
    return;
  std::vector<Pt3d> base = input;
  TS_ASSERT_EQUALS_VEC(base, out.m_pts);
  std::vector<size_t> baseLoop = {7, 2, 3, 4, 5, 6};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[0]);
  baseLoop = {0, 1, 8, 9, 10, 11};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[1]);
  std::vector<int> baseType = {0, 0};
  TS_ASSERT_EQUALS_VEC(baseType, out.m_loopTypes);
} // MePolyCleanerUnitTests::testCase7b
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCase7c()
{
  // x =     0                4                8
  //
  // y=6    10---------------------------------11
  //        |                                 |
  //        |                                 |
  //        |                                 |
  //        |                1                |
  // y=4     |                -----------------0
  //        |                |
  //        9----------------| 8
  //                         |
  //        6----------------| 7
  //        |                |
  // y=2     |                ------------------3
  //        |                2                 |
  //        |                                  |
  //        |                                  |
  //        |                                  |
  // y=0     5----------------------------------4
  MePolyCleanerImpl p;
  std::vector<Pt3d> input = {{8, 4, 0},   {4, 4, 0},   {4, 2, 0},   {8, 2, 0},
                             {8, 0, 0},   {0, 0, 0},   {0, 2.5, 0}, {4, 2.5, 0},
                             {4, 3.5, 0}, {0, 3.5, 0}, {0, 6, 0},   {8, 6, 0}};
  MePolyOffsetterOutput out;
  p.CleanPolyOffset(input, 0, 1e-9, out);
  TS_ASSERT_EQUALS(2, out.m_loops.size());
  if (out.m_loops.size() != 2)
    return;
  std::vector<Pt3d> base = input;
  TS_ASSERT_EQUALS_VEC(base, out.m_pts);
  std::vector<size_t> baseLoop = {7, 2, 3, 4, 5, 6};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[0]);
  baseLoop = {0, 1, 8, 9, 10, 11};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[1]);
  std::vector<int> baseType = {0, 0};
  TS_ASSERT_EQUALS_VEC(baseType, out.m_loopTypes);
} // MePolyCleanerUnitTests::testCase7c
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCase8()
{
  // x =     0                4                 8
  //
  // y=6    10----------------------------------9
  //        |                                  |
  //        |                                  |
  //        |                                  |
  //        |                                  |
  // y=4     |                 6        3       8
  //        |               /   \     / \     /
  //        |             /       \  /   \   /
  //        |           5          /\      \/
  //           |             \       /  \    / \
//        |               \   /     \ /     \
//y=2     |                 4        7       2
  //        |                                  |
  //        |                                  |
  //        |                                  |
  //        |                                  |
  // y=0     0----------------------------------1
  MePolyCleanerImpl p;
  std::vector<Pt3d> input = {{0, 0, 0}, {8, 0, 0}, {8, 2, 0}, {6, 4, 0}, {4, 2, 0}, {3, 3, 0},
                             {4, 4, 0}, {6, 2, 0}, {8, 4, 0}, {8, 6, 0}, {0, 6, 0}};
  MePolyOffsetterOutput out;
  p.CleanPolyOffset(input, 1, 1e-9, out);
  TS_ASSERT_EQUALS(2, out.m_loops.size());
  if (out.m_loops.size() != 2)
    return;
  std::vector<Pt3d> base = input;
  base.push_back(Pt3d(7, 3, 0));
  base.push_back(Pt3d(5, 3, 0));
  TS_ASSERT_EQUALS_VEC(base, out.m_pts);
  std::vector<size_t> baseLoop = {12, 4, 5, 6};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[0]);
  baseLoop = {0, 1, 2, 11, 8, 9, 10};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[1]);
  std::vector<int> baseType = {2, 1};
  TS_ASSERT_EQUALS_VEC(baseType, out.m_loopTypes);
} // MePolyCleanerUnitTests::testCase8
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCase8a()
{
  // x =     0                4                 8
  //
  // y=6     1----------------------------------0
  //        |                                  |
  //        |                                  |
  //        |                                  |
  //        |                                  |
  // y=4     |                 8        5       10
  //        |               /   \     / \     /
  //        |             /       \  /   \   /
  //        |           7          /\      \/
  //        |             \       /  \    / \
  //        |               \   /     \ /     \
  //y=2     |                 6        9       4
  //        |                                  |
  //        |                                  |
  //        |                                  |
  //        |                                  |
  // y=0     2----------------------------------3
  MePolyCleanerImpl p;
  std::vector<Pt3d> input = {{8, 6, 0}, {0, 6, 0}, {0, 0, 0}, {8, 0, 0}, {8, 2, 0}, {6, 4, 0},
                             {4, 2, 0}, {3, 3, 0}, {4, 4, 0}, {6, 2, 0}, {8, 4, 0}};
  MePolyOffsetterOutput out;
  p.CleanPolyOffset(input, 1, 1e-9, out);
  TS_ASSERT_EQUALS(2, out.m_loops.size());
  if (out.m_loops.size() != 2)
    return;
  std::vector<Pt3d> base = input;
  base.push_back(Pt3d(7, 3, 0));
  base.push_back(Pt3d(5, 3, 0));
  TS_ASSERT_EQUALS_VEC(base, out.m_pts);
  std::vector<size_t> baseLoop = {12, 6, 7, 8};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[0]);
  baseLoop = {0, 1, 2, 3, 4, 11, 10};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[1]);
  std::vector<int> baseType = {2, 1};
  TS_ASSERT_EQUALS_VEC(baseType, out.m_loopTypes);
} // MePolyCleanerUnitTests::testCase8a
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCase8b()
{
  // x =     0                4                 8
  //
  // y=6     7----------------------------------6
  //        |                                  |
  //        |                                  |
  //        |                                  |
  //        |                                  |
  // y=4     |                 3        0       5
  //        |               /   \     / \     /
  //        |             /       \  /   \   /
  //        |           2          /\      \/
  //        |             \       /  \    / \
  //        |               \   /     \ /     \
  //y=2     |                 1        4       10
  //        |                                  |
  //        |                                  |
  //        |                                  |
  //        |                                  |
  // y=0     8----------------------------------9
  MePolyCleanerImpl p;
  std::vector<Pt3d> input = {{6, 4, 0}, {4, 2, 0}, {3, 3, 0}, {4, 4, 0}, {6, 2, 0}, {8, 4, 0},
                             {8, 6, 0}, {0, 6, 0}, {0, 0, 0}, {6, 0, 0}, {8, 2, 0}};
  MePolyOffsetterOutput out;
  p.CleanPolyOffset(input, 1, 1e-9, out);
  TS_ASSERT_EQUALS(2, out.m_loops.size());
  if (out.m_loops.size() != 2)
    return;
  std::vector<Pt3d> base = input;
  base.push_back(Pt3d(5, 3, 0));
  base.push_back(Pt3d(7, 3, 0));
  TS_ASSERT_EQUALS_VEC(base, out.m_pts);
  std::vector<size_t> baseLoop = {11, 1, 2, 3};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[0]);
  baseLoop = {12, 5, 6, 7, 8, 9, 10};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[1]);
  std::vector<int> baseType = {2, 1};
  TS_ASSERT_EQUALS_VEC(baseType, out.m_loopTypes);
} // MePolyCleanerUnitTests::testCase8b
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCase8c()
{
  // x =     0                4                 8
  //
  // y=6     7----------------------------------6
  //        |                                  |
  //        |                                  |
  //        |                                  |
  //        |                                  |
  // y=4     |                 3        0       5
  //        |               /   \     /|      /
  //        |             /       \  / |     /
  //        |           2          /\  |    /
  //        |             \       /  \ |  /
  //        |               \   /     \|/
  // y=2     |                 1        4(10)
  //        |                            \
  //        |                              \
  //        |                                \
  //        |                                 \
  //y=0     8----------------------------------9
  MePolyCleanerImpl p;
  std::vector<Pt3d> input = {{6, 4, 0}, {4, 2, 0}, {3, 3, 0}, {4, 4, 0}, {6, 2, 0}, {8, 4, 0},
                             {8, 6, 0}, {0, 6, 0}, {0, 0, 0}, {6, 0, 0}, {6, 2, 0}};
  MePolyOffsetterOutput out;
  p.CleanPolyOffset(input, 1, 1e-9, out);
  TS_ASSERT_EQUALS(2, out.m_loops.size());
  if (2 != out.m_loops.size())
    return;
  std::vector<Pt3d> basePts = {{6, 4, 0}, {4, 2, 0}, {3, 3, 0}, {4, 4, 0}, {6, 2, 0}, {8, 4, 0},
                               {8, 6, 0}, {0, 6, 0}, {0, 0, 0}, {6, 0, 0}, {6, 2, 0}, {5, 3, 0}};
  TS_ASSERT_EQUALS_VEC(basePts, out.m_pts);
  std::vector<size_t> baseLoop = {11, 1, 2, 3};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[0]);
  baseLoop = {4, 5, 6, 7, 8, 9};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[1]);
  std::vector<int> basePolyType = {2, 1};
  TS_ASSERT_EQUALS_VEC(basePolyType, out.m_loopTypes);
} // MePolyCleanerUnitTests::testCase8c
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCase9()
{
  // x =     0                4                 8
  //
  // y=6     1----------------------------------0
  //        |                                  |
  //        |                                  |
  //        |                                  |
  //        |                                  |
  // y=4     |                                  |
  //        |                                  |
  //        |               5                  |
  //        |               /\                 |
  //        |              /  \                |
  //        |             /    \               |
  // y=2     |          6 /      \4             |
  //        |            \      /              |
  //        |             \    /               |
  //        |              \  /                |
  //        |              7\/                 |
  // y=0     2----------------3-----------------8
  MePolyCleanerImpl p;
  std::vector<Pt3d> input = {{8, 6, 0}, {0, 6, 0}, {0, 0, 0}, {4, 0, 0}, {5, 2, 0},
                             {3, 4, 0}, {3, 2, 0}, {4, 0, 0}, {8, 0, 0}};
  MePolyOffsetterOutput out;
  p.CleanPolyOffset(input, 1, 1e-9, out);
  TS_ASSERT_EQUALS(1, out.m_loops.size());
  if (1 != out.m_loops.size())
    return;
  std::vector<Pt3d> basePts = {{8, 6, 0}, {0, 6, 0}, {0, 0, 0}, {4, 0, 0}, {5, 2, 0},
                               {3, 4, 0}, {3, 2, 0}, {4, 0, 0}, {8, 0, 0}};
  TS_ASSERT_EQUALS_VEC(basePts, out.m_pts);
  std::vector<size_t> baseLoop = {0, 1, 2, 3, 8};
  TS_ASSERT_EQUALS_VEC(baseLoop, out.m_loops[0]);
  std::vector<int> basePolyType = {1};
  TS_ASSERT_EQUALS_VEC(basePolyType, out.m_loopTypes);
} // MePolyCleanerUnitTests::testCase9
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCleanIn0()
{
  // x =     0                4                 8
  //
  //
  // y=2     3----------2     7--------6
  //        |          |     |        |
  //        |          |     |        |
  //        |          |     |        |
  //        |          |     |        |
  // y=0     0----------1     4--------5
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{0, 0, 0}, {3, 0, 0}, {3, 2, 0}, {0, 2, 0},
             {4, 0, 0}, {6, 0, 0}, {6, 2, 0}, {4, 2, 0}};
  std::vector<size_t> v0 = {0, 1, 2, 3};
  o.m_loops.push_back(v0);
  std::vector<size_t> v1 = {4, 5, 6, 7};
  o.m_loops.push_back(v1);
  o.m_loopTypes.assign(2, MePolyOffsetter::INSIDE_POLY);
  MePolyCleanerImpl p;
  std::vector<MePolyOffsetterOutput> vO(1, o);
  p.IntersectCleanInPolys(vO, o1, 1e-9);
  TS_ASSERT_EQUALS(2, o1.m_loops.size());
  if (2 != o1.m_loops.size())
    return;
  TS_ASSERT_EQUALS_VEC(o.m_pts, o1.m_pts);
  TS_ASSERT_EQUALS_VEC(o.m_loops[0], o1.m_loops[0]);
  TS_ASSERT_EQUALS_VEC(o.m_loops[1], o1.m_loops[1]);
  TS_ASSERT_EQUALS_VEC(o.m_loopTypes, o1.m_loopTypes);
} // MePolyCleanerUnitTests::testCleanIn0
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCleanIn1()
{
  // x =     0                4                 8
  //
  // y=4     3----------------------------------2
  //        |                                  |
  //        |                                  |
  //        |                                  |
  //        |                                  |
  //        |                                  |
  // y=2     |        7------6                  |
  //        |        |      |                  |
  //        |        |      |                  |
  //        |        4------5                  |
  //        |                                  |
  // y=0     0----------------------------------1
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{0, 0, 0}, {8, 0, 0}, {8, 4, 0}, {0, 4, 0},
             {2, 1, 0}, {4, 1, 0}, {4, 2, 0}, {2, 2, 0}};
  std::vector<size_t> v0 = {0, 1, 2, 3};
  o.m_loops.push_back(v0);
  std::vector<size_t> v1 = {4, 5, 6, 7};
  o.m_loops.push_back(v1);
  o.m_loopTypes.assign(2, MePolyOffsetter::INSIDE_POLY);
  MePolyCleanerImpl p;
  std::vector<MePolyOffsetterOutput> vO(1, o);
  p.IntersectCleanInPolys(vO, o1, 1e-9);
  TS_ASSERT_EQUALS(1, o1.m_loops.size());
  if (1 != o1.m_loops.size())
    return;
  TS_ASSERT_EQUALS_VEC(o.m_pts, o1.m_pts);
  TS_ASSERT_EQUALS_VEC(o.m_loops[0], o1.m_loops[0]);
  std::vector<int> baseLoopType = {1};
  TS_ASSERT_EQUALS_VEC(baseLoopType, o1.m_loopTypes);
} // MePolyCleanerUnitTests::testCleanIn1
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCleanIn1a()
{
  // x =     0                4                 8
  //
  // y=4     7----------------------------------6
  //        |                                  |
  //        |                                  |
  //        |                                  |
  //        |                                  |
  //        |                                  |
  // y=2     |        3------2                  |
  //        |        |      |                  |
  //        |        |      |                  |
  //        |        0------1                  |
  //        |                                  |
  // y=0     4----------------------------------5
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{2, 1, 0}, {4, 1, 0}, {4, 2, 0}, {2, 2, 0},
             {0, 0, 0}, {8, 0, 0}, {8, 4, 0}, {0, 4, 0}};
  std::vector<size_t> v0 = {0, 1, 2, 3};
  o.m_loops.push_back(v0);
  std::vector<size_t> v1 = {4, 5, 6, 7};
  o.m_loops.push_back(v1);
  o.m_loopTypes.assign(2, MePolyOffsetter::INSIDE_POLY);
  MePolyCleanerImpl p;
  std::vector<MePolyOffsetterOutput> vO(1, o);
  p.IntersectCleanInPolys(vO, o1, 1e-9);
  TS_ASSERT_EQUALS(1, o1.m_loops.size());
  if (1 != o1.m_loops.size())
    return;
  TS_ASSERT_EQUALS_VEC(o.m_pts, o1.m_pts);
  TS_ASSERT_EQUALS_VEC(o.m_loops[1], o1.m_loops[0]);
  std::vector<int> baseLoopType = {1};
  TS_ASSERT_EQUALS_VEC(baseLoopType, o1.m_loopTypes);
} // MePolyCleanerUnitTests::testCleanIn1a
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCleanIn1b()
{
  // x =     0                4                 8
  //
  // y=4     7----------------------------------6
  //        |                                  |
  //        |                    15------14    |
  //        |                     |      |     |
  //        |                    12------13    |
  //        |                                  |
  // y=2     |        3------2     8------11    |
  //        |        |      |     |      |     |
  //        |        |      |     |      |     |
  //        |        0------1     9------10    |
  //        |                                  |
  // y=0     4----------------------------------5
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{2, 1, 0}, {4, 1, 0}, {4, 2, 0},   {2, 2, 0},  {0, 0, 0}, {8, 0, 0},
             {8, 4, 0}, {0, 4, 0}, {5, 1, 0},   {7, 1, 0},  {7, 2, 0}, {5, 2, 0},
             {5, 3, 0}, {7, 3, 0}, {7, 3.5, 0}, {5, 3.5, 0}};
  std::vector<size_t> v0 = {0, 1, 2, 3};
  o.m_loops.push_back(v0);
  std::vector<size_t> v1 = {4, 5, 6, 7};
  o.m_loops.push_back(v1);
  std::vector<size_t> v2 = {8, 9, 10, 11};
  o.m_loops.push_back(v2);
  std::vector<size_t> v3 = {12, 13, 14, 15};
  o.m_loops.push_back(v3);
  o.m_loopTypes.assign(4, MePolyOffsetter::INSIDE_POLY);
  MePolyCleanerImpl p;
  std::vector<MePolyOffsetterOutput> vO(1, o);
  p.IntersectCleanInPolys(vO, o1, 1e-9);
  TS_ASSERT_EQUALS(1, o1.m_loops.size());
  if (1 != o1.m_loops.size())
    return;
  TS_ASSERT_EQUALS_VEC(o.m_pts, o1.m_pts);
  TS_ASSERT_EQUALS_VEC(o.m_loops[1], o1.m_loops[0]);
  std::vector<int> baseLoopType = {1};
  TS_ASSERT_EQUALS_VEC(baseLoopType, o1.m_loopTypes);
} // MePolyCleanerUnitTests::testCleanIn1b
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCleanIn2()
{
  // x =     0                4                 8
  //
  // y=4     5----------------------------------4
  //        |                                  |
  //        |                                  |
  //        |        2-------------------------3
  //        |        |
  //        |        |
  // y=2     |        |       9------8
  //        |        |       |      |
  //        |        |       |      |
  //        |        |       6------7
  //        |        |
  // y=0     0--------1
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{0, 0, 0}, {2, 0, 0}, {2, 3, 0}, {8, 3, 0}, {8, 4, 0},
             {0, 4, 0}, {4, 1, 0}, {5, 1, 0}, {5, 2, 0}, {4, 2, 0}};
  std::vector<size_t> v0 = {0, 1, 2, 3, 4, 5};
  o.m_loops.push_back(v0);
  std::vector<size_t> v1 = {6, 7, 8, 9};
  o.m_loops.push_back(v1);
  o.m_loopTypes.assign(2, MePolyOffsetter::INSIDE_POLY);
  MePolyCleanerImpl p;
  std::vector<MePolyOffsetterOutput> vO(1, o);
  p.IntersectCleanInPolys(vO, o1, 1e-9);
  TS_ASSERT_EQUALS(2, o1.m_loops.size());
  if (2 != o1.m_loops.size())
    return;
  TS_ASSERT_EQUALS_VEC(o.m_pts, o1.m_pts);
  TS_ASSERT_EQUALS_VEC(o.m_loops[0], o1.m_loops[0]);
  TS_ASSERT_EQUALS_VEC(o.m_loops[1], o1.m_loops[1]);
  std::vector<int> baseLoopType = {1, 1};
  TS_ASSERT_EQUALS_VEC(baseLoopType, o1.m_loopTypes);
} // MePolyCleanerUnitTests::testCleanIn2
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCleanIn2a()
{
  // x =     0                4                 8
  //
  // y=4     5----------------4
  //        |                |
  //        |                |
  //        |        2-------3(9)---8
  //        |        |       |      |
  //        |        |       |      |
  // y=2     |        |       6------7
  //        |        |
  //        |        |
  //        |        |
  //        |        |
  // y=0     0--------1
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{0, 0, 0}, {2, 0, 0}, {2, 3, 0}, {4, 3, 0}, {4, 4, 0},
             {0, 4, 0}, {4, 2, 0}, {5, 2, 0}, {5, 3, 0}, {4, 3, 0}};
  std::vector<size_t> v0 = {0, 1, 2, 3, 4, 5};
  o.m_loops.push_back(v0);
  std::vector<size_t> v1 = {6, 7, 8, 9};
  o.m_loops.push_back(v1);
  o.m_loopTypes.assign(2, MePolyOffsetter::INSIDE_POLY);
  MePolyCleanerImpl p;
  std::vector<MePolyOffsetterOutput> vO(1, o);
  p.IntersectCleanInPolys(vO, o1, 1e-9);
  TS_ASSERT_EQUALS(2, o1.m_loops.size());
  if (2 != o1.m_loops.size())
    return;
  TS_ASSERT_EQUALS_VEC(o.m_pts, o1.m_pts);
  TS_ASSERT_EQUALS_VEC(o.m_loops[0], o1.m_loops[0]);
  std::vector<size_t> baseLoop = {6, 7, 8, 3};
  TS_ASSERT_EQUALS_VEC(baseLoop, o1.m_loops[1]);
  std::vector<int> baseLoopType = {1, 1};
  TS_ASSERT_EQUALS_VEC(baseLoopType, o1.m_loopTypes);
} // MePolyCleanerUnitTests::testCleanIn2a
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCleanIn2b()
{
  // x =     0                4                 8
  //
  // y=4     5----------------4
  //        |                |
  //        |         (3)    |
  //        |        2-------3(2)
  //           |        |\       \
//        |        | \       \
//y=2     |        |  \       \
//        |        |   \       \
//        |        |    \       \
//        |        |     0-------1
  //        |        |
  // y=0     0--------1
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{0, 0, 0}, {2, 0, 0}, {2, 3, 0}, {4, 3, 0}, {4, 4, 0}, {0, 4, 0}};
  std::vector<size_t> v0 = {0, 1, 2, 3, 4, 5};
  o.m_loops.push_back(v0);
  o.m_loopTypes.assign(1, MePolyOffsetter::INSIDE_POLY);
  std::vector<MePolyOffsetterOutput> vO(1, o);
  o.m_pts = {{3.5, 1, 0}, {4.5, 1, 0}, {4, 3, 0}, {2, 3, 0}};
  std::vector<size_t> v1 = {0, 1, 2, 3};
  o.m_loops[0] = v1;
  vO.push_back(o);
  MePolyCleanerImpl p;
  p.IntersectCleanInPolys(vO, o1, 1e-9);
  TS_ASSERT_EQUALS(1, o1.m_loops.size());
  if (2 != o1.m_loops.size())
    return;
  std::vector<Pt3d> basePts = {{0, 0, 0}, {2, 0, 0},   {2, 3, 0},   {4, 3, 0}, {4, 4, 0},
                               {0, 4, 0}, {3.5, 1, 0}, {4.5, 1, 0}, {4, 3, 0}, {2, 3, 0}};
  TS_ASSERT_EQUALS_VEC(basePts, o1.m_pts);
  std::vector<size_t> baseLoop = {0, 1, 2, 6, 7, 3, 4, 5};
  TS_ASSERT_EQUALS_VEC(baseLoop, o1.m_loops[0]);
  std::vector<int> baseLoopType = {1};
  TS_ASSERT_EQUALS_VEC(baseLoopType, o1.m_loopTypes);
} // MePolyCleanerUnitTests::testCleanIn2b
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCleanIn3()
{
  // x =     0                4                 8
  //
  // y=4     3---------------------2
  //        |                     |
  //        |                     |
  //        |        (7)3---------|------------2(6)
  //        |           |         |(8)         |
  //        |           |         |            |
  // y=2     |           |         |            |
  //        |           |         |            |
  //        |        (9)|         |            |
  //        0---------------------1            |
  //                    |                      |
  // y=0              (4)0----------------------1(5)
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{0, 1, 0}, {5, 1, 0}, {5, 4, 0}, {0, 4, 0}};
  std::vector<size_t> v0 = {0, 1, 2, 3};
  o.m_loops.push_back(v0);
  o.m_loopTypes.assign(1, MePolyOffsetter::INSIDE_POLY);
  std::vector<MePolyOffsetterOutput> vO(1, o);
  o.m_pts = {{3, 0, 0}, {8, 0, 0}, {8, 3, 0}, {3, 3, 0}};
  vO.push_back(o);
  MePolyCleanerImpl p;
  p.IntersectCleanInPolys(vO, o1, 1e-9);
  TS_ASSERT_EQUALS(1, o1.m_loops.size());
  if (1 != o1.m_loops.size())
    return;
  std::vector<Pt3d> basePts = {{0, 1, 0}, {5, 1, 0}, {5, 4, 0}, {0, 4, 0}, {3, 0, 0},
                               {8, 0, 0}, {8, 3, 0}, {3, 3, 0}, {5, 3, 0}, {3, 1, 0}};
  TS_ASSERT_EQUALS_VEC(basePts, o1.m_pts);
  std::vector<size_t> baseLoop = {9, 4, 5, 6, 8, 2, 3, 0};
  TS_ASSERT_EQUALS_VEC(baseLoop, o1.m_loops[0]);
  std::vector<int> baseLoopType = {1};
  TS_ASSERT_EQUALS_VEC(baseLoopType, o1.m_loopTypes);
} // MePolyCleanerUnitTests::testCleanIn3
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCleanIn3a()
{
  // x =     0                4                 8
  //
  // y=4  (4)0---------------------3(7)
  //        |                     |
  //        |                  (9)|
  //        |           0---------|------------3
  //        |           |         |            |
  //        |           |         |            |
  // y=2     |           |         |            |
  //        |           |         |            |
  //        |        (8)|         |            |
  //     (5)1---------------------2(6)         |
  //                    |                      |
  // y=0                 1----------------------2
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{3, 3, 0}, {3, 0, 0}, {8, 0, 0}, {8, 3, 0}};
  std::vector<size_t> v0 = {0, 1, 2, 3};
  o.m_loops.push_back(v0);
  o.m_loopTypes.assign(1, MePolyOffsetter::INSIDE_POLY);
  std::vector<MePolyOffsetterOutput> vO(1, o);
  o.m_pts = {{0, 4, 0}, {0, 1, 0}, {5, 1, 0}, {5, 4, 0}};
  vO.push_back(o);
  MePolyCleanerImpl p;
  p.IntersectCleanInPolys(vO, o1, 1e-9);
  TS_ASSERT_EQUALS(1, o1.m_loops.size());
  if (1 != o1.m_loops.size())
    return;
  std::vector<Pt3d> basePts = {{3, 3, 0}, {3, 0, 0}, {8, 0, 0}, {8, 3, 0}, {0, 4, 0},
                               {0, 1, 0}, {5, 1, 0}, {5, 4, 0}, {3, 1, 0}, {5, 3, 0}};
  TS_ASSERT_EQUALS_VEC(basePts, o1.m_pts);
  std::vector<size_t> baseLoop = {9, 7, 4, 5, 8, 1, 2, 3};
  TS_ASSERT_EQUALS_VEC(baseLoop, o1.m_loops[0]);
  std::vector<int> baseLoopType = {1};
  TS_ASSERT_EQUALS_VEC(baseLoopType, o1.m_loopTypes);
} // MePolyCleanerUnitTests::testCleanIn3a
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCleanIn4()
{
  // x =     0                4                 8
  //
  // y=4     0---------------------3
  //        |                     |
  //        |                     |
  //        |           0---------|------------3
  //        |           |         |            |
  //        |           |         |            |
  // y=2     |           |         |            |
  //        |           |         |            |
  //        |           |         |            |
  //        |           1---------|------------2
  //        |                     |
  // y=0     1---------------------2
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{3, 3, 0}, {3, 1, 0}, {8, 1, 0}, {8, 3, 0}};
  std::vector<size_t> v0 = {0, 1, 2, 3};
  o.m_loops.push_back(v0);
  o.m_loopTypes.assign(1, MePolyOffsetter::INSIDE_POLY);
  std::vector<MePolyOffsetterOutput> vO(1, o);
  o.m_pts = {{0, 4, 0}, {0, 0, 0}, {5, 0, 0}, {5, 4, 0}};
  vO.push_back(o);
  MePolyCleanerImpl p;
  p.IntersectCleanInPolys(vO, o1, 1e-9);
  TS_ASSERT_EQUALS(1, o1.m_loops.size());
  if (1 != o1.m_loops.size())
    return;
  std::vector<Pt3d> basePts = {{3, 3, 0}, {3, 1, 0}, {8, 1, 0}, {8, 3, 0}, {0, 4, 0},
                               {0, 0, 0}, {5, 0, 0}, {5, 4, 0}, {5, 1, 0}, {5, 3, 0}};
  TS_ASSERT_EQUALS_VEC(basePts, o1.m_pts);
  std::vector<size_t> baseLoop = {9, 7, 4, 5, 6, 8, 2, 3};
  TS_ASSERT_EQUALS_VEC(baseLoop, o1.m_loops[0]);
  std::vector<int> baseLoopType = {1};
  TS_ASSERT_EQUALS_VEC(baseLoopType, o1.m_loopTypes);
} // MePolyCleanerUnitTests::testCleanIn4
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCleanIn4a()
{
  // x =     0                4                 8
  //
  // y=4  (4)0---------------------3(7)
  //        |                     |
  //        |                     |
  //        |           0---------|------------3
  //        |           |         |            |
  //        |           |         |            |
  // y=2     |           |         |            |
  //        |           |         |            |
  //        |           |         |            |
  //        |           |         |            |
  //        |           |         |(6)         |
  // y=0  (5)1-----------1---------2------------2
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{3, 3, 0}, {3, 0, 0}, {8, 0, 0}, {8, 3, 0}};
  std::vector<size_t> v0 = {0, 1, 2, 3};
  o.m_loops.push_back(v0);
  o.m_loopTypes.assign(1, MePolyOffsetter::INSIDE_POLY);
  std::vector<MePolyOffsetterOutput> vO(1, o);
  o.m_pts = {{0, 4, 0}, {0, 0, 0}, {5, 0, 0}, {5, 4, 0}};
  vO.push_back(o);
  MePolyCleanerImpl p;
  p.IntersectCleanInPolys(vO, o1, 1e-9);
  TS_ASSERT_EQUALS(1, o1.m_loops.size());
  if (1 != o1.m_loops.size())
    return;
  std::vector<Pt3d> basePts = {{3, 3, 0}, {3, 0, 0}, {8, 0, 0}, {8, 3, 0}, {0, 4, 0},
                               {0, 0, 0}, {5, 0, 0}, {5, 4, 0}, {5, 3, 0}};
  TS_ASSERT_EQUALS_VEC(basePts, o1.m_pts);
  std::vector<size_t> baseLoop = {8, 7, 4, 5, 1, 6, 2, 3};
  TS_ASSERT_EQUALS_VEC(baseLoop, o1.m_loops[0]);
  std::vector<int> baseLoopType = {1};
  TS_ASSERT_EQUALS_VEC(baseLoopType, o1.m_loopTypes);
} // MePolyCleanerUnitTests::testCleanIn4a
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCleanIn4b()
{
  // x =     0                4                 8
  //
  // y=4  (4)0-----------0---------3(7)---------3
  //        |           |         |            |
  //        |           |         |            |
  //        |           |         |            |
  //        |           |         |            |
  //        |           |         |            |
  // y=2     |           |         |            |
  //        |           |         |            |
  //        |           |         |            |
  //        |           |         |            |
  //        |           |         |(6)         |
  // y=0  (5)1-----------1---------2------------2
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{3, 4, 0}, {3, 0, 0}, {8, 0, 0}, {8, 4, 0}};
  std::vector<size_t> v0 = {0, 1, 2, 3};
  o.m_loops.push_back(v0);
  o.m_loopTypes.assign(1, MePolyOffsetter::INSIDE_POLY);
  std::vector<MePolyOffsetterOutput> vO(1, o);
  o.m_pts = {{0, 4, 0}, {0, 0, 0}, {5, 0, 0}, {5, 4, 0}};
  vO.push_back(o);
  MePolyCleanerImpl p;
  p.IntersectCleanInPolys(vO, o1, 1e-9);
  TS_ASSERT_EQUALS(1, o1.m_loops.size());
  if (1 != o1.m_loops.size())
    return;
  std::vector<Pt3d> basePts = {{3, 4, 0}, {3, 0, 0}, {8, 0, 0}, {8, 4, 0},
                               {0, 4, 0}, {0, 0, 0}, {5, 0, 0}, {5, 4, 0}};
  TS_ASSERT_EQUALS_VEC(basePts, o1.m_pts);
  std::vector<size_t> baseLoop = {7, 0, 4, 5, 1, 6, 2, 3};
  TS_ASSERT_EQUALS_VEC(baseLoop, o1.m_loops[0]);
  std::vector<int> baseLoopType = {1};
  TS_ASSERT_EQUALS_VEC(baseLoopType, o1.m_loopTypes);
} // MePolyCleanerUnitTests::testCleanIn4b
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCleanIn4c()
{
  // x =     0                4                 8
  //
  // y=4  (4)0---------------------3(7)
  //        |                     |
  //        |                     |(9)
  //        0---------------------|------------3
  //        |                     |            |
  //        |                     |            |
  // y=2     |                     |            |
  //        |                     |(8)         |
  //        1---------------------|------------2
  //        |                     |
  //        |                     |
  // y=0  (5)1---------------------2(6)
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{0, 3, 0}, {0, 1, 0}, {8, 1, 0}, {8, 3, 0}};
  std::vector<size_t> v0 = {0, 1, 2, 3};
  o.m_loops.push_back(v0);
  o.m_loopTypes.assign(1, MePolyOffsetter::INSIDE_POLY);
  std::vector<MePolyOffsetterOutput> vO(1, o);
  o.m_pts = {{0, 4, 0}, {0, 0, 0}, {5, 0, 0}, {5, 4, 0}};
  vO.push_back(o);
  MePolyCleanerImpl p;
  p.IntersectCleanInPolys(vO, o1, 1e-9);
  TS_ASSERT_EQUALS(1, o1.m_loops.size());
  if (1 != o1.m_loops.size())
    return;
  std::vector<Pt3d> basePts = {{0, 3, 0}, {0, 1, 0}, {8, 1, 0}, {8, 3, 0}, {0, 4, 0},
                               {0, 0, 0}, {5, 0, 0}, {5, 4, 0}, {5, 1, 0}, {5, 3, 0}};
  TS_ASSERT_EQUALS_VEC(basePts, o1.m_pts);
  std::vector<size_t> baseLoop = {9, 7, 4, 0, 1, 5, 6, 8, 2, 3};
  TS_ASSERT_EQUALS_VEC(baseLoop, o1.m_loops[0]);
  std::vector<int> baseLoopType = {1};
  TS_ASSERT_EQUALS_VEC(baseLoopType, o1.m_loopTypes);
} // MePolyCleanerUnitTests::testCleanIn4c
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCleanIn4d()
{
  // x =     0                4                 8
  //
  // y=4                      1(5)
  //                        /|
  //                       / |
  //                     /   |
  //                    /    |
  //               (7)/      |
  // y=2     0----------------3
  //        |      /         |
  //        |    /           |
  //        |  /             |
  //        |/               |
  // y=0     1----------------2
  //        2(6)             0(4)
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{0, 2, 0}, {0, 0, 0}, {4, 0, 0}, {4, 2, 0}};
  std::vector<size_t> v0 = {0, 1, 2, 3};
  o.m_loops.push_back(v0);
  o.m_loopTypes.assign(1, MePolyOffsetter::INSIDE_POLY);
  std::vector<MePolyOffsetterOutput> vO(1, o);
  o.m_pts = {{4, 0, 0}, {4, 4, 0}, {0, 0, 0}};
  o.m_loops[0].pop_back();
  vO.push_back(o);
  MePolyCleanerImpl p;
  p.IntersectCleanInPolys(vO, o1, 1e-9);
  TS_ASSERT_EQUALS(1, o1.m_loops.size());
  if (1 != o1.m_loops.size())
    return;
  std::vector<Pt3d> basePts = {{0, 2, 0}, {0, 0, 0}, {4, 0, 0}, {4, 2, 0},
                               {4, 0, 0}, {4, 4, 0}, {0, 0, 0}, {2, 2, 0}};
  TS_ASSERT_EQUALS_VEC(basePts, o1.m_pts);
  std::vector<size_t> baseLoop = {1, 2, 3, 5, 7, 0};
  TS_ASSERT_EQUALS_VEC(baseLoop, o1.m_loops[0]);
  std::vector<int> baseLoopType = {1};
  TS_ASSERT_EQUALS_VEC(baseLoopType, o1.m_loopTypes);
} // MePolyCleanerUnitTests::testCleanIn4d
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCleanIn5()
{
  // x =     0                4                 8
  //
  // y=4                 0----------------------3
  //                    |                      |
  //                    |                      |
  //              0     |                      |
  //              | \   |                      |
  //              |   \ |                      |
  // y=2           |     |\                     |
  //              |     |  \                   |
  //              |     |     \                |
  //              |     1----------------------2
  //                 |               \
//y=0           1-----------------2
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{3, 4, 0}, {3, 1, 0}, {8, 1, 0}, {8, 4, 0}};
  std::vector<size_t> v0 = {0, 1, 2, 3};
  o.m_loops.push_back(v0);
  o.m_loopTypes.assign(1, MePolyOffsetter::INSIDE_POLY);
  std::vector<MePolyOffsetterOutput> vO(1, o);
  o.m_pts = {{2, 3, 0}, {2, 0, 0}, {5, 0, 0}};
  o.m_loops[0].pop_back();
  vO.push_back(o);
  MePolyCleanerImpl p;
  p.IntersectCleanInPolys(vO, o1, 1e-9);
  TS_ASSERT_EQUALS(1, o1.m_loops.size());
  if (1 != o1.m_loops.size())
    return;
  std::vector<Pt3d> basePts = {{3, 4, 0}, {3, 1, 0}, {8, 1, 0}, {8, 4, 0}, {2, 3, 0},
                               {2, 0, 0}, {5, 0, 0}, {3, 2, 0}, {4, 1, 0}};
  TS_ASSERT_EQUALS_VEC(basePts, o1.m_pts);
  std::vector<size_t> baseLoop = {8, 2, 3, 0, 7, 4, 5, 6};
  TS_ASSERT_EQUALS_VEC(baseLoop, o1.m_loops[0]);
  std::vector<int> baseLoopType = {1};
  TS_ASSERT_EQUALS_VEC(baseLoopType, o1.m_loopTypes);
} // MePolyCleanerUnitTests::testCleanIn5
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCleanIn6()
{
  // x =     0                4                 8
  //
  // y=8
  //      (4)                     (8)
  //       0----------------------4
  //       |                   /
  //       |            0----------------------3
  //       |            | /  (12)              |
  // y=4    |            |                      |
  //       |          / |(9)                   |
  //       |        /   |                      |
  //       |      3(7)  |                      |
  //       |        \   |                      |
  //       |          \ |                      |
  // y=2    |            |\(10)                 |
  //       |            |   \                  |
  //       |            |      \ (11)          |
  //       |            1----------------------2
  //          |                      \
//y=0    1------------------------2
  //      (5)                      (6)
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{3, 5, 0}, {3, 1, 0}, {8, 1, 0}, {8, 5, 0}};
  std::vector<size_t> v0 = {0, 1, 2, 3};
  o.m_loops.push_back(v0);
  o.m_loopTypes.assign(1, MePolyOffsetter::INSIDE_POLY);
  std::vector<MePolyOffsetterOutput> vO(1, o);
  o.m_pts = {{0, 6, 0}, {0, 0, 0}, {5, 0, 0}, {2, 3, 0}, {5, 6, 0}};
  o.m_loops[0] = {0, 1, 2, 3, 4};
  vO.push_back(o);
  MePolyCleanerImpl p;
  p.IntersectCleanInPolys(vO, o1, 1e-9);
  TS_ASSERT_EQUALS(2, o1.m_loops.size());
  if (2 != o1.m_loops.size())
    return;
  std::vector<Pt3d> basePts = {{3, 5, 0}, {3, 1, 0}, {8, 1, 0}, {8, 5, 0}, {0, 6, 0},
                               {0, 0, 0}, {5, 0, 0}, {2, 3, 0}, {5, 6, 0}, {3, 4, 0},
                               {3, 2, 0}, {4, 1, 0}, {4, 5, 0}};
  TS_ASSERT_EQUALS_VEC(basePts, o1.m_pts);
  std::vector<size_t> baseLoop = {9, 10, 7};
  TS_ASSERT_EQUALS_VEC(baseLoop, o1.m_loops[0]);
  baseLoop = {12, 8, 4, 5, 6, 11, 2, 3};
  TS_ASSERT_EQUALS_VEC(baseLoop, o1.m_loops[1]);
  std::vector<int> baseLoopType = {2, 1};
  TS_ASSERT_EQUALS_VEC(baseLoopType, o1.m_loopTypes);
} // MePolyCleanerUnitTests::testCleanIn6
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCleanIn7()
{
  // x =     0                4                 8
  //
  // y=8
  //
  //                      3(7)
  //                (15)/    \ (8)
  //              0----------------3
  //              |  /          \  |
  // y=4           |/              \|
  //            / |(14)         (9)| \
//          /   |                |    \
//     (4)0     |                |      2(6)
  //         \    |                |     /
  //           \  |                |   /
  // y=2          \|(13)        (10)|/
  //              |\              /|
  //              |  \ (12) (11)/  |
  //              1----------------2
  //                    \    /
  // y=0                    1(5)
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{1, 5, 0}, {1, 1, 0}, {5, 1, 0}, {5, 5, 0}};
  std::vector<size_t> v0 = {0, 1, 2, 3};
  o.m_loops.push_back(v0);
  o.m_loopTypes.assign(1, MePolyOffsetter::INSIDE_POLY);
  std::vector<MePolyOffsetterOutput> vO(1, o);
  o.m_pts = {{0, 3, 0}, {3, 0, 0}, {6, 3, 0}, {3, 6, 0}};
  vO.push_back(o);
  MePolyCleanerImpl p;
  p.IntersectCleanInPolys(vO, o1, 1e-9);
  TS_ASSERT_EQUALS(1, o1.m_loops.size());
  if (1 != o1.m_loops.size())
    return;
  std::vector<Pt3d> basePts = {{1, 5, 0}, {1, 1, 0}, {5, 1, 0}, {5, 5, 0}, {0, 3, 0}, {3, 0, 0},
                               {6, 3, 0}, {3, 6, 0}, {4, 5, 0}, {5, 4, 0}, {5, 2, 0}, {4, 1, 0},
                               {2, 1, 0}, {1, 2, 0}, {1, 4, 0}, {2, 5, 0}};
  TS_ASSERT_EQUALS_VEC(basePts, o1.m_pts);
  std::vector<size_t> baseLoop = {15, 0, 14, 4, 13, 1, 12, 5, 11, 2, 10, 6, 9, 3, 8, 7};
  TS_ASSERT_EQUALS_VEC(baseLoop, o1.m_loops[0]);
  std::vector<int> baseLoopType = {1};
  TS_ASSERT_EQUALS_VEC(baseLoopType, o1.m_loopTypes);
} // MePolyCleanerUnitTests::testCleanIn7
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCleanIn8()
{
  // x =    0                 4                 8
  //
  // y=4            0--------------3
  //               |              |
  //           (12)|  (9)   (5)   |(13)
  //   (10)2-------|---1     1----|------------0(4)
  //       |       |   |     |    |            |
  //       |       |   |     |    |            |
  // y=2    |       |   |     |    |            |
  //       |       |   |     |    |            |
  //       |   (15)|   |     |    |(14)        |
  //   (11)3-------|---0     2----|------------3(7)
  //               |  (8)   (6)   |
  // y=0            1--------------2
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{2, 4, 0}, {2, 0, 0}, {5, 0, 0}, {5, 4, 0}};
  std::vector<size_t> v0{0, 1, 2, 3};
  o.m_loops.push_back(v0);
  o.m_loopTypes.assign(1, MePolyOffsetter::INSIDE_POLY);
  std::vector<MePolyOffsetterOutput> vO(1, o);
  o.m_pts = {{8, 3, 0}, {4, 3, 0}, {4, 1, 0}, {8, 1, 0}};
  vO.push_back(o);
  o.m_pts = {{3, 1, 0}, {3, 3, 0}, {0, 3, 0}, {0, 1, 0}};
  vO.push_back(o);
  MePolyCleanerImpl p;
  p.IntersectCleanInPolys(vO, o1, 1e-9);
  TS_ASSERT_EQUALS(1, o1.m_loops.size());
  if (1 != o1.m_loops.size())
    return;
  std::vector<Pt3d> basePts = {{2, 4, 0}, {2, 0, 0}, {5, 0, 0}, {5, 4, 0}, {8, 3, 0}, {4, 3, 0},
                               {4, 1, 0}, {8, 1, 0}, {3, 1, 0}, {3, 3, 0}, {0, 3, 0}, {0, 1, 0},
                               {2, 3, 0}, {5, 3, 0}, {5, 1, 0}, {2, 1, 0}};
  TS_ASSERT_EQUALS_VEC(basePts, o1.m_pts);
  std::vector<size_t> baseLoop = {15, 1, 2, 14, 7, 4, 13, 3, 0, 12, 10, 11};
  TS_ASSERT_EQUALS_VEC(baseLoop, o1.m_loops[0]);
  std::vector<int> baseLoopType = {1};
  TS_ASSERT_EQUALS_VEC(baseLoopType, o1.m_loopTypes);
} // MePolyCleanerUnitTests::testCleanIn8

//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCleanIn9()
{
  // x =    0        2        4        6        8
  //
  // y=6
  //
  //
  //                  (5)1-----------------0(4)
  //                     |                 |
  //                     |(14)             |
  // y=4    2--------------------------1    |
  //       |             |            |    |
  //       |             |      (11)  |    |
  //       |   2(10)     |        0   |    |
  //       |   |\        |       /|   |    |
  //       |   |  \(16)  |(15) /  |   |    |
  // y=2    3--------------------------0    |
  //       (17)|     \   |  /     |        |
  //           |    (9)\ |/       |        |
  //           0---------1        |        |
  //         (8)       / |        |        |
  // y=0          (12)1---2--------2--------3(7)
  //                    (6)      (13)
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{6, 2, 0}, {6, 4, 0}, {0, 4, 0}, {0, 2, 0}};
  std::vector<size_t> v0 = {0, 1, 2, 3};
  o.m_loops.push_back(v0);
  o.m_loopTypes.assign(1, MePolyOffsetter::INSIDE_POLY);
  std::vector<MePolyOffsetterOutput> vO(1, o);
  o.m_pts = {{7, 5, 0}, {3, 5, 0}, {3, 0, 0}, {7, 0, 0}};
  vO.push_back(o);
  o.m_pts = {{1, 1, 0}, {3, 1, 0}, {1, 3, 0}};
  o.m_loops[0].pop_back();
  vO.push_back(o);
  o.m_pts = {{5, 3, 0}, {2, 0, 0}, {5, 0, 0}};
  vO.push_back(o);
  MePolyCleanerImpl p;
  p.IntersectCleanInPolys(vO, o1, 1e-9);

#if BOOST_VERSION > 106100
  TS_ASSERT_EQUALS(2, o1.m_loops.size());

  if (2 != o1.m_loops.size())
  {
    // iDumpOutput(o1);
    return;
  }
  std::vector<Pt3d> basePts = {{6, 2, 0}, {6, 4, 0}, {0, 4, 0}, {0, 2, 0}, {7, 5, 0}, {3, 5, 0},
                               {3, 0, 0}, {7, 0, 0}, {1, 1, 0}, {3, 1, 0}, {1, 3, 0}, {5, 3, 0},
                               {2, 0, 0}, {5, 0, 0}, {2, 2, 0}, {3, 2, 0}, {3, 4, 0}, {1, 2, 0}};
  TS_ASSERT_EQUALS_VEC(basePts, o1.m_pts);
  std::vector<size_t> loop0 = {9, 14, 15};
  TS_ASSERT_EQUALS_VEC(loop0, o1.m_loops[0])
  std::vector<size_t> loop1 = {17, 8, 9, 12, 6, 13, 7, 4, 5, 16, 2, 3};
  TS_ASSERT_EQUALS_VEC(loop1, o1.m_loops[1])
  std::vector<int> baseLoopTypes = {2, 1};
  TS_ASSERT_EQUALS_VEC(baseLoopTypes, o1.m_loopTypes);
#else
  TS_ASSERT_EQUALS(1, o1.m_loops.size());
  if (1 != o1.m_loops.size())
    return;
  std::vector<Pt3d> basePts = {{6, 2, 0}, {6, 4, 0}, {0, 4, 0}, {0, 2, 0}, {7, 5, 0}, {3, 5, 0},
                               {3, 0, 0}, {7, 0, 0}, {1, 1, 0}, {3, 1, 0}, {1, 3, 0}, {5, 3, 0},
                               {2, 0, 0}, {5, 0, 0}, {3, 4, 0}, {3, 2, 0}, {2, 2, 0}, {1, 2, 0}};
  TS_ASSERT_EQUALS_VEC(basePts, o1.m_pts);
  std::vector<size_t> baseLoop = {17, 8, 9, 16, 15, 9, 12, 6, 13, 7, 4, 5, 14, 2, 3};
  TS_ASSERT_EQUALS_VEC(baseLoop, o1.m_loops[0]);
  std::vector<int> baseLoopType = {1};
  TS_ASSERT_EQUALS_VEC(baseLoopType, o1.m_loopTypes);
#endif
} // MePolyCleanerUnitTests::testCleanIn9
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCleanIn10()
{
  // x =    0        2        4        6        8
  //
  // y=6
  //
  //
  //       1-------------------------------0
  //       |                               |
  //       |                               |
  // y=4    |                               |
  //       |                               |
  //       |                               |
  //       |   2(8)                        |
  //       |   |\                          |
  //       |   |  \(9)                     |
  // y=2    2-------------3                 |
  //       (10)|     \   |                 |
  //           |       \ |                 |
  //           0---------1(7)              |
  //         (6)         |                 |
  // y=0                  4-----------------5
  //
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{7, 7, 0}, {0, 7, 0}, {0, 2, 0}, {3, 2, 0}, {3, 0, 0}, {7, 0, 0}};
  std::vector<size_t> v0 = {0, 1, 2, 3, 4, 5};
  o.m_loops.push_back(v0);
  o.m_loopTypes.assign(1, MePolyOffsetter::INSIDE_POLY);
  std::vector<MePolyOffsetterOutput> vO(1, o);
  o.m_pts = {{1, 1, 0}, {3, 1, 0}, {1, 3, 0}};
  o.m_loops[0] = {0, 1, 2};
  vO.push_back(o);
  MePolyCleanerImpl p;
  p.IntersectCleanInPolys(vO, o1, 1e-9);

#if BOOST_VERSION > 106100
  TS_ASSERT_EQUALS(2, o1.m_loops.size());
  if (2 != o1.m_loops.size())
  {
    // iDumpOutput(o1);
    return;
  }
  std::vector<Pt3d> basePts = {{7, 7, 0}, {0, 7, 0}, {0, 2, 0}, {3, 2, 0}, {3, 0, 0}, {7, 0, 0},
                               {1, 1, 0}, {3, 1, 0}, {1, 3, 0}, {2, 2, 0}, {1, 2, 0}};
  TS_ASSERT_EQUALS_VEC(basePts, o1.m_pts);
  std::vector<size_t> loop0 = {7, 9, 3};
  TS_ASSERT_EQUALS_VEC(loop0, o1.m_loops[0]);
  std::vector<size_t> loop1 = {7, 4, 5, 0, 1, 2, 10, 6};
  TS_ASSERT_EQUALS_VEC(loop1, o1.m_loops[1]);
  std::vector<int> baseLoopType = {2, 1};
  TS_ASSERT_EQUALS_VEC(baseLoopType, o1.m_loopTypes);
#else
  TS_ASSERT_EQUALS(1, o1.m_loops.size());
  if (1 != o1.m_loops.size())
    return;
  std::vector<Pt3d> basePts = {{7, 7, 0}, {0, 7, 0}, {0, 2, 0}, {3, 2, 0}, {3, 0, 0}, {7, 0, 0},
                               {1, 1, 0}, {3, 1, 0}, {1, 3, 0}, {2, 2, 0}, {1, 2, 0}};
  TS_ASSERT_EQUALS_VEC(basePts, o1.m_pts);
  std::vector<size_t> baseLoop = {10, 6, 7, 9, 3, 4, 5, 0, 1, 2};
  TS_ASSERT_EQUALS_VEC(baseLoop, o1.m_loops[0]);
  std::vector<int> baseLoopType = {1};
  TS_ASSERT_EQUALS_VEC(baseLoopType, o1.m_loopTypes);
#endif
} // MePolyCleanerUnitTests::testCleanIn10
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCleanIn11()
{
  // x =     0                4                 8
  //
  // y=4     7----------------------------------6
  //        |                                  |
  //        |                    15------14    |
  //        |   1-----------2     |      |     |
  //        |   |           |    12------13    |
  //        |   |  19---18  |                  |
  // y=2     |   |   |   |   |     8------11    |
  //        |   |  16---17  |     |      |     |
  //        |   |           |     |      |     |
  //        |   0-----------3     9------10    |
  //        |                                  |
  // y=0     4----------------------------------5
  //
  // Poly 0,1,2,3 is a NEWOUT_POLY all the others are INSIDE_POLY
  //
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{1, 1, 0},   {1, 3, 0},     {4, 3, 0},     {4, 1, 0},     {0, 0, 0},
             {8, 0, 0},   {8, 4, 0},     {0, 4, 0},     {5, 1, 0},     {7, 1, 0},
             {7, 2, 0},   {5, 2, 0},     {5, 3, 0},     {7, 3, 0},     {7, 3.5, 0},
             {5, 3.5, 0}, {1.5, 1.5, 0}, {3.5, 1.5, 0}, {3.5, 2.5, 0}, {1.5, 2.5, 0}};
  std::vector<size_t> v0{0, 1, 2, 3};
  o.m_loops.push_back(v0);
  std::vector<size_t> v1{4, 5, 6, 7};
  o.m_loops.push_back(v1);
  std::vector<size_t> v2{8, 9, 10, 11};
  o.m_loops.push_back(v2);
  std::vector<size_t> v3{12, 13, 14, 15};
  o.m_loops.push_back(v3);
  std::vector<size_t> v4{16, 17, 18, 19};
  o.m_loops.push_back(v4);
  o.m_loopTypes.assign(5, MePolyOffsetter::INSIDE_POLY);
  o.m_loopTypes[0] = MePolyOffsetter::NEWOUT_POLY;
  MePolyCleanerImpl p;
  std::vector<MePolyOffsetterOutput> vO(1, o);
  p.IntersectCleanInPolys(vO, o1, 1e-9);
  TS_ASSERT_EQUALS(3, o1.m_loops.size());
  if (3 != o1.m_loops.size())
    return;
  TS_ASSERT_EQUALS_VEC(o.m_pts, o1.m_pts);
  TS_ASSERT_EQUALS_VEC(o.m_loops[0], o1.m_loops[0]);
  TS_ASSERT_EQUALS_VEC(o.m_loops[1], o1.m_loops[1]);
  TS_ASSERT_EQUALS_VEC(o.m_loops[4], o1.m_loops[2]);
  std::vector<int> baseLoopType = {2, 1, 1};
  TS_ASSERT_EQUALS_VEC(baseLoopType, o1.m_loopTypes);
} // MePolyCleanerUnitTests::testCleanIn11
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCleanInOut0()
{
  // x =     0                4                 8
  //
  // y=4     7----------------6
  //        |                |
  //        |   IN           |
  //        4----------------5
  //
  //
  // y=2     1----------------2
  //        |                |
  //        |                |
  //        |    OUT         |
  //        |                |
  // y=0     0----------------3
  //
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{0, 0, 0}, {0, 2, 0}, {4, 2, 0}, {4, 0, 0},
             {0, 3, 0}, {0, 4, 0}, {4, 4, 0}, {4, 3, 0}};
  std::vector<size_t> v0 = {0, 1, 2, 3};
  o.m_loops.push_back(v0);
  std::vector<size_t> v1 = {4, 5, 6, 7};
  o.m_loops.push_back(v1);
  o.m_loopTypes = {MePolyOffsetter::OUTSIDE_POLY, MePolyOffsetter::INSIDE_POLY};
  MePolyCleanerImpl p;
  p.IntersectCleanInOutPolys(o, o1, 1e-9);
  TS_ASSERT_EQUALS(2, o1.m_loops.size());
  if (2 != o1.m_loops.size())
    return;
  TS_ASSERT_EQUALS_VEC(o.m_pts, o1.m_pts);
  TS_ASSERT_EQUALS_VEC(o.m_loops[1], o1.m_loops[0]);
  TS_ASSERT_EQUALS_VEC(o.m_loops[0], o1.m_loops[1]);
  std::vector<int> baseLoops = {MePolyOffsetter::INSIDE_POLY, MePolyOffsetter::OUTSIDE_POLY};
  TS_ASSERT_EQUALS_VEC(baseLoops, o1.m_loopTypes);
} // MePolyCleanerUnitTests::testCleanInOut0
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCleanInOut0a()
{
  // x =     0                4                 8
  //
  // y=4     1----------------------------------2
  //        |                                  |
  //        |                                  |
  //        |      4---------------------------3
  //        |      |
  //        |      |
  // y=2     |      |         6-------9
  //        |      |         |       |
  //        |      |         |       |
  //        | OUT  |         |  IN   |
  //        |      |         |       |
  // y=0     0------5         7-------8
  //
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{0, 0, 0}, {0, 4, 0}, {8, 4, 0}, {8, 3, 0}, {2, 3, 0},
             {2, 0, 0}, {4, 2, 0}, {4, 0, 0}, {6, 0, 0}, {6, 2, 0}};
  std::vector<size_t> v0{0, 1, 2, 3, 4, 5};
  o.m_loops.push_back(v0);
  std::vector<size_t> v1{6, 7, 8, 9};
  o.m_loops.push_back(v1);
  o.m_loopTypes = {MePolyOffsetter::OUTSIDE_POLY, MePolyOffsetter::INSIDE_POLY};
  MePolyCleanerImpl p;
  p.IntersectCleanInOutPolys(o, o1, 1e-9);
  TS_ASSERT_EQUALS(2, o1.m_loops.size());
  if (2 != o1.m_loops.size())
    return;
  TS_ASSERT_EQUALS_VEC(o.m_pts, o1.m_pts);
  TS_ASSERT_EQUALS_VEC(o.m_loops[1], o1.m_loops[0]);
  TS_ASSERT_EQUALS_VEC(o.m_loops[0], o1.m_loops[1]);
  std::vector<int> baseLoops = {MePolyOffsetter::INSIDE_POLY, MePolyOffsetter::OUTSIDE_POLY};
  TS_ASSERT_EQUALS_VEC(baseLoops, o1.m_loopTypes);
} // MePolyCleanerUnitTests::testCleanInOut0a
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCleanInOut0b()
{
  // x =     0                4                 8
  //
  // y=4     1----------------------------------2
  //        |                                  |
  //        |                                  |
  //        |                                  |
  //        |                                  |
  //        |                                  |
  // y=2     |        7-------6                 |
  //        |        |       |                 |
  //        |        |  IN   |                 |
  //        |        4-------5                 |
  //        | OUT                              |
  // y=0     0----------------------------------3
  //
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{0, 0, 0}, {0, 4, 0}, {8, 4, 0}, {8, 0, 0},
             {2, 1, 0}, {4, 1, 0}, {4, 2, 0}, {2, 2, 0}};
  std::vector<size_t> v0{0, 1, 2, 3};
  o.m_loops.push_back(v0);
  std::vector<size_t> v1{4, 5, 6, 7};
  o.m_loops.push_back(v1);
  o.m_loopTypes = {MePolyOffsetter::OUTSIDE_POLY, MePolyOffsetter::INSIDE_POLY};
  MePolyCleanerImpl p;
  p.IntersectCleanInOutPolys(o, o1, 1e-9);
  TS_ASSERT_EQUALS(2, o1.m_loops.size());
  if (2 != o1.m_loops.size())
    return;
  TS_ASSERT_EQUALS_VEC(o.m_pts, o1.m_pts);
  TS_ASSERT_EQUALS_VEC(o.m_loops[1], o1.m_loops[0]);
  TS_ASSERT_EQUALS_VEC(o.m_loops[0], o1.m_loops[1]);
  std::vector<int> baseLoops = {MePolyOffsetter::INSIDE_POLY, MePolyOffsetter::OUTSIDE_POLY};
  TS_ASSERT_EQUALS_VEC(baseLoops, o1.m_loopTypes);
} // MePolyCleanerUnitTests::testCleanInOut0b
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCleanInOut1()
{
  // x =    0                 4                 8
  //
  // y=4            0--------------1
  //               |              |
  //           (12)|              |(13)
  //      10-------|---9     5----|------------4
  //       |       |   |     |    |            |
  //       |       |   |     |    |            |
  // y=2    |       |   |     |    |            |
  //       |       |   |     |    |            |
  //       |   (15)|   |     |    |(14)        |
  //      11-------|---8     6----|------------7
  //               |              |
  // y=0            3--------------2
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{2, 4, 0}, {5, 4, 0}, {5, 0, 0}, {2, 0, 0}, {8, 3, 0}, {4, 3, 0},
             {4, 1, 0}, {8, 1, 0}, {3, 1, 0}, {3, 3, 0}, {0, 3, 0}, {0, 1, 0}};
  std::vector<size_t> v0 = {0, 1, 2, 3};
  o.m_loops.push_back(v0);
  std::vector<size_t> v1 = {4, 5, 6, 7};
  o.m_loops.push_back(v1);
  std::vector<size_t> v2 = {8, 9, 10, 11};
  o.m_loops.push_back(v2);

  o.m_loopTypes.assign(3, MePolyOffsetter::INSIDE_POLY);
  o.m_loopTypes[0] = MePolyOffsetter::OUTSIDE_POLY;
  MePolyCleanerImpl p;
  p.IntersectCleanInOutPolys(o, o1, 1e-9);
  TS_ASSERT_EQUALS(1, o1.m_loops.size());
  if (1 != o1.m_loops.size())
    return;
  std::vector<Pt3d> basePts = {{2, 4, 0}, {5, 4, 0}, {5, 0, 0}, {2, 0, 0}, {8, 3, 0}, {4, 3, 0},
                               {4, 1, 0}, {8, 1, 0}, {3, 1, 0}, {3, 3, 0}, {0, 3, 0}, {0, 1, 0},
                               {2, 3, 0}, {5, 3, 0}, {5, 1, 0}, {2, 1, 0}};
  TS_ASSERT_EQUALS_VEC(basePts, o1.m_pts);
  std::vector<size_t> baseLoop{12, 0, 1, 13, 5, 6, 14, 2, 3, 15, 8, 9};
  TS_ASSERT_EQUALS_VEC(baseLoop, o1.m_loops[0]);
  std::vector<int> baseLoopType = {MePolyOffsetter::OUTSIDE_POLY};
  TS_ASSERT_EQUALS_VEC(baseLoopType, o1.m_loopTypes);
} // MePolyCleanerUnitTests::testCleanInOut1
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCleanInOut1a()
{
  // x =    0                 4                 8
  //
  // y=4            0--------------3
  //               |              |
  //           (15)|              |(13)
  //      10-------|---11    7----|------------4
  //       |       |   |     |    |            |
  //       |       |   |     |    |            |
  // y=2    |       |   |     |    |            |
  //       |       |   |     |    |            |
  //       |   (14)|   |     |    |(12)        |
  //       9-------|---8     6----|------------5
  //               |              |
  // y=0            1--------------2
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{2, 4, 0}, {2, 0, 0}, {5, 0, 0}, {5, 4, 0}, {8, 3, 0}, {8, 1, 0},
             {4, 1, 0}, {4, 3, 0}, {3, 1, 0}, {0, 1, 0}, {0, 3, 0}, {3, 3, 0}};
  std::vector<size_t> v0 = {0, 1, 2, 3};
  o.m_loops.push_back(v0);
  std::vector<size_t> v1 = {4, 5, 6, 7};
  o.m_loops.push_back(v1);
  std::vector<size_t> v2 = {8, 9, 10, 11};
  o.m_loops.push_back(v2);

  o.m_loopTypes.assign(3, MePolyOffsetter::OUTSIDE_POLY);
  o.m_loopTypes[0] = MePolyOffsetter::INSIDE_POLY;
  MePolyCleanerImpl p;
  p.IntersectCleanInOutPolys(o, o1, 1e-9);
  TS_ASSERT_EQUALS(2, o1.m_loops.size());
  if (2 != o1.m_loops.size())
    return;
  std::vector<Pt3d> basePts = {{2, 4, 0}, {2, 0, 0}, {5, 0, 0}, {5, 4, 0}, {8, 3, 0}, {8, 1, 0},
                               {4, 1, 0}, {4, 3, 0}, {3, 1, 0}, {0, 1, 0}, {0, 3, 0}, {3, 3, 0},
                               {5, 1, 0}, {5, 3, 0}, {2, 1, 0}, {2, 3, 0}};
  TS_ASSERT_EQUALS_VEC(basePts, o1.m_pts);
  std::vector<size_t> baseLoop{12, 13, 4, 5};
  TS_ASSERT_EQUALS_VEC(baseLoop, o1.m_loops[0]);
  baseLoop = {14, 9, 10, 15};
  TS_ASSERT_EQUALS_VEC(baseLoop, o1.m_loops[1]);
  std::vector<int> baseLoopType(2, MePolyOffsetter::OUTSIDE_POLY);
  TS_ASSERT_EQUALS_VEC(baseLoopType, o1.m_loopTypes);
} // MePolyCleanerUnitTests::testCleanInOut1a
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyCleanerUnitTests::testCleanInOut1b()
{
  // x =    0                 4                 8
  //
  // y=4            0--------------3
  //               |              |
  //            (9)|              |(11)
  //       7-------|--------------|------------4
  //       |       |              |            |
  //       |       |              |            |
  // y=2    |       |              |            |
  //       |       |              |            |
  //       |    (8)|              |(10)        |
  //       6-------|--------------|------------5
  //               |              |
  // y=0            1--------------2
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{2, 4, 0}, {2, 0, 0}, {5, 0, 0}, {5, 4, 0},
             {8, 3, 0}, {8, 1, 0}, {0, 1, 0}, {0, 3, 0}};
  std::vector<size_t> v0 = {0, 1, 2, 3};
  o.m_loops.push_back(v0);
  std::vector<size_t> v1 = {4, 5, 6, 7};
  o.m_loops.push_back(v1);

  o.m_loopTypes.assign(2, MePolyOffsetter::OUTSIDE_POLY);
  o.m_loopTypes[0] = MePolyOffsetter::INSIDE_POLY;
  MePolyCleanerImpl p;
  p.IntersectCleanInOutPolys(o, o1, 1e-9);
  TS_ASSERT_EQUALS(2, o1.m_loops.size());
  if (2 != o1.m_loops.size())
    return;
  std::vector<Pt3d> basePts = {{2, 4, 0}, {2, 0, 0}, {5, 0, 0}, {5, 4, 0}, {8, 3, 0}, {8, 1, 0},
                               {0, 1, 0}, {0, 3, 0}, {2, 1, 0}, {2, 3, 0}, {5, 1, 0}, {5, 3, 0}};
  TS_ASSERT_EQUALS_VEC(basePts, o1.m_pts);
  std::vector<size_t> baseLoop{8, 6, 7, 9};
  TS_ASSERT_EQUALS_VEC(baseLoop, o1.m_loops[0]);
  baseLoop = {10, 11, 4, 5};
  TS_ASSERT_EQUALS_VEC(baseLoop, o1.m_loops[1]);
  std::vector<int> baseLoopType(2, MePolyOffsetter::OUTSIDE_POLY);
  TS_ASSERT_EQUALS_VEC(baseLoopType, o1.m_loopTypes);
} // MePolyCleanerUnitTests::testCleanInOut1b

//} // namespace xms
#endif // CXX_TEST
