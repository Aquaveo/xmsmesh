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
#include <xmsmesh/meshing/detail/MePolyPaverToMeshPts.h>

// 3. Standard library headers
#include <cfloat>
#include <list>

// 4. External library headers
#include <boost/unordered_set.hpp>

// 5. Shared code headers
#include <xmscore/points/pt.h>
#include <xmsinterp/geometry/geoms.h>
#include <xmsinterp/geometry/GmPtSearch.h>
#include <xmsmesh/meshing/detail/MeIntersectPolys.h>
#include <xmsmesh/meshing/detail/MePolyCleaner.h>
#include <xmsmesh/meshing/detail/MePolyOffsetter.h>
#include <xmsmesh/meshing/MePolyRedistributePts.h>
#include <xmscore/misc/Observer.h>
#include <xmscore/misc/XmError.h>
#include <xmscore/misc/XmConst.h>

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
namespace
{
class Poly
{
public:
  Poly()
  : m_iter(0)
  {
  }
  std::vector<Pt3d> m_outside;
  std::vector<std::vector<Pt3d>> m_inside;
  int m_iter;
};

class MePolyPaverToMeshPtsImpl : public MePolyPaverToMeshPts
{
public:
  MePolyPaverToMeshPtsImpl()
  : m_meshPts(nullptr)
  , m_xyTol(1e-9)
  , m_bias(1)
  , m_polyEnvelopeArea(0)
  , m_polyOffsetIter(1)
  , m_prog()
  {
  }

  bool PolyToMeshPts(const std::vector<Pt3d>& a_outPoly,
                     const std::vector<std::vector<Pt3d>>& a_inPolys,
                     double a_bias,
                     double a_xyTol,
                     std::vector<Pt3d>& a_meshPts);
  //------------------------------------------------------------------------------
  /// \brief
  //------------------------------------------------------------------------------
  void SetRedistributor(BSHP<MePolyRedistributePts> a_) override { m_externalRedist = a_; }
  void SetObserver(BSHP<Observer> a_) override { m_prog = a_; }
  void Setup();
  void TearDown();
  void ProcessStack();
  void AddPolygonToMeshPoints(const Poly& a_poly, bool a_first);
  void DoPave(const Poly& a_poly);
  void CleanPave();
  void RedistributePts();
  void ClassifyPolys();
  double AreaFromPolyStack();

  BSHP<VecPt3d> m_meshPts;
  std::list<Poly> m_polyStack;
  BSHP<MePolyOffsetter> m_offsetter;
  BSHP<MePolyCleaner> m_cleaner;
  BSHP<MePolyRedistributePts> m_redist;
  BSHP<MePolyRedistributePts> m_externalRedist;
  double m_xyTol;
  double m_bias;
  double m_polyEnvelopeArea;
  int m_polyOffsetIter;
  BSHP<Observer> m_prog;
  boost::unordered_set<std::pair<double, double>> m_ptHash;

  std::vector<MePolyOffsetterOutput> m_offsetOutputs;
};

} // unnamed namespace

//------------------------------------------------------------------------------
/// \brief Creates a new instance of this class
/// \return MePolyPaverToMeshPts.
//------------------------------------------------------------------------------
BSHP<MePolyPaverToMeshPts> MePolyPaverToMeshPts::New()
{
  BSHP<MePolyPaverToMeshPts> ret(new MePolyPaverToMeshPtsImpl);
  return ret;
} // MePolyPaverToMeshPts::New
//------------------------------------------------------------------------------
/// \brief
//------------------------------------------------------------------------------
MePolyPaverToMeshPts::MePolyPaverToMeshPts()
{
}
//------------------------------------------------------------------------------
/// \brief
//------------------------------------------------------------------------------
MePolyPaverToMeshPts::~MePolyPaverToMeshPts()
{
}

namespace
{
////////////////////////////////////////////////////////////////////////////////
/// \class MePolyPaverToMeshPtsImpl
/// \brief Generates seed points for a mesh from a polygon
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief Creates a vector of point locations for a mesh by paving a polygon.
/// The polygon may have interior holes.
/// \param a_outPoly Vector of points that define the outer loop of the polygon.
/// The loop closes on itself and the closing point is NOT repeated
/// \param a_insidePolys Vector of vectors of points that define the loops of
/// inside polygons.
/// \param a_bias Factor for transitioning between areas of high refinement to less
/// refinement.
/// \param a_xyTol Tolerance used for floating point geometry comparisons.
/// \param a_meshPts Vector of points generated by this method.
//------------------------------------------------------------------------------
bool MePolyPaverToMeshPtsImpl::PolyToMeshPts(const std::vector<Pt3d>& a_outPoly,
                                             const std::vector<std::vector<Pt3d>>& a_inPolys,
                                             double a_bias,
                                             double a_xyTol,
                                             std::vector<Pt3d>& a_meshPts)
{
  m_bias = a_bias;
  m_xyTol = a_xyTol;
  // Error checks
  if (a_outPoly.empty())
    return false;
  for (size_t i = 0; i < a_inPolys.size(); ++i)
  {
    if (a_inPolys[i].empty())
      return false;
  }

  m_polyStack.push_back(Poly());
  m_polyStack.back().m_outside = a_outPoly;
  m_polyStack.back().m_inside = a_inPolys;
  m_polyStack.back().m_iter = 1;

  Setup();
  ProcessStack();
  // fill the output variable
  a_meshPts.swap(*m_meshPts);
  TearDown();
  return true;
} // MePolyPaverToMeshPtsImpl::PolyToMeshPts
//------------------------------------------------------------------------------
/// \brief
//------------------------------------------------------------------------------
void MePolyPaverToMeshPtsImpl::Setup()
{
  // Calculate tolerance for point comparison
  m_meshPts = BSHP<VecPt3d>(new VecPt3d());
  m_offsetter = MePolyOffsetter::New();
  m_cleaner = MePolyCleaner::New();
  m_polyEnvelopeArea = AreaFromPolyStack();
  if (!m_externalRedist)
  {
    m_redist = MePolyRedistributePts::New();
    m_redist->SetSizeFuncFromPoly(m_polyStack.front().m_outside, m_polyStack.front().m_inside,
                                  m_bias);
  }
  else
    m_redist = m_externalRedist;
} // MePolyPaverToMeshPtsImpl::Setup
//------------------------------------------------------------------------------
/// \brief
//------------------------------------------------------------------------------
void MePolyPaverToMeshPtsImpl::TearDown()
{
  m_meshPts.reset();
  m_offsetter.reset();
  m_cleaner.reset();
  m_redist.reset();
  m_ptHash.clear();
} // MePolyPaverToMeshPtsImpl::TearDown
//------------------------------------------------------------------------------
/// \brief Processes a stack of polygons. Starts with the polygon that was
/// passed into the PolyToMeshPts method. On each paving iteration new polygons
/// are created and the points from the previous polygons are moved to the
/// vector of mesh node locations.
//------------------------------------------------------------------------------
void MePolyPaverToMeshPtsImpl::ProcessStack()
{
  if (m_prog)
    m_prog->BeginOperationString("Paving Polygon");

  double area;
  std::list<Poly>::iterator it(m_polyStack.begin());
  bool first(true);
  while (it != m_polyStack.end())
  {
    Poly& p(*it);

    m_polyOffsetIter = p.m_iter;
    // Pave the polygon
    DoPave(p);
    // clean the results from the pave
    CleanPave();
    // redistribute points on the polygons
    RedistributePts();
    // clean again after redistributing the points
    CleanPave();
    // classify the newly created polys into polygons defined by the "Poly"
    // class and put them on the stack
    ClassifyPolys();
    // add the points from this polygon to the mesh points
    AddPolygonToMeshPoints(p, first);
    first = false;

    // remove this poly from the stack and process the next one
    m_polyStack.erase(it);
    it = m_polyStack.begin();

    // do progress
    area = AreaFromPolyStack();
    if (m_prog)
      m_prog->ProgressStatus(area / m_polyEnvelopeArea);
  }
  if (m_prog)
    m_prog->EndOperation();
} // MePolyPaverToMeshPtsImpl::ProcessStack
//------------------------------------------------------------------------------
/// \brief Takes the points on a_poly and moves them to the output of mesh
/// node locations
//------------------------------------------------------------------------------
void MePolyPaverToMeshPtsImpl::AddPolygonToMeshPoints(const Poly& a_poly, bool // a_first
)
{
  auto itEnd = m_ptHash.end();
  std::pair<double, double> pt;
  for (size_t i = 0; i < a_poly.m_outside.size(); ++i)
  {
    pt.first = a_poly.m_outside[i].x;
    pt.second = a_poly.m_outside[i].y;
    if (m_ptHash.find(pt) == itEnd)
    {
      m_ptHash.insert(pt);
      m_meshPts->push_back(a_poly.m_outside[i]);
    }
  }
  for (size_t i = 0; i < a_poly.m_inside.size(); ++i)
  {
    for (size_t j = 0; j < a_poly.m_inside[i].size(); ++j)
    {
      pt.first = a_poly.m_inside[i][j].x;
      pt.second = a_poly.m_inside[i][j].y;
      if (m_ptHash.find(pt) == itEnd)
      {
        m_ptHash.insert(pt);
        m_meshPts->push_back(a_poly.m_inside[i][j]);
      }
    }
  }

  // if (a_first)
  //{
  //  BSHP<GmPtSearch> ptSearch = GmPtSearch::New(true);
  //  ptSearch->VectorThatGrowsToSearch(m_meshPts);
  //  int idx;
  //  for (size_t i=0; i<a_poly.m_outside.size(); ++i)
  //  {
  //    if (!ptSearch->AddPtToVectorIfUnique(a_poly.m_outside[i], m_xyTol, idx))
  //      idx = -2;
  //  }
  //  for (size_t i=0; i<a_poly.m_inside.size(); ++i)
  //  {
  //    for (size_t j=0; j<a_poly.m_inside[i].size(); ++j)
  //    {
  //      if (!ptSearch->AddPtToVectorIfUnique(a_poly.m_inside[i][j], m_xyTol, idx))
  //        idx = -2;
  //    }
  //  }
  //}
  // else
  //{
  //  m_meshPts->insert(m_meshPts->end(), a_poly.m_outside.begin(),
  //    a_poly.m_outside.end());
  //  for (size_t i=0; i<a_poly.m_inside.size(); ++i)
  //  {
  //    m_meshPts->insert(m_meshPts->end(), a_poly.m_inside[i].begin(),
  //      a_poly.m_inside[i].end());
  //  }
  //}

} // MePolyPaverToMeshPtsImpl::AddPolygonToMeshPoints
//------------------------------------------------------------------------------
/// \brief Paves inward from OUTSIDE_POLY and paves outward from INSIDE_POLY
//------------------------------------------------------------------------------
void MePolyPaverToMeshPtsImpl::DoPave(const Poly& a_poly)
{
  m_offsetOutputs.clear();

  MePolyOffsetterOutput offsetOut;
  MePolyOffsetter::polytype pIn(MePolyOffsetter::INSIDE_POLY), pOut(MePolyOffsetter::OUTSIDE_POLY);
  std::vector<Pt3d> pts;
  // pave in from the outer poly
  m_offsetter->Offset(a_poly.m_outside, pOut, offsetOut, m_xyTol);
  m_offsetOutputs.push_back(offsetOut);

  // pave out from the inner polys
  for (size_t i = 0; i < a_poly.m_inside.size(); ++i)
  {
    m_offsetter->Offset(a_poly.m_inside[i], pIn, offsetOut, m_xyTol);
    m_offsetOutputs.push_back(offsetOut);
  }
} // MePolyPaverToMeshPtsImpl::DoPave
//------------------------------------------------------------------------------
/// \brief Cleans up the results from the paving operation. There may be newly
/// created polygons that intersect.
//------------------------------------------------------------------------------
void MePolyPaverToMeshPtsImpl::CleanPave()
{
  // intersect the new inner polys with the other inner polys
  MePolyOffsetterOutput out, out2;
  m_cleaner->IntersectCleanInPolys(m_offsetOutputs, out, m_xyTol);

  // intersect the new inner polys with the outer polys
  m_cleaner->IntersectCleanInOutPolys(out, out2, m_xyTol);
  m_offsetOutputs.clear();
  m_offsetOutputs.push_back(out2);
} // MePolyPaverToMeshPtsImpl::CleanPave
//------------------------------------------------------------------------------
/// \brief Redistributes the points on polygons. A sizing function maybe used
/// that comes from a set of points.
//------------------------------------------------------------------------------
void MePolyPaverToMeshPtsImpl::RedistributePts()
{
  XM_ASSERT(m_redist);
  if (!m_redist)
    return;
  MePolyOffsetterOutput o;
  m_redist->Redistribute(m_offsetOutputs[0], o, m_polyOffsetIter);
  m_offsetOutputs[0] = o;
} // MePolyPaverToMeshPtsImpl::RedistributePts
//------------------------------------------------------------------------------
/// \brief Create new polygons to put onto the processing stack. These are the
/// polygons that are left after paving, cleaning, redistributing points, and
/// cleaning again.
/// \param [in] a_iter Tells what iteration this is from the boundary of the
/// polygon.
//------------------------------------------------------------------------------
void MePolyPaverToMeshPtsImpl::ClassifyPolys()
{
  if (m_offsetOutputs[0].m_pts.empty())
    return;

  std::vector<std::vector<size_t>> polys;
  MeIntersectPolys ip;
  ip.ClassifyPolys(m_offsetOutputs[0], polys);

  // put these on the stack
  std::vector<Pt3d>& pts(m_offsetOutputs[0].m_pts);
  Poly p;
  p.m_iter = m_polyOffsetIter + 1;
  for (size_t i = 0; i < polys.size(); ++i)
  {
    p.m_outside.resize(0);
    p.m_inside.resize(0);
    // get the outside loop
    std::vector<size_t>& oLoop(m_offsetOutputs[0].m_loops[polys[i][0]]);
    for (size_t j = 0; j < oLoop.size(); ++j)
    {
      p.m_outside.push_back(pts[oLoop[j]]);
    }
    // do inside loops
    for (size_t j = 1; j < polys[i].size(); ++j)
    {
      std::vector<size_t>& iLoop(m_offsetOutputs[0].m_loops[polys[i][j]]);
      p.m_inside.push_back(std::vector<Pt3d>());
      for (size_t k = 0; k < iLoop.size(); ++k)
      {
        p.m_inside.back().push_back(pts[iLoop[k]]);
      }
    }

    m_polyStack.push_back(p);
  }

} // MePolyPaverToMeshPtsImpl::ClassifyPolys
//------------------------------------------------------------------------------
/// \brief Computes the envelope area from a vector of points
//------------------------------------------------------------------------------
static double iEnvelopeArea(const std::vector<Pt3d>& a_pts)
{
  double area;
  Pt3d pMin(XM_DBL_HIGHEST),
    pMax(XM_DBL_LOWEST); // changed to XM_DBL_LOWEST because it is an arbitrarily low value
  for (size_t i = 0; i < a_pts.size(); ++i)
  {
    gmAddToExtents(a_pts[i], pMin, pMax);
  }
  area = (pMax.x - pMin.x) * (pMax.y - pMin.y);
  return area;
} // iEnvelopeArea
//------------------------------------------------------------------------------
/// \brief Computes the area of the polygon envelopes that are on the stack
//------------------------------------------------------------------------------
double MePolyPaverToMeshPtsImpl::AreaFromPolyStack()
{
  double area(0);
  // out polys are positive and in polys are negative
  std::list<Poly>::iterator it(m_polyStack.begin()), itEnd(m_polyStack.end());
  for (; it != itEnd; ++it)
  {
    Poly& p(*it);
    area += iEnvelopeArea(p.m_outside);
  }
  return area;
} // MePolyPaverToMeshPtsImpl::AreaFromPolyStack

} // unnamed namespace

} // namespace xms

#ifdef CXX_TEST
////////////////////////////////////////////////////////////////////////////////

#include <xmsmesh/meshing/detail/MePolyPaverToMeshPts.t.h>

#include <xmscore/testing/TestTools.h>

// namespace xms {
using namespace xms;

////////////////////////////////////////////////////////////////////////////////
/// \class MePolyPaverToMeshPtsUnitTests
/// \brief tester for the MePolyPaverToMeshPts class
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief tests creating a class
//------------------------------------------------------------------------------
void MePolyPaverToMeshPtsUnitTests::testCreateClass()
{
  BSHP<MePolyPaverToMeshPts> p = MePolyPaverToMeshPts::New();
  TS_ASSERT(p);
} // MeIntersectPolysTest::testCreateClass
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyPaverToMeshPtsUnitTests::testCase1()
{
  // x =     0                2                 4
  //
  // y=4     3--------4-------5-------6---------7
  //        |                                  |
  //        |                                  |
  //        2                                  8
  //        |                                  |
  //        |                                  |
  // y=2     1                                  9
  //        |                                  |
  //        |                                  |
  //        0                                  10
  //        |                                  |
  //        |                                  |
  // y=0    15--------14------13------12--------11
  //
  std::vector<Pt3d> outPoly = {{0, 1, 0}, {0, 2, 0}, {0, 3, 0}, {0, 4, 0}, {1, 4, 0}, {2, 4, 0},
                               {3, 4, 0}, {4, 4, 0}, {4, 3, 0}, {4, 2, 0}, {4, 1, 0}, {4, 0, 0},
                               {3, 0, 0}, {2, 0, 0}, {1, 0, 0}, {0, 0, 0}};
  std::vector<std::vector<Pt3d>> inPoly;
  MePolyPaverToMeshPtsImpl paver;
  double bias(1), tol(1e-9);
  std::vector<Pt3d> outPts;
  paver.PolyToMeshPts(outPoly, inPoly, bias, tol, outPts);
  std::vector<Pt3d> basePts = {
    {0.0, 1.0, 0.0},       {0.0, 2.0, 0.0},       {0.0, 3.0, 0.0},       {0.0, 4.0, 0.0},
    {1.0, 4.0, 0.0},       {2.0, 4.0, 0.0},       {3.0, 4.0, 0.0},       {4.0, 4.0, 0.0},
    {4.0, 3.0, 0.0},       {4.0, 2.0, 0.0},       {4.0, 1.0, 0.0},       {4.0, 0.0, 0.0},
    {3.0, 0.0, 0.0},       {2.0, 0.0, 0.0},       {1.0, 0.0, 0.0},       {0.0, 0.0, 0.0},
    {0.8660, 0.8660, 0.0}, {0.8660, 1.8740, 0.0}, {0.8660, 2.8820, 0.0}, {1.6220, 3.1340, 0.0},
    {2.6300, 3.1340, 0.0}, {3.1340, 2.6300, 0.0}, {3.1340, 1.6220, 0.0}, {2.8820, 0.8660, 0.0},
    {1.8740, 0.8660, 0.0}, {2.3148, 1.7390, 0.0}, {1.7390, 2.0535, 0.0}, {2.3323, 2.3802, 0.0}};
  TS_ASSERT_DELTA_VECPT3D(basePts, outPts, 1e-4);
} // MePolyPaverToMeshPtsUnitTests::testCase1
//------------------------------------------------------------------------------
/// \brief tests the figure below
//------------------------------------------------------------------------------
void MePolyPaverToMeshPtsUnitTests::testCase2()
{
  // x =     0                2                 4
  //
  // y=4     3--------4-------5-------6---------7
  //        |                                  |
  //        |                                  |
  //        2                                  8
  //        |                18--17            |
  //        |                 |  |             |
  // y=2     1                19--16            9
  //        |                                  |
  //        |                                  |
  //        0                                  10
  //        |                                  |
  //        |                                  |
  // y=0    15--------14------13------12--------11
  //
  std::vector<Pt3d> outPoly = {{0, 1, 0}, {0, 2, 0}, {0, 3, 0}, {0, 4, 0}, {1, 4, 0}, {2, 4, 0},
                               {3, 4, 0}, {4, 4, 0}, {4, 3, 0}, {4, 2, 0}, {4, 1, 0}, {4, 0, 0},
                               {3, 0, 0}, {2, 0, 0}, {1, 0, 0}, {0, 0, 0}};
  std::vector<std::vector<Pt3d>> inPoly;
  std::vector<Pt3d> tmp = {{2.25, 2, 0}, {2.25, 2.25, 0}, {2, 2.25, 0}, {2, 2, 0}};
  inPoly.push_back(tmp);
  MePolyPaverToMeshPtsImpl paver;
  double bias(1), tol(1e-9);
  std::vector<Pt3d> outPts;
  paver.PolyToMeshPts(outPoly, inPoly, bias, tol, outPts);
  std::vector<Pt3d> basePts = {
    {0.0, 1.0, 0.0},       {0.0, 2.0, 0.0},       {0.0, 3.0, 0.0},       {0.0, 4.0, 0.0},
    {1.0, 4.0, 0.0},       {2.0, 4.0, 0.0},       {3.0, 4.0, 0.0},       {4.0, 4.0, 0.0},
    {4.0, 3.0, 0.0},       {4.0, 2.0, 0.0},       {4.0, 1.0, 0.0},       {4.0, 0.0, 0.0},
    {3.0, 0.0, 0.0},       {2.0, 0.0, 0.0},       {1.0, 0.0, 0.0},       {0.0, 0.0, 0.0},
    {2.25, 2.0, 0.0},      {2.25, 2.25, 0.0},     {2.0, 2.25, 0.0},      {2.0, 2.0, 0.0},
    {0.8660, 0.8660, 0.0}, {0.8660, 1.7951, 0.0}, {0.8660, 2.7070, 0.0}, {1.3659, 3.1340, 0.0},
    {2.2505, 3.1340, 0.0}, {3.1340, 3.1224, 0.0}, {3.1340, 2.2316, 0.0}, {3.1340, 1.3448, 0.0},
    {2.6848, 0.8660, 0.0}, {1.7732, 0.8660, 0.0}, {1.7835, 2.1250, 0.0}, {2.0069, 1.7912, 0.0},
    {2.4260, 1.8985, 0.0}, {2.4335, 2.3440, 0.0}, {2.0093, 2.4589, 0.0}};
  TS_ASSERT_DELTA_VECPT3D(basePts, outPts, 1e-4);
} // MePolyPaverToMeshPtsUnitTests::testCase1

//} // namespace xms
#endif
