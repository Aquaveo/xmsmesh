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
#include <xmsmesh/meshing/detail/MeRelaxer.h>

// 3. Standard library headers
#include <cfloat>

// 4. External library headers

// 5. Shared code headers
#include <xmscore/points/pt.h>
#include <xmscore/stl/vector.h>
#include <xmscore/math/math.h>
#include <xmscore/misc/carray.h>
#include <xmscore/misc/Observer.h>
#include <xmscore/misc/XmError.h>
#include <xmscore/misc/xmstype.h> // for XM_PI
#include <xmscore/misc/XmConst.h>
#include <xmsinterp/geometry/geoms.h>
#include <xmsinterp/triangulate/TrTin.h>
#include <xmsinterp/triangulate/triangles.h>

// 6. Non-shared code headers

//----- Forward declarations ---------------------------------------------------

//----- External globals -------------------------------------------------------

//----- Namespace declaration --------------------------------------------------

namespace xms
{
//----- Constants / Enumerations -----------------------------------------------

//----- Classes / Structs ------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
class MeRelaxerImpl : public MeRelaxer
{
public:
  /// Point flags
  enum RelaxFlagEnum { RELAX_RELAX = 1, RELAX_MATSLIDE = 2, RELAX_BNDSLIDE = 8 };
  /// Types of relaxation, or the relaxation algorithms
  enum RelaxTypeEnum { RELAXTYPE_AREA, RELAXTYPE_DENSITY, RELAXTYPE_ANGLE, RELAXTYPE_CIRCLE };

  MeRelaxerImpl();
  virtual ~MeRelaxerImpl();

  // protected:
  //  MeRelaxerImpl(const MeRelaxerImpl&);

private:
  XM_DISALLOW_COPY_AND_ASSIGN(MeRelaxerImpl)

public:
  virtual void Relax(/*MeshPolyEnum a_meshPolyEnum,*/
                     const VecInt& a_fixedPoints,
                     BSHP<TrTin> a_tin) override;
  //------------------------------------------------------------------------------
  /// \brief Sets the observer to report progress on the process.
  /// \param a_: The observer.
  //------------------------------------------------------------------------------
  virtual void SetObserver(BSHP<Observer> a_) override { m_observer = a_; }

private:
  void ComputeCentroids();
  void RelaxMarkedPoints(RelaxTypeEnum a_relaxType, int a_iteration, int a_numiterations);
  bool AreaRelax(int a_point, Pt3d& a_newLocation);
  bool AngleRelax(int a_point, Pt3d& a_newLocation);
  bool RelaxOnBoundary(int a_point, Pt3d& a_newLocation, int a_relaxType);

private:
  BSHP<TrTin> m_tin;           ///< triangles connecting points (mesh elements)
  const VecInt* m_fixedPoints; ///< locations that may not be relaxed
  VecPt3d m_centroids;         ///< The triangle centroids
  double m_slideangle;         ///< Not sure how to document this one
  VecInt m_flags;              ///< Flags for points of type RelaxFlagEnum
  BSHP<Observer> m_observer;   ///< observer class to report progress

}; // class MeRelaxerImpl

////////////////////////////////////////////////////////////////////////////////
/// \class MeRelaxerImpl
/// \brief Relaxes mesh points. Moves them around to form a better mesh.
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief
//------------------------------------------------------------------------------
MeRelaxerImpl::MeRelaxerImpl()
: m_tin()
/*, m_meshPolyEnum*/
, m_fixedPoints(nullptr)
, m_centroids()
, m_slideangle(5.0)
, m_flags()
, m_observer()
{
} // MeRelaxerImpl::MeRelaxerImpl
//------------------------------------------------------------------------------
/// \brief
//------------------------------------------------------------------------------
MeRelaxerImpl::~MeRelaxerImpl()
{
} // MeRelaxerImpl::~MeRelaxerImpl
//------------------------------------------------------------------------------
/// \brief Moves interior mesh points to create a better mesh.
///
/// Iteratively relaxes the interior nodes.  Also swaps edges to maintain the
/// Delauney criterion after each set of relaxations. Also merges triangles.
/// Compare to myiRelaxSwapAndMerge.
/// \param a_tin: The tin with triangles.
/// \param a_fixedPoints: Points that shouldn't be moved (boundaries,
///                       breaklines). 0-based point indices.
//------------------------------------------------------------------------------
void MeRelaxerImpl::Relax(/*MeshPolyEnum a_meshPolyEnum,*/
                          const VecInt& a_fixedPoints,
                          BSHP<TrTin> a_tin)
{
  XM_ENSURE_TRUE_VOID_NO_ASSERT(a_tin);

  // Set up member variables
  /// \todo Save the MeshPolyEnum
  // m_meshPolyEnum = a_meshPolyEnum;
  m_tin = a_tin;
  m_fixedPoints = &a_fixedPoints;

  ComputeCentroids();

  // Set up iterations
  int numiterations;
  RelaxTypeEnum relaxtype;
  /// \todo Different options for scalar paving
  // if (m_meshPolyEnum == MESH_SPAVE) {
  //  numiterations = 5; // 2 AREA, then 3 DENSITY
  //  relaxtype = RELAXTYPE_DENSITY;
  //}
  // else {
  numiterations = 3; // 3 AREA
  relaxtype = RELAXTYPE_AREA;
  //}

  m_flags.assign(m_tin->Points().size(), 0);

  // Make sure we have adjaceny info
  if (m_tin->TrisAdjToPts().empty())
  {
    m_tin->BuildTrisAdjToPts();
  }

  // Mark the points to be relaxed
  for (size_t i = 0; i < m_flags.size(); ++i)
  {
    /// \todo don't relax on boundary
    // Do we really need to check if it's on the boundary? Shouldn't that be
    // included in a_fixedPoints?
    if (/*elBoundaryNode(point) ||*/ m_tin->TrisAdjToPts()[i].empty())
      m_flags[i] = 0;
    else
      m_flags[i] = RELAX_RELAX;
  }

  // Mark fixed points to not be relaxed (breaklines)
  for (size_t i = 0; i < a_fixedPoints.size(); ++i)
  {
    XM_ASSERT(a_fixedPoints[i] >= 0 && a_fixedPoints[i] < (int)m_flags.size());
    m_flags[a_fixedPoints[i]] = 0;
  }

  int iteration = 0;

  // if (m_meshPolyEnum == MESH_SPAVE) {
  /// \todo Do two more area relaxes for density meshing
  //}

  // Perform relaxation and swap
  for (int i = 0; i < 3; ++i)
  {
    // Relax marked nodes
    if (i < 2)
    {
      ++iteration;
      RelaxMarkedPoints(relaxtype, iteration, numiterations);
    }
    else
    { // i == 2
      ++iteration;
      RelaxMarkedPoints(relaxtype, iteration, numiterations);
    }

    m_tin->OptimizeTriangulation();
  }

} // MeRelaxerImpl::Relax
//------------------------------------------------------------------------------
/// \brief Computes the centroids of each triangle.
//------------------------------------------------------------------------------
void MeRelaxerImpl::ComputeCentroids()
{
  size_t numTri = m_tin->NumTriangles();
  m_centroids.assign(numTri, Pt3d());
  for (size_t t = 0; t < numTri; ++t)
  {
    m_centroids[t] = m_tin->TriangleCentroid((int)t);
  }
} // MeRelaxerImpl::ComputeCentroids
//------------------------------------------------------------------------------
/// \brief Relaxes the points marked by m_flags. Compare to rlRelaxMarkedNodes.
/// \param a_relaxType: RelaxTypeEnum.
/// \param a_iteration: Current iteration (for progress).
/// \param a_numiterations: Total iterations we will do (for progress).
//------------------------------------------------------------------------------
void MeRelaxerImpl::RelaxMarkedPoints(RelaxTypeEnum a_relaxType,
                                      int a_iteration,
                                      int a_numiterations)
{
  XM_ASSERT(a_iteration > 0 && a_numiterations > 0 && a_iteration <= a_numiterations);

  /// \todo  a_relaxType == RELAXTYPE_DENSITY or RELAXTYPE_CIRCLE

  // Progress
  double iterationPercent = a_iteration / (double)a_numiterations;

  // Relax the marked nodes
  VecPt3d& points = m_tin->Points(); // for convenience
  VecInt2d& trisAdjToPts = m_tin->TrisAdjToPts();
  size_t nPoints = points.size();
  bool ok;
  Pt3d origlocation, newlocation;
  for (size_t p = 0; p < nPoints; ++p)
  {
    if (m_flags[p])
    {
      ok = true;
      origlocation = points[p];

      // do general point relax
      if (m_flags[p] & RELAX_RELAX)
      {
        switch (a_relaxType)
        {
        case RELAXTYPE_AREA:
          ok = AreaRelax((int)p, newlocation);
          break;
        case RELAXTYPE_DENSITY:
          /// \todo DensityRelax
          XM_ASSERT(false);
          break;
        case RELAXTYPE_ANGLE:
          ok = AngleRelax((int)p, newlocation);
          break;
        case RELAXTYPE_CIRCLE:
          /// \todo CircleRelax
          XM_ASSERT(false);
          break;
        default:
          XM_ASSERT(false);
          break;
        }
      }

      // Allow nodes on material boundaries and on the mesh boundary to slide
      else if (m_flags[p] & RELAX_MATSLIDE || m_flags[p] & RELAX_BNDSLIDE)
      {
        ok = RelaxOnBoundary((int)p, newlocation, a_relaxType);
      }

      // Update the z value of relaxed point
      if (ok)
      {
        newlocation.z = points[p].z;
      }

      // Change the points location and check for bad triangles
      if (ok)
      {
        points[p] = newlocation;
        if (a_relaxType != RELAXTYPE_CIRCLE)
        {
          const VecInt& adjTris = trisAdjToPts[p];
          for (size_t i = 0; i < adjTris.size() && ok; ++i)
          {
            // Is triangle ill-formed? See that it has some area
            ok = (GT_EPS(m_tin->TriangleArea(adjTris[i]), 0.0, FLT_EPSILON));
          }
        }
      }

      if (!ok)
      {
        // Reset the point and elements
        points[p] = origlocation;
      }

      if (ok && (a_relaxType == RELAXTYPE_CIRCLE || a_relaxType == RELAXTYPE_DENSITY))
      {
        /// \todo recalculate weight size for point since it was moved
        XM_ASSERT(false);
      }
    } // if (m_flags[p])

    // Update progress
    if (m_observer)
    {
      double pcntdone = iterationPercent * p / (double)nPoints;
      m_observer->ProgressStatus(pcntdone);
    }
  } // for (int p = 0; p < nPoints; ++p)
} // MeRelaxerImpl::RelaxMarkedPoints
//------------------------------------------------------------------------------
/// \brief Relax a point using area of surrounding triangles. Trys to move
///        point to the center of the area. Compare to rliAreaRelax.
/// \param[in] a_point: The point index to be relaxed.
/// \param[out] a_newLocation: The location after point is relaxed (moved).
/// \return true if there were no errors.
//------------------------------------------------------------------------------
bool MeRelaxerImpl::AreaRelax(int a_point, Pt3d& a_newLocation)
{
  double sumx(0.0), sumy(0.0), area, sumarea(0.0);
  const VecInt& adjTris = m_tin->TrisAdjToPts()[a_point];
  for (size_t i = 0; i < adjTris.size(); ++i)
  {
    area = m_tin->TriangleArea(adjTris[i]);
    sumx += area * m_centroids[adjTris[i]].x;
    sumy += area * m_centroids[adjTris[i]].y;
    sumarea += area;
  }

  if (sumarea < XM_DBL_LOWEST)
  {
    return false;
  }
  a_newLocation.x = sumx / sumarea;
  a_newLocation.y = sumy / sumarea;
  return true;
} // MeRelaxerImpl::AreaRelax
//------------------------------------------------------------------------------
/// \brief Relax a point using angles. Moves point to equalize the angles
///        between the edges touching the point. Compare to rliAngleRelax.
/// \param[in] a_point: The point index to be relaxed.
/// \param[out] a_newLocation: The location after point is relaxed (moved).
/// \return true if there were no errors.
//------------------------------------------------------------------------------
bool MeRelaxerImpl::AngleRelax(int a_point, Pt3d& a_newLocation)
{
  int id;
  double targetangle, mag, a, c, d;
  double da, db, dx, dy, sx, sy, radius;
  Pt2d sum, pt, p1, p2, p3, centroid;
  int point1, point2;

  // get the target angle
  const VecInt& adjTris = m_tin->TrisAdjToPts()[a_point];
  VecPt3d& points = m_tin->Points();
  size_t num = adjTris.size();
  XM_ASSERT(num > 0);

  targetangle = 2.0 * XM_PI / (double)num;
  // loop through adjacent elems and find the
  // circle with the two adjacent points and
  // the new point
  sum.x = sum.y = 0.0;
  for (size_t i = 0; i < adjTris.size(); ++i)
  {
    int tri = adjTris[i];
    id = m_tin->LocalIndex(tri, a_point);
    point1 = m_tin->GlobalIndex(tri, trIncrementIndex(id));
    point2 = m_tin->GlobalIndex(tri, trDecrementIndex(id));
    // set up points
    pt.x = points[a_point].x;
    pt.y = points[a_point].y;
    p1.x = points[point1].x;
    p1.y = points[point1].y;
    p2.x = points[point2].x;
    p2.y = points[point2].y;
    // find length of line bisecting p1 and p2
    dx = p1.x - p2.x;
    dy = p1.y - p2.y;
    c = sqrt(sqr(dx) + sqr(dy));
    sx = p1.x + dx / 2.0;
    sy = p1.y + dy / 2.0;
    // compute length of triangle sides
    a = (c / 2.0) / sin(targetangle / 2.0);
    d = a * cos(targetangle / 2.0);
    radius = (a * a) / (2.0 * d);
    // find candidate point p3 and circle centroid
    dx /= c;
    dy /= c;
    p3.x = sx + d * dy;
    p3.y = sy - d * dx;
    da = sqr(p3.x - pt.x) + sqr(p3.y - pt.y);
    p3.x = sx - d * dy;
    p3.y = sy + d * dx;
    db = sqr(p3.x - pt.x) + sqr(p3.y - pt.y);
    if (da < db)
    {
      centroid.x = sx + (d - radius) * dy;
      centroid.y = sy - (d - radius) * dx;
    }
    else
    {
      centroid.x = sx - (d - radius) * dy;
      centroid.y = sy + (d - radius) * dx;
    }
    // find p3 along circle, radius away from centroid
    dx = pt.x - centroid.x;
    dy = pt.y - centroid.y;
    mag = sqrt(dx * dx + dy * dy);
    dx /= mag;
    dy /= mag;
    // find new point
    p3.x = centroid.x + radius * dx;
    p3.y = centroid.y + radius * dy;
    // sum new point's x and y
    sum.x += p3.x;
    sum.y += p3.y;
  }
  a_newLocation.x = (sum.x / (double)num);
  a_newLocation.y = (sum.y / (double)num);
  return true;
} // MeRelaxerImpl::AngleRelax
//------------------------------------------------------------------------------
/// \brief Relax along a material or outer boundary. Compare to
///        rliRelaxOnBoundary.
/// \param[in] a_point: The point index to be relaxed.
/// \param[out] a_newLocation: The location after point is relaxed (moved).
/// \param[in] a_relaxType: RelaxTypeEnum.
/// \return true if there were no errors.
//------------------------------------------------------------------------------
bool MeRelaxerImpl::RelaxOnBoundary(int a_point, Pt3d& a_newLocation, int a_relaxType)
{
  double angle, minangle, pcnt, dxp, dyp, dxn, dyn, dx1, dx2, dy1, dy2, totdist, dist1, dist2;

  VecPt3d& points = m_tin->Points();

  minangle = m_slideangle * (XM_PI / 180.0);
  int nextPoint(-1), prevPoint(-1);

  if (m_flags[a_point] == RELAX_MATSLIDE)
  {
    /// \todo make sure point only touches two materials
    XM_ASSERT(false);
  }
  else
  {
    prevPoint = m_tin->PreviousBoundaryPoint(a_point);
    nextPoint = m_tin->NextBoundaryPoint(a_point);
  }

  if (prevPoint != -1 && nextPoint != -1)
  {
    dxp = points[prevPoint].x - points[a_point].x;
    dyp = points[prevPoint].y - points[a_point].y;
    dxn = points[nextPoint].x - points[a_point].x;
    dyn = points[nextPoint].y - points[a_point].y;
    angle = gmAngleBetween2DVectors(dxp, dyp, dxn, dyn);
    angle -= XM_PI;
    if (fabs(angle) > minangle)
    {
      return false;
    }

    dx1 = points[nextPoint].x - points[a_point].x;
    dy1 = points[nextPoint].y - points[a_point].y;
    dx2 = points[prevPoint].x - points[a_point].x;
    dy2 = points[prevPoint].y - points[a_point].y;
    dist1 = sqrt(sqr(dx1) + sqr(dy1));
    dist2 = sqrt(sqr(dx2) + sqr(dy2));
    totdist = dist1 + dist2;

    // size func based relax
    switch (a_relaxType)
    {
    case RELAXTYPE_DENSITY:
    case RELAXTYPE_CIRCLE:
      /// \todo handle RELAXTYPE_DENSITY && RELAXTYPE_CIRCLE
      XM_ASSERT(false);
      break;
    case RELAXTYPE_AREA:
    case RELAXTYPE_ANGLE:
      if (dist1 > dist2)
      {
        pcnt = (0.5 * totdist - dist2) / dist1;
        a_newLocation.x = pcnt * dx1 + points[a_point].x;
        a_newLocation.y = pcnt * dy1 + points[a_point].y;
      }
      else if (dist2 > dist1)
      {
        pcnt = (0.5 * totdist - dist1) / dist2;
        a_newLocation.x = pcnt * dx2 + points[a_point].x;
        a_newLocation.y = pcnt * dy2 + points[a_point].y;
      }
      break;
    default:
      XM_ASSERT(0);
      return false;
      break;
    }
  }

  return true;
} // MeRelaxerImpl::RelaxOnBoundary

////////////////////////////////////////////////////////////////////////////////
/// \class MeRelaxer
/// \brief Relaxes mesh points.
/// \see MeRelaxerImpl
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief Creates a new instance of the class
/// \return MeRelaxer.
//------------------------------------------------------------------------------
BSHP<MeRelaxer> MeRelaxer::New()
{
  BSHP<MeRelaxer> polyMesher(new MeRelaxerImpl);
  // MeRelaxerImpl->SetBSHP(polyMesher); // If MeRelaxerImpl needs ptr to MeRelaxer
  return polyMesher;
} // MeRelaxer::New
//------------------------------------------------------------------------------
/// \brief Constructor
//------------------------------------------------------------------------------
MeRelaxer::MeRelaxer()
{
} // MeRelaxer::MeRelaxer
//------------------------------------------------------------------------------
/// \brief Destructor
//------------------------------------------------------------------------------
MeRelaxer::~MeRelaxer()
{
} // MeRelaxer::~MeRelaxer

} // namespace xms

#if CXX_TEST
////////////////////////////////////////////////////////////////////////////////
// UNIT TESTS
////////////////////////////////////////////////////////////////////////////////

#include <xmsmesh/meshing/detail/MeRelaxer.t.h>

#include <boost/assign.hpp>

#include <xmscore/testing/TestTools.h>
#include <xmsinterp/triangulate/TrTriangulatorPoints.h>

//----- Namespace declaration --------------------------------------------------

// namespace xms {
using namespace xms;

namespace
{
//------------------------------------------------------------------------------
/// \brief Convert an array to a vector of Pt3d.
//------------------------------------------------------------------------------
static VecPt3d iArrayToVecPt3d(double* a_array, int a_size)
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
/// \class MeRelaxerUnitTests
/// \brief Tests for MeRelaxer.
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief Tests the relaxation code as if we were creating a mesh.
/// \verbatim
//  40-   20-----21-----22-----23-----24
//    |   |\2   1|\2   1|\0   2|\1   0|
//    |   |0\ 12 |1\ 26 |2\ 25 |1\ 28 |
//    |   |  \   |  \   |  \   |  \   |
//    |   | 8 \  |   \  |   \  |   \  |
//    |   |    \0| 11 \0| 23 \1| 27 \2|
//    |   |1   2\|2   0\|0   1\|2   0\|
//  30-   15-----16-----17-----18-----19
//    |   |\2   1|\2   1|0   2/|\1   0|
//    |   |0\    |1\ 31 | 24 /2|2\ 29 |
//    |   |  \ 9 |  \   |   /  |  \   |
//    |   | 7 \  |   \  |  /   |   \  |
//    |   |    \0| 10 \0|1/ 22 | 30 \2|
//    |   |1   2\|2   0\|/0   1|0   1\|
//  20-   10-----11-----12-----13-----14
//    |   |0   2/|\1   0|\0   2|\0   2|
//    |   | 13 /2|2\ 14 |2\ 18 |2\ 21 |
//    |   |   /  |  \   |  \   |  \   |
//    |   |  / 6 | 5 \  |   \  |   \  |
//    |   |1/    |    \2| 15 \1| 19 \1|
//    |   |/0   1|0   1\|0   1\|0   1\|
//  10-   5------6------7------8------9
//    |   |\0   2|\0   2|\1   0|\1   0|
//    |   |2\    |1\    |0\ 16 |0\ 20 |
//    |   |  \ 1 |  \ 2 |  \   |  \   |
//    |   | 0 \  | 3 \  | 4 \  |   \  |
//    |   |    \1|    \1|    \2| 17 \2|
//    |   |0   1\|2   0\|1   2\|1   2\|
//  0 -   0------1------2------3------4
//
//        |------|------|------|------|
//        0     10     20     30     40
/// \endverbatim
//------------------------------------------------------------------------------
void MeRelaxerUnitTests::testRelaxWhileMeshing()
{
  // Create tin and get some convenience variables
  BSHP<TrTin> tin = TrTin::New();
  VecInt2d& trisAdjToPts = tin->TrisAdjToPts();

  // Set up to our example tin points
  tin->Points().reserve(5 * 5);
  for (size_t i = 0; i < 5; ++i)
  {
    for (size_t j = 0; j < 5; ++j)
    {
      tin->Points().push_back(Pt3d(j * 10.0, i * 10.0, 0.0));
    }
  }

  // Triangulate the points
  TrTriangulatorPoints client(tin->Points(), tin->Triangles(), &trisAdjToPts);
  client.Triangulate();

  // Set up fixed points (outer boundary)
  int outPolyA[] = {0, 1, 2, 3, 4, 9, 14, 19, 24, 23, 22, 21, 20, 15, 10, 5};
  VecInt outPoly(&outPolyA[0], &outPolyA[XM_COUNTOF(outPolyA)]);

  BSHP<MeRelaxer> relaxer = MeRelaxer::New();
  relaxer->Relax(outPoly, tin);
  const double kDelta = 1e-5;
  double pointsAfterA[] = {0.0,
                           0.0,
                           10.0,
                           0.0,
                           20.0,
                           0.0,
                           30.0,
                           0.0,
                           40.0,
                           0.0,
                           0.0,
                           10.0,
                           11.095340626589,
                           8.8783699489561,
                           20.021946561235,
                           10.198352688628,
                           29.943135400129,
                           9.9222857580709,
                           40.0,
                           10.0,
                           0.0,
                           20.0,
                           9.3320588918781,
                           19.085690838232,
                           20.677562368755,
                           20.864280263397,
                           31.148854873347,
                           18.861555902424,
                           40.0,
                           20.0,
                           0.0,
                           30.0,
                           9.9438849472478,
                           29.753838310984,
                           18.833441072562,
                           31.18892136462,
                           29.137598779841,
                           29.132561821603,
                           40.0,
                           30.0,
                           0.0,
                           40.0,
                           10.0,
                           40.0,
                           20.0,
                           40.0,
                           30.0,
                           40.0,
                           40.0,
                           40.0};
  VecPt3d pointsAfter = iArrayToVecPt3d(pointsAfterA, XM_COUNTOF(pointsAfterA));
  TS_ASSERT_DELTA_VECPT3D(pointsAfter, tin->Points(), kDelta);

} // MeRelaxerUnitTests::testRelaxWhileMeshing

//} // namespace xms
#endif // CXX_TEST
