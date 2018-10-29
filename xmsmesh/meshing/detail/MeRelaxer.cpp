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
#include <xmscore/math/math.h>
#include <xmscore/misc/StringUtil.h>
#include <xmscore/misc/XmConst.h>
#include <xmscore/misc/XmError.h>
#include <xmscore/misc/carray.h>
#include <xmscore/misc/xmstype.h> // for XM_PI
#include <xmscore/points/pt.h>
#include <xmscore/stl/vector.h>
#include <xmsinterp/geometry/geoms.h>
#include <xmsinterp/interpolate/InterpBase.h>
#include <xmsinterp/triangulate/TrTin.h>
#include <xmsinterp/triangulate/triangles.h>
#include <xmsmesh/meshing/MePolyRedistributePts.h>

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
  enum RelaxFlagEnum { RELAX_RELAX = 1 };
  /// Types of relaxation, or the relaxation algorithms
  enum RelaxTypeEnum { RELAXTYPE_AREA, RELAXTYPE_ANGLE, RELAXTYPE_SPRING };

  MeRelaxerImpl();
  virtual ~MeRelaxerImpl();

  // protected:
  //  MeRelaxerImpl(const MeRelaxerImpl&);

  XM_DISALLOW_COPY_AND_ASSIGN(MeRelaxerImpl)

  virtual void Relax(/*MeshPolyEnum a_meshPolyEnum,*/
                     const VecInt& a_fixedPoints,
                     BSHP<TrTin> a_tin) override;
  virtual bool SetRelaxationMethod(const std::string& a_relaxMethod) override;
  //------------------------------------------------------------------------------
  /// \brief Sets size function used by the spring relaxation method
  /// \param a_sizer: The size function class
  //------------------------------------------------------------------------------
  virtual void SetPointSizer(BSHP<MePolyRedistributePts> a_sizer) override { m_sizer = a_sizer; }

  void ComputeCentroids();
  void RelaxMarkedPoints(RelaxTypeEnum a_relaxType, int a_iteration, int a_numiterations);
  void AreaRelax(int a_point, Pt3d& a_newLocation);
  void AngleRelax(int a_point, Pt3d& a_newLocation);
  void SpringRelaxSinglePoint(int a_point, Pt3d& a_newLocation);
  void SetupNeighbors();
  void SetupPointSizes();
  bool NewLocationIsValid(size_t a_idx, Pt3d& a_newLocation);
  bool AllTrianglesHavePositiveArea(BSHP<TrTin> a_tin);

  BSHP<TrTin> m_tin;           ///< triangles connecting points (mesh elements)
  const VecInt* m_fixedPoints; ///< locations that may not be relaxed
  VecPt3d m_centroids;         ///< The triangle centroids
  double m_slideangle;         ///< Used on boundary relaxation. Not sure how to document this one
  VecInt m_flags;              ///< Flags for points of type RelaxFlagEnum
  RelaxTypeEnum m_relaxType;   ///< the type of relaxation to perform. See RelaxTypeEnum
  BSHP<MePolyRedistributePts> m_sizer; ///< size function used by the spring relax method
  VecDbl m_pointSizes;                 ///< sizer size at each mesh point
  VecInt2d m_pointNeighbors;           ///< neighbor points for spring relaxation
  VecInt m_pointsToDelete;             ///< indexes of points that must be removed
};                                     // class MeRelaxerImpl

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
, m_relaxType(RELAXTYPE_AREA)
, m_sizer()
, m_pointSizes()
, m_pointsToDelete()
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
  XM_ENSURE_TRUE_VOID(!a_tin->TrisAdjToPts().empty());
  XM_ENSURE_TRUE_VOID(AllTrianglesHavePositiveArea(a_tin));
  // Set up member variables
  /// \todo Save the MeshPolyEnum
  // m_meshPolyEnum = a_meshPolyEnum;
  m_tin = a_tin;
  m_fixedPoints = &a_fixedPoints;

  // Mark fixed points to not be relaxed (breaklines)
  m_flags.assign(m_tin->Points().size(), RELAX_RELAX);
  for (size_t i = 0; i < a_fixedPoints.size(); ++i)
  {
    XM_ASSERT(a_fixedPoints[i] >= 0 && a_fixedPoints[i] < (int)m_flags.size());
    m_flags[a_fixedPoints[i]] = 0;
  }

  ComputeCentroids();

  // Set up iterations
  VecPt3d& pts(m_tin->Points());
  int numiterations(3); // default to 3 for area relaxation
  RelaxTypeEnum relaxtype(m_relaxType);
  if (relaxtype == RELAXTYPE_SPRING)
  {
    if (!m_sizer)
    {
      std::string msg =
        "No size function specified with spring relaxation "
        "method. Relaxation method has been set to AREA "
        "relaxation.";
      XM_LOG(xmlog::warning, msg);
      relaxtype = RELAXTYPE_AREA;
    }
    else
    {
      if (m_pointSizes.empty())
      {
        SetupPointSizes();
      }
      SetupNeighbors();
      XM_ASSERT(pts.size() == m_pointNeighbors.size());
      XM_ASSERT(pts.size() == m_pointSizes.size());
    }
  }
  XM_ASSERT(pts.size() == m_flags.size());

  int iteration = 0;

  // Perform relaxation and swap
  for (int i = 0; i < numiterations; ++i)
  {
    ++iteration;
    RelaxMarkedPoints(relaxtype, iteration, numiterations);
    // delete points and triangles
    if (!m_pointsToDelete.empty())
    {
      m_tin->Triangles().clear();
      m_tin->TrisAdjToPts().clear();
      VecPt3d& pts(m_tin->Points());
      VecDbl& szs(m_pointSizes);
      for (auto it = m_pointsToDelete.rbegin(); it != m_pointsToDelete.rend(); ++it)
      {
        int idx = *it;
        pts[idx] = pts.back();
        pts.pop_back();
        szs[idx] = szs.back();
        szs.pop_back();
        m_flags[idx] = m_flags.back();
        m_flags.pop_back();
      }
      m_pointsToDelete.clear();
      break;
    }
    bool trianglesChanged(false);
    trianglesChanged = m_tin->OptimizeTriangulation();

    if (trianglesChanged && RELAXTYPE_SPRING == relaxtype)
    { // set up spring relax if the triangles have changed
      SetupNeighbors();
    }
  }

} // MeRelaxerImpl::Relax
//------------------------------------------------------------------------------
/// \brief Sets the relaxation method
/// \param[in] a_relaxType: Case insensitive string that identifies the type of
/// relaxation. Acceptable values include: "spring_relax" \return true if the
/// relaxtype could be set
//------------------------------------------------------------------------------
bool MeRelaxerImpl::SetRelaxationMethod(const std::string& a_relaxType)
{
  bool rval(false);
  std::string type = stToLowerCopy(a_relaxType);
  if (type == "spring_relaxation")
  {
    m_relaxType = RELAXTYPE_SPRING;
    rval = true;
  }
  return rval;
} // MeRelaxerImpl::SetRelaxationMethod
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

  // Progress
  double iterationPercent = a_iteration / (double)a_numiterations;

  // Relax the marked nodes
  VecPt3d& points = m_tin->Points(); // for convenience
  VecInt2d& trisAdjToPts = m_tin->TrisAdjToPts();
  size_t nPoints = points.size();
  Pt3d origlocation, newlocation;
  for (size_t p = 0; p < nPoints; ++p)
  {
    if (m_flags[p])
    {
      origlocation = points[p];
      // do general point relax
      switch (a_relaxType)
      {
      case RELAXTYPE_AREA:
        AreaRelax((int)p, newlocation);
        break;
      case RELAXTYPE_ANGLE:
        AngleRelax((int)p, newlocation);
        break;
      case RELAXTYPE_SPRING:
        SpringRelaxSinglePoint((int)p, newlocation);
        break;
      default:
        XM_ASSERT(false);
        break;
      }

      // Change the points location and check for bad triangles
      if (NewLocationIsValid(p, newlocation))
      {
        points[p].x = newlocation.x;
        points[p].y = newlocation.y;
        if (m_sizer)
        {
          // point moved and we must update its target size when doing
          // spring_relax
          m_pointSizes[p] = m_sizer->SizeFromLocation(points[p]);
        }
      }
      else if (RELAXTYPE_SPRING == a_relaxType)
      { // mark point for deletion
        m_pointsToDelete.push_back(static_cast<int>(p));
      }

    } // if (m_flags[p])

  } // for (int p = 0; p < nPoints; ++p)
} // MeRelaxerImpl::RelaxMarkedPoints
//------------------------------------------------------------------------------
/// \brief Relax a point using area of surrounding triangles. Trys to move
///        point to the center of the area. Compare to rliAreaRelax.
/// \param[in] a_point: The point index to be relaxed.
/// \param[out] a_newLocation: The location after point is relaxed (moved).
//------------------------------------------------------------------------------
void MeRelaxerImpl::AreaRelax(int a_point, Pt3d& a_newLocation)
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

  a_newLocation.x = sumx / sumarea;
  a_newLocation.y = sumy / sumarea;
} // MeRelaxerImpl::AreaRelax
//------------------------------------------------------------------------------
/// \brief Relax a point using angles. Moves point to equalize the angles
///        between the edges touching the point. Compare to rliAngleRelax.
/// \param[in] a_point: The point index to be relaxed.
/// \param[out] a_newLocation: The location after point is relaxed (moved).
//------------------------------------------------------------------------------
void MeRelaxerImpl::AngleRelax(int a_point, Pt3d& a_newLocation)
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
} // MeRelaxerImpl::AngleRelax
//------------------------------------------------------------------------------
/// \brief Relax a point using spring method
/// \param[in] a_point: The point index to be relaxed.
/// \param[out] a_newLocation: The location after point is relaxed (moved).
/// \return true if there were no errors.
//------------------------------------------------------------------------------
void MeRelaxerImpl::SpringRelaxSinglePoint(int a_point, Pt3d& a_newLocation)
{
  double fx(0.0);
  double fy(0.0);

  double springStiffness(1.0);
  double stepSize(0.1);
  VecPt3d& points = m_tin->Points();
  VecInt& neighbors(m_pointNeighbors[a_point]);
  double numPtsFactor(1.7 / (double)neighbors.size());
  Pt3d& p0(points[a_point]);
  double size0(m_pointSizes[a_point]);
  for (size_t i = 0; i < neighbors.size(); ++i)
  {
    int neighborIdx = neighbors[i];
    Pt3d& p1(points[neighborIdx]);
    double size1(m_pointSizes[neighborIdx]);
    Pt3d vec = gmCreateVector(p0, p1);
    double dist = Mdist(p0.x, p0.y, p1.x, p1.y);
    double targetDist = 0.5 * (size0 + size1); /// \todo: * this->get_adjust_factor(p0, p1)
    double force(springStiffness * (dist - targetDist));
    double _fx = force * vec.x / dist;
    double _fy = force * vec.y / dist;
    fx += _fx;
    fy += _fy;
  }
  fx *= numPtsFactor;
  fy *= numPtsFactor;
  a_newLocation.x = p0.x + fx;
  a_newLocation.y = p0.y + fy;
  a_newLocation.z = p0.z;
} // MeRelaxerImpl::SpringRelax
//------------------------------------------------------------------------------
/// \brief Set up neighbor point information to be used by the
/// spring relax algorithm
//------------------------------------------------------------------------------
void MeRelaxerImpl::SetupNeighbors()
{
  m_pointNeighbors.clear();
  VecPt3d& points = m_tin->Points();
  size_t nPts(points.size());
  if (m_flags.empty())
  {
    m_flags.assign(nPts, RELAX_RELAX);
  }
  m_pointNeighbors.assign(nPts, VecInt());
  VecInt& tris = m_tin->Triangles();
  VecInt2d& adjTris2d = m_tin->TrisAdjToPts();
  for (size_t i = 0; i < points.size(); ++i)
  {
    if (m_flags[i] == 0) // We don't relax fixed points, so we don't need its neighbors.
    {
      continue;
    }
    const VecInt& adjTris = adjTris2d[i];
    size_t numTri = adjTris.size();
    // find the other points that make up the edges that a_point is connected to
    std::set<int> setPoints;
    for (size_t j = 0; j < numTri; ++j)
    {
      int idx = adjTris[j] * 3;
      for (int t = 0; t < 3; ++t)
      {
        int ptIdx = tris[idx + t];
        XM_ASSERT(ptIdx < nPts);
        setPoints.insert(tris[idx + t]);
      }
    }
    setPoints.erase(static_cast<int>(i));
    m_pointNeighbors[i] = VecInt(setPoints.begin(), setPoints.end());
  }
} // MeRelaxerImpl::SetupNeighbors()
//------------------------------------------------------------------------------
/// \brief Set up point sizes to be used by the
/// spring relax algorithm
//------------------------------------------------------------------------------
void MeRelaxerImpl::SetupPointSizes()
{
  m_pointSizes.clear();
  VecPt3d& points = m_tin->Points();
  size_t nPts(points.size());
  m_pointSizes.reserve(nPts);
  for (size_t i = 0; i < points.size(); ++i)
  {
    // compute point sizes prior to relaxing below
    Pt3d& p(m_tin->Points()[i]);
    m_pointSizes.push_back(m_sizer->SizeFromLocation(p));
  }
} // MeRelaxerImpl::SetupPointSizes
//------------------------------------------------------------------------------
/// \brief Checks if a new relaxed location results in valid triangle (area is positive)
/// \param[in] a_idx: the index to the point
/// \param[in] a_newLocation: x,y,z new coordinates
/// \return true if the new location will give valid triangles
//------------------------------------------------------------------------------
bool MeRelaxerImpl::NewLocationIsValid(size_t a_idx, Pt3d& a_newLocation)
{
  Pt3d originalLocation = m_tin->Points()[a_idx];
  m_tin->Points()[a_idx] = a_newLocation;
  const VecInt& adjTris = m_tin->TrisAdjToPts()[a_idx];
  bool rval(true);
  for (size_t i = 0; rval && i < adjTris.size(); ++i)
  {
    // Is triangle ill-formed? See that it has some area
    if (!GT_EPS(m_tin->TriangleArea(adjTris[i]), 0.0, FLT_EPSILON))
    {
      rval = false;
    }
  }
  m_tin->Points()[a_idx] = originalLocation;
  return rval;
} // MeRelaxerImpl::NewLocationIsValid
//------------------------------------------------------------------------------
/// \brief Checks if all triangles have a positive area
/// \param[in] a_tin: the tin
/// \return true if the triangles all have positive area
//------------------------------------------------------------------------------
bool MeRelaxerImpl::AllTrianglesHavePositiveArea(BSHP<TrTin> a_tin)
{
  for (size_t i = 0; i < a_tin->NumTriangles(); ++i)
  {
    if (a_tin->TriangleArea(static_cast<int>(i)) < 0.0)
      return false;
  }
  return true;
} // MeRelaxerImpl::AllTrianglesHavePositiveArea

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
  // MeRelaxerImpl->SetBSHP(polyMesher); // If MeRelaxerImpl needs ptr to
  // MeRelaxer
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
} // namespace

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
  double pointsAfterA[] = {
    0.0, 0.0,  10.0,      0.0,       20.0,      0.0,       30.0,      0.0,       40.0, 0.0,
    0.0, 10.0, 11.095340, 8.8783699, 20.021946, 10.198352, 29.943135, 9.9222857, 40.0, 10.0,
    0.0, 20.0, 9.3320588, 19.085690, 20.677562, 20.864280, 31.148854, 18.861555, 40.0, 20.0,
    0.0, 30.0, 9.9438849, 29.753838, 18.833441, 31.188921, 29.137598, 29.132561, 40.0, 30.0,
    0.0, 40.0, 10.0,      40.0,      20.0,      40.0,      30.0,      40.0,      40.0, 40.0};
  VecPt3d pointsAfter = iArrayToVecPt3d(pointsAfterA, XM_COUNTOF(pointsAfterA));
  TS_ASSERT_DELTA_VECPT3D(pointsAfter, tin->Points(), kDelta);

} // MeRelaxerUnitTests::testRelaxWhileMeshing
//------------------------------------------------------------------------------
/// \brief Tests the spring relax setup. Uses the TIN below
/// \code
///
///  x is a missing triangle
///
///  6-9--7----8
///  | // |  / |
///  3----4----5
///  | /  |  \ |
///  0----1----2
///
/// \endcode
//------------------------------------------------------------------------------
void MeRelaxerUnitTests::testSpringRelaxSetup()
{
  VecPt3d points = {{-10, -10, 0}, {0, -10, 0},  {10, -10, 0}, {-10, 0, 0}, {0, 0, 0},
                    {10, 0, 0},    {-10, 10, 0}, {0, 10, 0},   {10, 10, 0}, {-5, 10, 0}};
  VecInt tris = {0, 4, 3, 0, 1, 4, 1, 2, 4, 2, 5, 4, 3, 7, 9, 4, 8, 7, 4, 5, 8, 3, 9, 6, 3, 4, 7};
  // create a test tin
  BSHP<TrTin> tin = TrTin::New();
  tin->Points() = points;
  tin->Triangles() = tris;
  tin->BuildTrisAdjToPts();
  // create a sizer that returns a constant value
  BSHP<MePolyRedistributePts> sizer = MePolyRedistributePts::New();
  sizer->SetConstantSizeFunc(3.0);

  MeRelaxerImpl r;
  r.m_tin = tin;
  r.m_sizer = sizer;
  r.SetupPointSizes();
  r.SetupNeighbors();

  VecDbl expectedPointSizes(points.size(), 3);
  TS_ASSERT_EQUALS_VEC(expectedPointSizes, r.m_pointSizes);
  VecInt2d expectedPointNeighbors = {
    {1, 3, 4}, {0, 2, 4}, {1, 4, 5},    {0, 4, 6, 7, 9}, {0, 1, 2, 3, 5, 7, 8},
    {2, 4, 8}, {3, 9},    {3, 4, 8, 9}, {4, 5, 7},       {3, 6, 7}};
  TS_ASSERT_EQUALS_VEC2D(expectedPointNeighbors, r.m_pointNeighbors);
} // MeRelaxerUnitTests::testSpringRelaxSetup
//------------------------------------------------------------------------------
/// \brief Tests the spring relax setup. Uses the TIN below. This tin has
/// missing triangles compared to the previous example \code
///
///  x is a missing triangle
///
///  6-9  7----8
///  |/x/x|  / |
///  3----4----5
///  | /  |  \ |
///  0----1----2
///
/// \endcode
//------------------------------------------------------------------------------
void MeRelaxerUnitTests::testSpringRelaxSetup2()
{
  VecPt3d points = {{-10, -10, 0}, {0, -10, 0},  {10, -10, 0}, {-10, 0, 0}, {0, 0, 0},
                    {10, 0, 0},    {-10, 10, 0}, {0, 10, 0},   {10, 10, 0}, {-5, 10, 0}};
  VecInt tris = {0, 4, 3, 0, 1, 4, 1, 2, 4, 2, 5, 4, 4, 8, 7, 4, 5, 8, 3, 9, 6};
  // create a test tin
  BSHP<TrTin> tin = TrTin::New();
  tin->Points() = points;
  tin->Triangles() = tris;
  tin->BuildTrisAdjToPts();
  // create a sizer that returns a constant value
  BSHP<MePolyRedistributePts> sizer = MePolyRedistributePts::New();
  sizer->SetConstantSizeFunc(3.0);

  MeRelaxerImpl r;
  r.m_tin = tin;
  r.m_sizer = sizer;
  r.SetupPointSizes();
  r.SetupNeighbors();

  VecDbl expectedPointSizes(points.size(), 3);
  TS_ASSERT_EQUALS_VEC(expectedPointSizes, r.m_pointSizes);
  VecInt2d expectedPointNeighbors = {
    {1, 3, 4}, {0, 2, 4}, {1, 4, 5}, {0, 4, 6, 9}, {0, 1, 2, 3, 5, 7, 8},
    {2, 4, 8}, {3, 9},    {4, 8},    {4, 5, 7},    {3, 6}};
  TS_ASSERT_EQUALS_VEC2D(expectedPointNeighbors, r.m_pointNeighbors);
} // MeRelaxerUnitTests::testSpringRelaxSetup
//------------------------------------------------------------------------------
/// \brief Tests the spring relax of a point. Uses the TIN below.
/// \code
///
///  x is a missing triangle
///
///  6----7----8
///  | \  |  / |
///  3-------4-5
///  | /  |  \ |
///  0----1----2
///
/// \endcode
//------------------------------------------------------------------------------
void MeRelaxerUnitTests::testSpringRelaxSinglePoint()
{
  VecPt3d points = {{-10, -10, 1}, {0, -10, 2},  {10, -10, 3}, {-10, 0, 4}, {8, 7, 5},
                    {10, 0, 6},    {-10, 10, 7}, {0, 10, 8},   {10, 10, 9}};
  VecInt tris = {0, 4, 3, 0, 1, 4, 1, 2, 4, 2, 5, 4, 3, 4, 6, 4, 7, 6, 4, 5, 8, 4, 8, 7};
  // create a test tin
  BSHP<TrTin> tin = TrTin::New();
  tin->Points() = points;
  tin->Triangles() = tris;
  tin->BuildTrisAdjToPts();
  // create a sizer that returns a constant value
  BSHP<MePolyRedistributePts> sizer = MePolyRedistributePts::New();
  sizer->SetConstantSizeFunc(10.0);

  MeRelaxerImpl r;
  r.m_tin = tin;
  r.m_sizer = sizer;
  r.SetupPointSizes();
  r.SetupNeighbors();

  Pt3d newloc, startPt(8, 7, 0), ptZero;
  double startDist(Mdist(0.0, 0.0, startPt.x, startPt.y));
  r.SpringRelaxSinglePoint(4, newloc);
  double newDist(Mdist(0.0, 0.0, newloc.x, newloc.y));
  TS_ASSERT(newDist < startDist);
  TS_ASSERT(r.NewLocationIsValid(4, newloc));

  Pt3d expectedNewLoc(0.905, 0.542, 5);
  TS_ASSERT_DELTA_PT3D(expectedNewLoc, newloc, 1e-2);
  TS_ASSERT(r.NewLocationIsValid(4, newloc));

  tin->Points()[4] = newloc;
  startDist = newDist;
  r.SpringRelaxSinglePoint(4, newloc);
  newDist = Mdist(0.0, 0.0, newloc.x, newloc.y);
  TS_ASSERT(newDist < startDist);
  expectedNewLoc = Pt3d(0.023, 0.016, 5);
  TS_ASSERT_DELTA_PT3D(expectedNewLoc, newloc, 1e-2);

  tin->Points()[4] = newloc;
  startDist = newDist;
  r.SpringRelaxSinglePoint(4, newloc);
  newDist = Mdist(0.0, 0.0, newloc.x, newloc.y);
  TS_ASSERT(newDist < startDist);
  expectedNewLoc = Pt3d(0, 0, 5);
  TS_ASSERT_DELTA_PT3D(expectedNewLoc, newloc, 1e-2);
  TS_ASSERT(r.NewLocationIsValid(4, newloc));
} // MeRelaxerUnitTests::testSpringRelaxSinglePoint
//------------------------------------------------------------------------------
/// \brief Tests the spring relax of a point. Uses the TIN below.
/// \code
///
///  x is a missing triangle
///
///  6----7----8
///  | \  |  /   \
///  3-------4-----5
///  | /  |  \   /
///  0----1----2
///
/// \endcode
//------------------------------------------------------------------------------
void MeRelaxerUnitTests::testSpringRelaxSinglePoint2()
{
  VecPt3d points = {{-10, -10, 0}, {0, -10, 0},  {10, -10, 0}, {-10, 0, 0}, {8, 7, 0},
                    {15, 0, 0},    {-10, 10, 0}, {0, 10, 0},   {10, 10, 0}};
  VecInt tris = {0, 4, 3, 0, 1, 4, 1, 2, 4, 2, 5, 4, 3, 4, 6, 4, 7, 6, 4, 5, 8, 4, 8, 7};
  // create a test tin
  BSHP<TrTin> tin = TrTin::New();
  tin->Points() = points;
  tin->Triangles() = tris;
  tin->BuildTrisAdjToPts();
  // create a sizer that returns a constant value
  BSHP<MePolyRedistributePts> sizer = MePolyRedistributePts::New();
  sizer->SetConstantSizeFunc(10.0);

  MeRelaxerImpl r;
  r.m_tin = tin;
  r.m_sizer = sizer;
  r.SetupPointSizes();
  r.SetupNeighbors();
  r.m_pointSizes[5] = 15;

  Pt3d newloc, startPt(8, 7, 0), ptZero;
  double startDist(Mdist(0.0, 0.0, startPt.x, startPt.y));
  r.SpringRelaxSinglePoint(4, newloc);
  double newDist(Mdist(0.0, 0.0, newloc.x, newloc.y));
  TS_ASSERT(newDist < startDist);
  TS_ASSERT(r.NewLocationIsValid(4, newloc));

  Pt3d expectedNewLoc(0.673, 0.377, 0);
  TS_ASSERT_DELTA_PT3D(expectedNewLoc, newloc, 1e-2);
  TS_ASSERT(r.NewLocationIsValid(4, newloc));

  tin->Points()[4] = newloc;
  startDist = newDist;
  r.SpringRelaxSinglePoint(4, newloc);
  newDist = Mdist(0.0, 0.0, newloc.x, newloc.y);
  TS_ASSERT(newDist < startDist);
  expectedNewLoc = Pt3d(0.547, -0.005, 0);
  TS_ASSERT_DELTA_PT3D(expectedNewLoc, newloc, 1e-2);
  TS_ASSERT(r.NewLocationIsValid(4, newloc));

  tin->Points()[4] = newloc;
  startDist = newDist;
  r.SpringRelaxSinglePoint(4, newloc);
  newDist = Mdist(0.0, 0.0, newloc.x, newloc.y);
  TS_ASSERT(newDist < startDist);
  expectedNewLoc = Pt3d(0.545, 0, 0);
  TS_ASSERT_DELTA_PT3D(expectedNewLoc, newloc, 1e-2);
  TS_ASSERT(r.NewLocationIsValid(4, newloc));

} // MeRelaxerUnitTests::testSpringRelaxSinglePoint2
//------------------------------------------------------------------------------
/// \brief Tests the spring relax of a point. Uses the TIN below. We are relaxing
/// point 4 and it has 3 connections to the upper edge, 2 connection to the middle
/// and 1 connection to the bottom edge. The result is that the points is "pulled"
/// toward the upper edge. This is not desirable for pre quad mesh generation.
/// \code
///
///  x is a missing triangle
///
///  6----7----8
///  | \  | /  |
///  3------4--5
///  | \  |  / |
///  0----1----2
///
/// \endcode
//------------------------------------------------------------------------------
void MeRelaxerUnitTests::testSpringRelaxSinglePoint3()
{
  VecPt3d points = {{-10, -10, 0}, {0, -10, 0},  {10, -10, 0}, {-10, 0, 0}, {8, 7, 0},
                    {10, 0, 0},    {-10, 10, 0}, {0, 10, 0},   {10, 10, 0}};
  VecInt tris = {0, 1, 3, 1, 4, 3, 1, 5, 4, 1, 2, 5, 3, 4, 6, 4, 7, 6, 4, 5, 8, 4, 8, 7};
  // create a test tin
  BSHP<TrTin> tin = TrTin::New();
  tin->Points() = points;
  tin->Triangles() = tris;
  tin->BuildTrisAdjToPts();
  // create a sizer that returns a constant value
  BSHP<MePolyRedistributePts> sizer = MePolyRedistributePts::New();
  sizer->SetConstantSizeFunc(100.0);

  MeRelaxerImpl r;
  r.m_tin = tin;
  r.m_sizer = sizer;
  r.SetupPointSizes();
  r.SetupNeighbors();

  Pt3d newloc, startPt(8, 7, 0), ptZero;
  double startDist(Mdist(0.0, 0.0, startPt.x, startPt.y));
  r.SpringRelaxSinglePoint(4, newloc);
  TS_ASSERT(!r.NewLocationIsValid(4, newloc));
} // MeRelaxerUnitTests::testSpringRelaxSinglePoint3
//------------------------------------------------------------------------------
/// \brief Tests that a new location for a point is within one of the triangles
/// adjacent to that point's current position.
/// \code
///
///  x is a missing triangle
///
///  6----7----8
///  | /  | \ /
///  3----4--5
///  | /  | \ \ 
///  0----1----2
///
/// \endcode
//------------------------------------------------------------------------------
void MeRelaxerUnitTests::testNewLocationIsValid()
{
  VecPt3d points = {{-10, -10, 0}, {0, -10, 0},  {10, -10, 0}, {-10, 0, 0}, {0, 0, 0},
                    {5, 0, 0},     {-10, 10, 0}, {0, 10, 0},   {10, 10, 0}};
  VecInt tris = {0, 4, 3, 0, 1, 4, 1, 2, 4, 2, 5, 4, 4, 5, 7, 5, 8, 7, 3, 4, 7, 3, 7, 6};
  // create a test tin
  BSHP<TrTin> tin = TrTin::New();
  tin->Points() = points;
  tin->Triangles() = tris;
  tin->BuildTrisAdjToPts();
  // create a sizer that returns a constant value
  BSHP<MePolyRedistributePts> sizer = MePolyRedistributePts::New();
  sizer->SetConstantSizeFunc(3.0);

  MeRelaxerImpl r;
  r.m_tin = tin;
  r.m_sizer = sizer;
  r.SetupPointSizes();
  r.SetupNeighbors();

  Pt3d newloc(-2, -2, 0), ptZero;
  TS_ASSERT(r.NewLocationIsValid(4, newloc));

  tin->Points()[4] = ptZero;
  newloc = Pt3d(-2, 2, 0);
  TS_ASSERT(r.NewLocationIsValid(4, newloc));

  tin->Points()[4] = ptZero;
  newloc = Pt3d(-7, 7, 0);
  TS_ASSERT(!r.NewLocationIsValid(4, newloc));

  tin->Points()[4] = ptZero;
  newloc = Pt3d(7, 0, 0);
  TS_ASSERT(!r.NewLocationIsValid(4, newloc));

  tin->Points()[4] = ptZero;
  newloc = Pt3d(5, 5, 0);
  TS_ASSERT(!r.NewLocationIsValid(4, newloc));

  tin->Points()[4] = ptZero;
  newloc = Pt3d(-10, 0, 0);
  TS_ASSERT(!r.NewLocationIsValid(4, newloc));

  tin->Points()[4] = ptZero;
  newloc = Pt3d(-11, 0, 0);
  TS_ASSERT(!r.NewLocationIsValid(4, newloc));

  tin->Points()[4] = ptZero;
  newloc = Pt3d(-2, -11, 0);
  TS_ASSERT(!r.NewLocationIsValid(4, newloc));

} // MeRelaxerUnitTests::testNewLocationIsValid
//------------------------------------------------------------------------------
/// \brief Tests the spring relax of a point. Uses the TIN below.
/// \code
///
///  x is a missing triangle
///
///  6----7----8
///  | \  |  / |
///  3----4----5
///  | /  |  \ |
///  0----1----2
///
/// \endcode
//------------------------------------------------------------------------------
void MeRelaxerUnitTests::testAllTrianglesHavePositiveArea()
{
  VecPt3d points = {{-10, -10, 0}, {0, -10, 0},  {10, -10, 0}, {-10, 0, 0}, {8, 7, 0},
                    {10, 0, 0},    {-10, 10, 0}, {0, 10, 0},   {10, 10, 0}};
  VecInt tris = {0, 4, 3, 0, 1, 4, 1, 2, 4, 2, 5, 4, 3, 4, 6, 4, 7, 6, 4, 5, 8, 4, 8, 7};
  // create a test tin
  BSHP<TrTin> tin = TrTin::New();
  tin->Points() = points;
  tin->Triangles() = tris;
  tin->BuildTrisAdjToPts();

  MeRelaxerImpl r;
  TS_ASSERT(r.AllTrianglesHavePositiveArea(tin));
  // create an invalid triangle
  tris[1] = 3;
  tris[2] = 4;
  tin->Triangles() = tris;
  TS_ASSERT(!r.AllTrianglesHavePositiveArea(tin));
} // MeRelaxerUnitTests::testAllTrianglesHavePositiveArea

//} // namespace xms
#endif // CXX_TEST
