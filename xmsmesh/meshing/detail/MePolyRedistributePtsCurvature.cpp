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
#include <xmsmesh/meshing/detail/MePolyRedistributePtsCurvature.h>

// 3. Standard library headers
#include <map>

// 4. External library headers

// 5. Shared code headers
#include <xmscore/math/math.h>
#include <xmscore/misc/XmError.h>
#include <xmscore/points/pt.h>
#include <xmscore/stl/vector.h>
#include <xmsinterp/geometry/geoms.h>

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

/// \brief Redistributes the point locations on a polyline or polygon based on curvature
class MePolyRedistributePtsCurvatureImpl : public MePolyRedistributePtsCurvature
{
public:
  MePolyRedistributePtsCurvatureImpl() {}

  virtual VecPt3d Redistribute(const VecPt3d& a_points,
                               double a_featureSize,
                               double a_mean_spacing,
                               double a_minimumCurvature = 0.001,
                               bool a_smooth = false);
  void Setup(const VecPt3d&);
  VecPt3d PlacePoints(double a_featureSize,
                      int a_numPoints,
                      double a_minimumCurvature,
                      bool a_smooth);
  void GetSignificantPoints(double a_featureSize);
  void CalculateCurvature(double a_featureSize, double a_minimumCurvature);
  double GetCurvatureFromParameter(double a_param, double a_interval);
  Pt3d GetPointFromParameter(double a_param);
  void GetParameterIFM(double a_param, double a_interval, double& a_ti, double& a_tm, double& a_tf);
  void ShiftAndAggregateOpen();
  void ShiftAndAggregateClosed();
  void DoSmoothing(bool a_smooth);
  VecPt3d NewPointsFromParamCurvs(int a_numPoints);

  VecPt3d m_points;                   ///< locations defining polyline/polygon
  VecDbl m_segmentLengths;            ///< length of each segment
  VecDbl m_accumulatedSegmentLengths; ///< accumulated length of polyline
  VecDbl m_parametricDistance;        ///< distance (0-1) from starting to ending point
  VecDbl m_curvature;                 ///< curvature at each point
  double m_length = 0.0;              ///< total length
  bool m_open = false;                ///< false means polygon, true mean polyline
  double m_tol = 1e-6;                ///< tolerance used for geometric calculations
  std::map<double, int> m_distMap;    ///< map of distances and indices
};
//------------------------------------------------------------------------------
/// \brief Creates an instance of this class
/// \return MePolyRedistributePts
//------------------------------------------------------------------------------
BSHP<MePolyRedistributePtsCurvature> MePolyRedistributePtsCurvature::New()
{
  BSHP<MePolyRedistributePtsCurvature> ret(new MePolyRedistributePtsCurvatureImpl());
  return ret;
} // MePolyRedistributePtsCurvature::New
//------------------------------------------------------------------------------
/// \brief destructor
//------------------------------------------------------------------------------
MePolyRedistributePtsCurvature::~MePolyRedistributePtsCurvature()
{
} // MePolyRedistributePts::~MePolyRedistributePts
//------------------------------------------------------------------------------
/// \brief constructor
//------------------------------------------------------------------------------
MePolyRedistributePtsCurvature::MePolyRedistributePtsCurvature()
{
} // MePolyRedistributePtsCurvature::MePolyRedistributePtsCurvature

//------------------------------------------------------------------------------
/// \brief Redistribute points according to curvature
/// \param[in] a_points: Points defining a closed polygon (if last point is the same as the first)
///   or open polyline.
/// \param[in] a_featureSize: The size of the smallest feature in the polyline to be detected.
///   Large values will generate point distributions that follow coarser curvatures.
/// \param[in] a_meanSpacing: The mean spacing between the distributed points.
/// \param[in] a_minimumCurvature: The value of the curvature to be used instead of 0 in straight
///   lines. It limits the maximum spacing between points. If not included, the default is 0.001.
/// \param[in] a_smooth: Detemines if the curvatures are to be averaged by a rolling 0.25-0.5-0.25
///   weighted rolling average.
/// \return a vector of the redistributed points.
//------------------------------------------------------------------------------
VecPt3d MePolyRedistributePtsCurvatureImpl::Redistribute(const VecPt3d& a_points,
                                                         double a_featureSize,
                                                         double a_meanSpacing,
                                                         double a_minimumCurvature,
                                                         bool a_smooth)
{
  XM_ENSURE_TRUE(!a_points.empty(), VecPt3d());

  Setup(a_points);
  int numPoints = int(m_length / a_meanSpacing); /// \todo:  Guarantee an even number of segments?
  // if open, numPoints should be odd for an even number of segments; if  clossed, it should be
  // even.
  return PlacePoints(a_featureSize, numPoints, a_minimumCurvature, a_smooth);
} // MePolyRedistributePtsCurvatureImpl::Redistribute
//------------------------------------------------------------------------------
/// \brief sets up the class to do a Redistribute operation
/// \param[in] a_points: The locations of the input polyline or polygon
//------------------------------------------------------------------------------
void MePolyRedistributePtsCurvatureImpl::Setup(const VecPt3d& a_points)
{
  m_parametricDistance.clear();
  m_curvature.clear();

  m_points = a_points;
  Pt3d pMin, pMax;
  gmExtents2D(m_points, pMin, pMax);
  m_tol = Mdist(pMin.x, pMin.y, pMax.x, pMax.y) * 1e-9;
  m_open = m_points.back() != m_points.front();
  size_t nLengths = m_open ? m_points.size() - 1 : m_points.size();
  m_segmentLengths.assign(m_points.size(), 0.0);
  m_accumulatedSegmentLengths.assign(m_points.size(), 0.0);

  double total(0.0);
  for (size_t i = 0; i < m_points.size() - 1; ++i)
  {
    const Pt3d& p0 = a_points[i];
    const Pt3d& p1 = a_points[i + 1];
    double length = Mdist(p0.x, p0.y, p1.x, p1.y);
    m_segmentLengths[i] = length;
    m_accumulatedSegmentLengths[i] = total;
    total += length;
  }
  m_length = total;
} // MePolyRedistributePtsCurvatureImpl::SetUp
//------------------------------------------------------------------------------
/// \brief Redistribute points according to curvature
/// \param[in] a_featureSize: The maximum distance around each point when computing curvature
/// at that point. In the future we may want to change the algorithm to use this distance
/// around each point instead of just using it as a maximum.
/// \param[in] a_numPoints: The number of points to be distributed along the polyline.
/// \param[in] a_minimumCurvature: The value of the curvature to be used instead of 0 in staight
///   lines. It limits the maximum spacing between points. If not included, the default is 0.001.
/// \param[in] a_smooth: Detemines if the curvatures are to be averaged by a rolling 0.25-0.5-0.25
///   weighted rolling average.
/// \return a vector of the redistributed points along the interior of the polyline. It does NOT
///   include the original first and last point.
//------------------------------------------------------------------------------
VecPt3d MePolyRedistributePtsCurvatureImpl::PlacePoints(double a_featureSize,
                                                        int a_numPoints,
                                                        double a_minimumCurvature,
                                                        bool a_smooth)
{
  GetSignificantPoints(a_featureSize);
  CalculateCurvature(a_featureSize, a_minimumCurvature);

  if (m_open)
  {
    ShiftAndAggregateOpen();
  }
  else
  {
    ShiftAndAggregateClosed();
  }
  DoSmoothing(a_smooth);

  return NewPointsFromParamCurvs(a_numPoints);
} // MePolyRedistributePtsCurvatureImpl::PlacePoints
//------------------------------------------------------------------------------
/// \brief Calculates the curvature at the points where it has not yet been
/// calculated.
/// \param[in] a_featureSize: The size of the smallest feature in the polyline to be detected.
///   Large values will generate point distributions that follow coarser curvatures.
/// \param[in] a_minimumCurvature: The value of the curvature to be used instead of 0 in staight
///   lines. It limits the maximum spacing between points. If not included, the default is 0.001.
//------------------------------------------------------------------------------
void MePolyRedistributePtsCurvatureImpl::CalculateCurvature(double a_featureSize,
                                                            double a_minimumCurvature)
{
  m_distMap.clear();
  int i(0);
  for (const auto& d : m_accumulatedSegmentLengths)
  {
    m_distMap.insert(std::make_pair(d, i));
    i++;
  }

  double interval = a_featureSize / m_length;
  for (size_t i = 0; i < m_parametricDistance.size(); ++i)
  {
    double curv(0.0);
    if (m_curvature[i] < 0) // Not calculated yet  .isnan()
    {
      curv = GetCurvatureFromParameter(m_parametricDistance[i], interval);
      // This puts a minimum limit to the curvature
      curv = std::max(a_minimumCurvature, fabs(curv));
      m_curvature[i] = curv;
    }
  }
} // MePolyRedistributePtsCurvatureImpl::CalculateCurvature
//------------------------------------------------------------------------------
/// \brief Redistribute points according to curvature
/// \param[in] a_featureSize: The size of the smallest feature in the polyline to be detected.
///   Large values will generate point distributions that follow coarser curvatures.
//------------------------------------------------------------------------------
void MePolyRedistributePtsCurvatureImpl::GetSignificantPoints(double a_featureSize)
{
  double sumDistance(0.0);
  double total = m_length;
  double halfFeatureSize = a_featureSize / 2.0;
  for (size_t i = 0; i < m_points.size() - 1; ++i)
  {
    double delta_d = m_segmentLengths[i];
    XM_ASSERT(sumDistance == m_accumulatedSegmentLengths[i]);
    m_parametricDistance.push_back(sumDistance / total);
    m_curvature.push_back(-1);
    if (delta_d > a_featureSize) // Add two params half a featureSize from each segment end
    {
      m_parametricDistance.push_back((sumDistance + halfFeatureSize) / total);
      m_curvature.push_back(0);

      m_parametricDistance.push_back((sumDistance + delta_d - halfFeatureSize) / total);
      m_curvature.push_back(0);
    }
    else
    {
      m_parametricDistance.push_back((sumDistance + 0.5 * delta_d) / total);
      m_curvature.push_back(0);
    }
    sumDistance += delta_d;
  }
  XM_ASSERT(sumDistance == total);

  m_parametricDistance.push_back(1);
  m_curvature.push_back(-1);
} // MePolyRedistributePtsCurvatureImpl::GetSignificantPoints
//------------------------------------------------------------------------------
/// \brief Redistribute points according to curvature
/// \param[in] a_param: Point, in parameterized station form, where the curvature will be
///   calculated.
/// \param[in] a_interval: Parameterized form of the feature_size. Determines the two other points
/// used to calculate the curvature at tc
/// \return The computed curvature.
//------------------------------------------------------------------------------
double MePolyRedistributePtsCurvatureImpl::GetCurvatureFromParameter(double a_param,
                                                                     double a_interval)
{
  double ti, tm, tf;
  GetParameterIFM(a_param, a_interval, ti, tm, tf);
  Pt3d pi = GetPointFromParameter(ti);
  Pt3d pm = GetPointFromParameter(tm);
  Pt3d pf = GetPointFromParameter(tf);

  // get curvature
  double xc, yc, r2;
  gmCircumcircleWithTol(&pi, &pm, &pf, &xc, &yc, &r2, m_tol);
  double r = sqrt(r2);
  return 1 / r;
} // MePolyRedistributePtsCurvatureImpl::GetCurvatureFromParameter
//------------------------------------------------------------------------------
/// \brief Get location based on parametric value a_param
/// \param[in] a_param: The parameterized position between [0, 1] along the curve.
/// \return The point for the paramater.
//------------------------------------------------------------------------------
Pt3d MePolyRedistributePtsCurvatureImpl::GetPointFromParameter(double a_param)
{
  int idx(-1);
  double t = std::min(1.0, std::max(0.0, a_param));
  double station = t * m_length;
  auto it = m_distMap.lower_bound(station);
  if (it != m_distMap.begin())
    it--;
  if (it != m_distMap.end())
  {
    idx = it->second;
  }

  if (idx > -1)
  {
    double d0 = m_accumulatedSegmentLengths[idx];
    double d1 = d0 + m_segmentLengths[idx];
    if (d0 <= station && station <= d1)
    {
      double fraction = (station - d0) / (d1 - d0);
      Pt3d& p0 = m_points[idx];
      Pt3d& p1 = m_points[idx + 1];
      Pt3d pt;
      pt.x = p0.x + fraction * (p1.x - p0.x);
      pt.y = p0.y + fraction * (p1.y - p0.y);
      return pt;
    }
  }
  return m_points.back();
} // MePolyRedistributePtsCurvatureImpl::GetPointFromParameter
//------------------------------------------------------------------------------
/// \brief Calculates the parameterized station values of the two point at each side of tc.
/// \param[in] a_param: A parameterized station [0,1] of the central point.
/// \param[in] a_interval: A parameterized form of the feature size. Determines the two other
///   points used to calculate the curvature at tc.
/// \param[out] a_ti: Parameter of the initial point.
/// \param[out] a_tm: Parameter of the mid point.
/// \param[out] a_tf: Parameter of the final point.
//------------------------------------------------------------------------------
void MePolyRedistributePtsCurvatureImpl::GetParameterIFM(double a_param,
                                                         double a_interval,
                                                         double& a_ti,
                                                         double& a_tm,
                                                         double& a_tf)
{
  a_ti = a_param - a_interval;
  a_tf = a_param + a_interval;
  a_tm = a_param;
  if (m_open)
  {
    if (a_ti < 0.0)
    {
      a_ti = 0.0;
      a_tm = 0.5 * (a_ti + a_tf);
    }
    else if (a_tf > 1.0)
    {
      a_tf = 1.0;
      a_tm = 0.5 * (a_ti + a_tf);
    }
  }
  else
  {
    if (a_ti < 0.0)
    {
      a_ti += 1.0;
    }
    if (a_tf > 1.0)
    {
      a_tf -= 1.0;
    }
  }
} // MePolyRedistributePtsCurvatureImpl::GetParameterIFM
//------------------------------------------------------------------------------
/// \brief Shifts half an interval and aggregates the curvature data in an open polyline
//------------------------------------------------------------------------------
void MePolyRedistributePtsCurvatureImpl::ShiftAndAggregateOpen()
{
  VecDbl shiftedParametricDistance(m_parametricDistance.size() + 1);
  VecDbl shiftedCurvature(m_curvature.size() + 1);
  size_t i = 0;

  // shifted is 1 longer than a_paramCurv so this has to happend once:
  double agg_curv = 0;
  for (size_t i = 1; i < shiftedCurvature.size() - 1; ++i)
  {
    double pStart = m_parametricDistance[i - 1];
    agg_curv += m_curvature[i - 1];
    double p(0.5 * (pStart + m_parametricDistance[i]));
    shiftedParametricDistance[i] = p;
    shiftedCurvature[i] = agg_curv;
  }
  shiftedParametricDistance.back() = 1.0;
  shiftedCurvature.back() = agg_curv;
  m_parametricDistance.swap(shiftedParametricDistance);
  m_curvature.swap(shiftedCurvature);
} // MePolyRedistributePtsCurvatureImpl::ShiftAndAggregateOpen
//------------------------------------------------------------------------------
/// \brief Shifts half an interval and aggregates the curvature data in closed polygon
//------------------------------------------------------------------------------
void MePolyRedistributePtsCurvatureImpl::ShiftAndAggregateClosed()
{
  // a is the point before the end.
  // b is the point at the end and the beginning.
  // c is the point after the beginning.
  XM_ASSERT(m_parametricDistance.size() >= 2);
  // XM_ASSERT(a_paramCurvs.back() == a_paramCurvs.front());
  VecDbl& pdist(m_parametricDistance);
  VecDbl& curve(m_curvature);
  int idxA(static_cast<int>(m_parametricDistance.size() - 2)), idxB(0), idxC(1);
  double p_a(pdist[idxA]), p_b(pdist[idxB]), p_c(pdist[idxC]);
  double curv_a(curve[idxA]), curv_b(curve[idxB]), curv_c(curve[idxC]);
  double s_ab(1.0 - p_a);
  double s_bc(p_c); // p_c - 0.0;
  double proportion = s_ab / (s_ab + s_bc);
  // Curvature at m_param = 0 and 1 after shifting.
  double curv_0 = curv_a + proportion * (curv_b - curv_a);

  VecDbl shiftedPdist(pdist.size() + 1);
  VecDbl shiftedCurve(curve.size() + 1);
  shiftedCurve[0] = curv_0;
  for (size_t i = 1; i < shiftedPdist.size() - 1; ++i)
  {
    double p(0.5 * (pdist[i - 1] + pdist[i]));
    shiftedPdist[i] = p;
    shiftedCurve[i] = curve[i - 1];
  }
  shiftedPdist.back() = 1.0;
  shiftedCurve.back() = curv_0;

  // Now Aggregate
  VecDbl aggregatedCurve(shiftedCurve.size());
  // we start with zero curvature
  aggregatedCurve[0] = 0.0;
  // the second point we set the aggregated to be just the delta between 0 and 1
  // because m_param=0 doesn't have a m_curv=0. After the second point we just
  // aggregate by summing. We have not thought of a way to make this more clear.
  double deltaCurvature = shiftedCurve[1] - shiftedCurve[0];
  aggregatedCurve[1] = deltaCurvature;
  double agg_curv(aggregatedCurve[1]);
  for (size_t i = 2; i < shiftedCurve.size(); ++i)
  {
    agg_curv += shiftedCurve[i];
    aggregatedCurve[i] = agg_curv;
  }

  m_parametricDistance.swap(shiftedPdist);
  m_curvature.swap(aggregatedCurve);
} // MePolyRedistributePtsCurvatureImpl::ShiftAndAggregateClosed
//------------------------------------------------------------------------------
/// \brief Shifts half an interval and aggregates the curvature data in closed polygon
/// \param[in] a_smooth: flag to indicate if smoothing should happen
//------------------------------------------------------------------------------
void MePolyRedistributePtsCurvatureImpl::DoSmoothing(bool a_smooth)
{
  if (!a_smooth)
    return;

  // size_t n = a_paramCurvs.size() - 1;
  size_t n = m_curvature.size() - 1;
  VecDbl smoothCurve(m_curvature.size());
  smoothCurve[0] = m_curvature[0];
  for (size_t i = 1; i < n; ++i)
  {
    double smooth_curv = 0.25 * (m_curvature[i - 1] + m_curvature[i + 1]) + 0.5 * m_curvature[i];
    smoothCurve[i] = smooth_curv;
  }
  smoothCurve[n] = m_curvature[n];
  m_curvature.swap(smoothCurve);
} // MePolyRedistributePtsCurvatureImpl::DoSmoothing
//------------------------------------------------------------------------------
/// \brief Shifts half an interval and aggregates the curvature data in closed polygon
/// \param[in] a_numPoints: the number of desired points
/// \return vector of points redistributed according to curvature
//------------------------------------------------------------------------------
VecPt3d MePolyRedistributePtsCurvatureImpl::NewPointsFromParamCurvs(int a_numPoints)
{
  // The height of the accumulated curvature is divided by the desired number of points
  double delta_threshold = m_curvature.back() / (a_numPoints - 1);
  double tol(1e-9), sumDelta(0);
  for (int i = 0; i < a_numPoints - 1; ++i)
  {
    sumDelta += delta_threshold;
  }
  tol = fabs(m_curvature.back() - sumDelta) * 10;
  double threshold(0.0);

  VecPt3d result;
  // size_t n = a_paramCurvs.size() - 1;
  size_t n = m_parametricDistance.size() - 1;
  for (size_t i = 0; i < n; ++i)
  {
    // ParamCurv& pc0 = a_paramCurvs[i];
    // double p0 = pc0.m_param;
    double p0 = m_parametricDistance[i];
    // double c0 = pc0.m_curv;
    double c0 = m_curvature[i];
    // ParamCurv& pc1 = a_paramCurvs[i + 1];
    // double p1 = pc1.m_param;
    double p1 = m_parametricDistance[i + 1];
    // double c1 = pc1.m_curv;
    double c1 = m_curvature[i + 1];
    while (LT_TOL(threshold, c1, tol))
    {
      double t = (threshold - c0) / (c1 - c0);
      double p = p0 + t * (p1 - p0);
      Pt3d pt = GetPointFromParameter(p);
      result.push_back(pt);
      threshold += delta_threshold;
    }
  }
  result.push_back(m_points.back());
  return result;
} // MePolyRedistributePtsCurvatureImpl::NewPointsFromParamCurvs

} // namespace xms

#ifdef CXX_TEST
////////////////////////////////////////////////////////////////////////////////

#include <xmsmesh/meshing/detail/MePolyRedistributePtsCurvature.t.h>

#include <fstream>

#include <xmscore/misc/StringUtil.h>
#include <xmscore/testing/TestTools.h>
#include <xmsinterp/geometry/geoms.h>

using namespace xms;

////////////////////////////////////////////////////////////////////////////////
/// \class MePolyRedistributePtsUnitTests
/// \brief tester for the MePolyRedistributePts class
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief tests creating the class
//------------------------------------------------------------------------------
void MePolyRedistributePtsCurvatureUnitTests::testCreateClass()
{
  BSHP<MePolyRedistributePtsCurvature> b = MePolyRedistributePtsCurvature::New();
  TS_ASSERT(b);
} // MePolyRedistributePtsUnitTests::testCreateClass#endif
//------------------------------------------------------------------------------
/// \brief tests a simple case shown below
/// \code
///
///
///          2       4
///        /   \   /   \
///      1       3      \
///    /                 \
///  0                    5
///
/// \endcode
//------------------------------------------------------------------------------
void MePolyRedistributePtsCurvatureUnitTests::testSetup()
{
  MePolyRedistributePtsCurvatureImpl r;
  VecPt3d pts = {{0, 0, 0}, {5, 5, 0}, {10, 10, 0}, {15, 5, 0}, {20, 10, 0}, {25, 0, 0}};
  r.Setup(pts);

  {
    TS_ASSERT(r.m_parametricDistance.empty());
    TS_ASSERT(r.m_curvature.empty());
    double tol(1e-3);
    VecDbl expectedAccumulatedSegmentLengths = {0.0, 7.071, 14.142, 21.213, 28.284, 0};
    TS_ASSERT_DELTA_VEC(expectedAccumulatedSegmentLengths, r.m_accumulatedSegmentLengths, tol);

    double expectedLength(39.464);
    TS_ASSERT_DELTA(expectedLength, r.m_length, tol);

    TS_ASSERT(true == r.m_open);

    TS_ASSERT(pts == r.m_points);

    VecDbl expectedSegmentLengths = {7.071, 7.071, 7.071, 7.071, 11.180, 0};
    TS_ASSERT_DELTA_VEC(expectedSegmentLengths, r.m_segmentLengths, tol);

    double expectTolerance(2.6925824035672520e-008);
    TS_ASSERT_DELTA(expectTolerance, r.m_tol, 1e-7);
  }

  pts.push_back(pts.front());
  r.Setup(pts);
  { // test with a closed loop
    TS_ASSERT(r.m_parametricDistance.empty());
    TS_ASSERT(r.m_curvature.empty());
    double tol(1e-3);
    VecDbl expectedAccumulatedSegmentLengths = {0.0, 7.071, 14.142, 21.213, 28.284, 39.464, 0};
    TS_ASSERT_DELTA_VEC(expectedAccumulatedSegmentLengths, r.m_accumulatedSegmentLengths, tol);

    double expectedLength(64.464);
    TS_ASSERT_DELTA(expectedLength, r.m_length, tol);

    TS_ASSERT(false == r.m_open);

    TS_ASSERT(pts == r.m_points);

    VecDbl expectedSegmentLengths = {7.071, 7.071, 7.071, 7.071, 11.180, 25.0, 0};
    TS_ASSERT_DELTA_VEC(expectedSegmentLengths, r.m_segmentLengths, tol);

    double expectTolerance(2.6925824035672520e-008);
    TS_ASSERT_DELTA(expectTolerance, r.m_tol, 1e-7);
  }

} // MePolyRedistributePtsCurvatureUnitTests::testSetup
//------------------------------------------------------------------------------
/// \brief tests a simple case shown below
/// \code
///
///
///          2       4
///        /   \   /   \
///      1       3      \
///    /                 \
///  0                    5
///
/// \endcode
//------------------------------------------------------------------------------
void MePolyRedistributePtsCurvatureUnitTests::testGetSignificantPoints()
{
  MePolyRedistributePtsCurvatureImpl r;
  VecPt3d pts = {{0, 0, 0}, {5, 5, 0}, {10, 10, 0}, {15, 5, 0}, {20, 10, 0}, {25, 0, 0}};
  r.Setup(pts);
  double featureSize(3.0);
  r.GetSignificantPoints(featureSize);

  {
    VecDbl expectedParam = {0,     0.038, 0.141, 0.179, 0.217, 0.320, 0.358, 0.396,
                            0.500, 0.538, 0.576, 0.679, 0.717, 0.755, 0.962, 1};
    TS_ASSERT_DELTA_VEC(expectedParam, r.m_parametricDistance, 1e-3);
    VecDbl expectedCurvature = {-1, 0, 0, -1, 0, 0, -1, 0, 0, -1, 0, 0, -1, 0, 0, -1};
    TS_ASSERT_DELTA_VEC(expectedCurvature, r.m_curvature, 1e-12);
  }
} // MePolyRedistributePtsCurvatureUnitTests::testGetSignificantPoints
//------------------------------------------------------------------------------
/// \brief tests a simple case shown below
/// \code
///
///
///          2       4-5
///        /   \   /    \
///      1       3       \
///    /                  \
///  0                     6
///
/// \endcode
//------------------------------------------------------------------------------
void MePolyRedistributePtsCurvatureUnitTests::testCalculateCurvature()
{
  MePolyRedistributePtsCurvatureImpl r;
  VecPt3d pts = {{0, 0, 0},   {5, 5, 0},   {10, 10, 0}, {15, 5, 0},
                 {20, 10, 0}, {21, 10, 0}, {25, 0, 0}};
  r.Setup(pts);
  double featureSize(3.0);
  double minCurvature(0.001);
  r.GetSignificantPoints(featureSize);
  r.CalculateCurvature(featureSize, minCurvature);
  {
    VecDbl expectedParam = {0,     0.037, 0.139, 0.177, 0.214, 0.316, 0.353, 0.391, 0.492,
                            0.530, 0.567, 0.669, 0.706, 0.719, 0.731, 0.769, 0.963, 1};
    TS_ASSERT_DELTA_VEC(expectedParam, r.m_parametricDistance, 1e-3);
    VecDbl expectedCurvature = {0.001, 0.000, 0.000, 0.001, 0.000, 0.000, 0.471, 0.000, 0.000,
                                0.471, 0.000, 0.000, 0.516, 0.000, 0.522, 0.000, 0.000, 0.001};
    TS_ASSERT_DELTA_VEC(expectedCurvature, r.m_curvature, 1e-3);
  }
  // now do a closed loop, make the last point the same as the first
  pts.push_back(pts.front());
  r.Setup(pts);
  r.GetSignificantPoints(featureSize);
  r.CalculateCurvature(featureSize, minCurvature);
  {
    VecDbl expectedParam = {0.000, 0.023, 0.086, 0.109, 0.132, 0.194, 0.217,
                            0.240, 0.303, 0.326, 0.349, 0.412, 0.435, 0.442,
                            0.450, 0.473, 0.593, 0.616, 0.639, 0.977, 1.000};
    TS_ASSERT_DELTA_VEC(expectedParam, r.m_parametricDistance, 1e-3);
    VecDbl expectedCurvature = {0.616, 0.000, 0.000, 0.001, 0.000, 0.000, 0.471,
                                0.000, 0.000, 0.471, 0.000, 0.000, 0.516, 0.000,
                                0.522, 0.000, 0.000, 0.552, 0.000, 0.000, 0.616};
    TS_ASSERT_DELTA_VEC(expectedCurvature, r.m_curvature, 1e-3);
  }
} // MePolyRedistributePtsCurvatureUnitTests::testCalculateCurvature
//------------------------------------------------------------------------------
/// \brief tests a simple case shown below
/// \code
///
///
///          2       4-5
///        /   \   /    \
///      1       3       \
///    /                  \
///  0                     6
///
/// \endcode
//------------------------------------------------------------------------------
void MePolyRedistributePtsCurvatureUnitTests::testShiftAndAggregateOpen()
{
  MePolyRedistributePtsCurvatureImpl r;
  VecPt3d pts = {{0, 0, 0},   {5, 5, 0},   {10, 10, 0}, {15, 5, 0},
                 {20, 10, 0}, {21, 10, 0}, {25, 0, 0}};
  r.Setup(pts);
  double featureSize(3.0);
  double minCurvature(0.001);
  r.GetSignificantPoints(featureSize);
  r.CalculateCurvature(featureSize, minCurvature);
  r.ShiftAndAggregateOpen();
  {
    VecDbl expectedParam = {0.000, 0.019, 0.088, 0.158, 0.195, 0.265, 0.334, 0.372, 0.441, 0.511,
                            0.548, 0.618, 0.687, 0.712, 0.725, 0.750, 0.866, 0.981, 1.000};
    TS_ASSERT_DELTA_VEC(expectedParam, r.m_parametricDistance, 1e-3);
    VecDbl expectedCurvature = {0.000, 0.001, 0.001, 0.001, 0.002, 0.002, 0.002,
                                0.473, 0.473, 0.473, 0.945, 0.945, 0.945, 1.461,
                                1.461, 1.983, 1.983, 1.983, 1.983};
    TS_ASSERT_DELTA_VEC(expectedCurvature, r.m_curvature, 1e-3);
  }
} // MePolyRedistributePtsCurvatureUnitTests::testShiftAndAggregateOpen
//------------------------------------------------------------------------------
/// \brief tests a simple case shown below
/// \code
///
///
///          2       4-5
///        /   \   /    \
///      1       3       \
///    /                  \
///  0                     6
///
/// \endcode
//------------------------------------------------------------------------------
void MePolyRedistributePtsCurvatureUnitTests::testShiftAndAggregateClosed()
{
  MePolyRedistributePtsCurvatureImpl r;
  VecPt3d pts = {{0, 0, 0},   {5, 5, 0},   {10, 10, 0}, {15, 5, 0},
                 {20, 10, 0}, {21, 10, 0}, {25, 0, 0},  {0, 0, 0}};
  r.Setup(pts);
  double featureSize(3.0);
  double minCurvature(0.001);
  r.GetSignificantPoints(featureSize);
  r.CalculateCurvature(featureSize, minCurvature);
  r.ShiftAndAggregateClosed();
  {
    VecDbl expectedParam = {0.000, 0.012, 0.054, 0.097, 0.120, 0.163, 0.206, 0.229,
                            0.272, 0.315, 0.338, 0.380, 0.423, 0.439, 0.446, 0.462,
                            0.533, 0.604, 0.627, 0.808, 0.988, 1.000};
    TS_ASSERT_DELTA_VEC(expectedParam, r.m_parametricDistance, 1e-3);
    VecDbl expectedCurvature = {0.000, 0.308, 0.308, 0.308, 0.309, 0.309, 0.309, 0.780,
                                0.780, 0.780, 1.252, 1.252, 1.252, 1.768, 1.768, 2.290,
                                2.290, 2.290, 2.842, 2.842, 2.842, 3.150};
    TS_ASSERT_DELTA_VEC(expectedCurvature, r.m_curvature, 1e-3);
  }
} // MePolyRedistributePtsCurvatureUnitTests::testShiftAndAggregateClosed
//------------------------------------------------------------------------------
/// \brief tests a simple case shown below
/// \code
///
///
///          2       4-5
///        /   \   /    \
///      1       3       \
///    /                  \
///  0                     6
///
/// \endcode
//------------------------------------------------------------------------------
void MePolyRedistributePtsCurvatureUnitTests::testDoSmoothing()
{
  MePolyRedistributePtsCurvatureImpl r;
  VecPt3d pts = {{0, 0, 0},   {5, 5, 0},   {10, 10, 0}, {15, 5, 0},
                 {20, 10, 0}, {21, 10, 0}, {25, 0, 0}};
  r.Setup(pts);
  double featureSize(3.0);
  double minCurvature(0.001);
  r.GetSignificantPoints(featureSize);
  r.CalculateCurvature(featureSize, minCurvature);
  r.ShiftAndAggregateOpen();
  bool smooth(true);
  r.DoSmoothing(smooth);
  {
    VecDbl expectedParam = {0.000, 0.019, 0.088, 0.158, 0.195, 0.265, 0.334, 0.372, 0.441, 0.511,
                            0.548, 0.618, 0.687, 0.712, 0.725, 0.750, 0.866, 0.981, 1.000};
    TS_ASSERT_DELTA_VEC(expectedParam, r.m_parametricDistance, 1e-3);
    VecDbl expectedCurvature = {0.00000, 0.00075, 0.00100, 0.00125, 0.00175, 0.00200, 0.11985,
                                0.35555, 0.47340, 0.59126, 0.82696, 0.94481, 1.07384, 1.33190,
                                1.59154, 1.85277, 1.98338, 1.98338, 1.98338};
    TS_ASSERT_DELTA_VEC(expectedCurvature, r.m_curvature, 1e-3);
  }
} // MePolyRedistributePtsCurvatureUnitTests::testDoSmoothing
//------------------------------------------------------------------------------
/// \brief tests a simple case shown below
/// \code
///
///
///          2       4-5
///        /   \   /    \
///      1       3       \
///    /                  \
///  0                     6
///
/// \endcode
//------------------------------------------------------------------------------
void MePolyRedistributePtsCurvatureUnitTests::testNewPointsFromParamCurvs()
{
  MePolyRedistributePtsCurvatureImpl r;
  VecPt3d pts = {{0, 0, 0},   {5, 5, 0},   {10, 10, 0}, {15, 5, 0},
                 {20, 10, 0}, {21, 10, 0}, {25, 0, 0}};
  r.Setup(pts);
  double featureSize(3.0);
  double minCurvature(0.001);
  r.GetSignificantPoints(featureSize);
  r.CalculateCurvature(featureSize, minCurvature);
  r.ShiftAndAggregateOpen();
  int numPoints(10);
  {
    VecPt3d outPts = r.NewPointsFromParamCurvs(numPoints);
    VecPt3d expectedPoints{{0, 0, 0},          {9.961, 9.961, 0},   {10.457, 9.543, 0},
                           {14.892, 5.108, 0}, {15.388, 5.388, 0},  {19.685, 9.685, 0},
                           {19.987, 9.987, 0}, {20.906, 10.000, 0}, {21.122, 9.695, 0},
                           {25.000, 0.000, 0}};
    TS_ASSERT_DELTA_VECPT3D(expectedPoints, outPts, 1e-3);
  }
  // now with smoothing
  {
    bool smooth(true);
    r.DoSmoothing(smooth);
    VecPt3d outPts = r.NewPointsFromParamCurvs(numPoints);
    VecPt3d expectedPoints{{0, 0, 0},           {9.922, 9.922, 0},   {11.954, 8.046, 0},
                           {14.784, 5.216, 0},  {16.442, 6.442, 0},  {19.546, 9.546, 0},
                           {20.213, 10.000, 0}, {20.656, 10.000, 0}, {21.151, 9.623, 0},
                           {25.000, 0.000, 0}};
    TS_ASSERT_DELTA_VECPT3D(expectedPoints, outPts, 1e-3);
  }

  // now do a closed loop
  pts.push_back(pts.front());
  r.Setup(pts);
  r.GetSignificantPoints(featureSize);
  r.CalculateCurvature(featureSize, minCurvature);
  r.ShiftAndAggregateClosed();
  {
    VecPt3d outPts = r.NewPointsFromParamCurvs(numPoints);
    VecPt3d expectedPoints{{0, 0, 0},          {9.562, 9.562, 0},  {10.350, 9.650, 0},
                           {15.077, 5.077, 0}, {19.673, 9.673, 0}, {20.216, 10.000, 0},
                           {21.143, 9.641, 0}, {24.883, 0.293, 0}, {24.364, 0.000, 0},
                           {0, 0, 0}};
    TS_ASSERT_DELTA_VECPT3D(expectedPoints, outPts, 1e-3);
  }
  // now with smoothing
  {
    bool smooth(true);
    r.DoSmoothing(smooth);
    VecPt3d outPts = r.NewPointsFromParamCurvs(numPoints);
    VecPt3d expectedPoints{{0, 0, 0},          {8.187, 8.187, 0},  {11.158, 8.842, 0},
                           {15.153, 5.153, 0}, {19.523, 9.523, 0}, {20.464, 10.000, 0},
                           {21.194, 9.515, 0}, {24.766, 0.586, 0}, {16.082, 0.000, 0},
                           {0, 0, 0}};
    TS_ASSERT_DELTA_VECPT3D(expectedPoints, outPts, 1e-3);
  }
} // MePolyRedistributePtsCurvatureUnitTests::testNewPointsFromParamCurvs
//------------------------------------------------------------------------------
/// \brief Tests redistribution along a coastline
//------------------------------------------------------------------------------
void MePolyRedistributePtsCurvatureIntermediateTests::testCoastline()
{
  std::string path(std::string(XMS_TEST_PATH) + "redistribution/");
  std::string infile(path + "Coastline.txt"), outFile(path + "Coastline_out.txt"),
    baseFile(path + "Coastline_base.txt");
  std::fstream is;
  is.open(infile, std::fstream::in);
  VecPt3d pts;
  while (is.good())
  {
    Pt3d p;
    is >> p.x >> p.y >> p.z;
    if (is.good())
      pts.push_back(p);
  }

  MePolyRedistributePtsCurvatureImpl r;
  double featureSize(200.0), meanSpacing(50.0), minimumCurvature(.001);
  bool smoothCurvature(true);
  pts = r.Redistribute(pts, featureSize, meanSpacing, minimumCurvature, smoothCurvature);
  {
    int flag = STR_FULLWIDTH;
    int width = 9;
    std::fstream os(outFile, std::fstream::out);
    for (auto& p : pts)
    {
      os << STRstd(p.x, -1, width, flag) << " " << STRstd(p.y, -1, width, flag) << "\n";
    }
  }
  TS_ASSERT_TXT_FILES_EQUAL(baseFile, outFile);
} // MePolyRedistributePtsCurvatureIntermediateTests::testCoastline
//------------------------------------------------------------------------------
/// \brief Tests redistribution along an island
//------------------------------------------------------------------------------
void MePolyRedistributePtsCurvatureIntermediateTests::testIsland()
{
  std::string path(std::string(XMS_TEST_PATH) + "redistribution/");
  std::string infile(path + "Island.txt"), outFile(path + "Island_out.txt"),
    baseFile(path + "Island_base.txt");
  std::fstream is;
  is.open(infile, std::fstream::in);
  VecPt3d pts;
  while (is.good())
  {
    Pt3d p;
    is >> p.x >> p.y >> p.z;
    if (is.good())
      pts.push_back(p);
  }

  MePolyRedistributePtsCurvatureImpl r;
  double featureSize(200.0), meanSpacing(50.0), minimumCurvature(.001);
  bool smoothCurvature(true);
  pts = r.Redistribute(pts, featureSize, meanSpacing, minimumCurvature, smoothCurvature);
  {
    int flag = STR_FULLWIDTH;
    int width = 9;
    std::fstream os(outFile, std::fstream::out);
    for (auto& p : pts)
    {
      os << STRstd(p.x, -1, width, flag) << " " << STRstd(p.y, -1, width, flag) << "\n";
    }
  }
  TS_ASSERT_TXT_FILES_EQUAL(baseFile, outFile);
} // MePolyRedistributePtsCurvatureIntermediateTests::testCoastline
#endif