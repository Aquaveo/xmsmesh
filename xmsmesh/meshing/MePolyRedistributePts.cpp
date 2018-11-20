//------------------------------------------------------------------------------
/// \file
/// \ingroup meshing
/// \copyright (C) Copyright Aquaveo 2018. Distributed under FreeBSD License
/// (See accompanying file LICENSE or https://aqaveo.com/bsd/license.txt)
//------------------------------------------------------------------------------

//----- Included files ---------------------------------------------------------

// 1. Precompiled header

// 2. My own header
#include <xmsmesh/meshing/MePolyRedistributePts.h>

// 3. Standard library headers
#include <cfloat>

// 4. External library headers
#include <boost/make_shared.hpp>
#include <boost/thread/mutex.hpp>

// 5. Shared code headers
#include <xmscore/points/pt.h>
#include <xmscore/math/math.h>
#include <xmscore/stl/vector.h>
#include <xmsinterp/geometry/GmMultiPolyIntersector.h>
#include <xmsinterp/geometry/GmMultiPolyIntersectionSorter.h>
#include <xmsinterp/geometry/GmMultiPolyIntersectionSorterTerse.h>
#include <xmsinterp/interpolate/InterpBase.h>
#include <xmsmesh/meshing/detail/MePolyOffsetter.h>
#include <xmsmesh/meshing/detail/MePolyRedistributePtsCurvature.h>
#include <xmscore/misc/xmstype.h>
#include <xmscore/misc/XmError.h>
#include <xmscore/misc/XmLog.h>
#include <xmsinterp/triangulate/TrTriangulatorPoints.h>
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

class MePolyRedistributePtsImpl : public MePolyRedistributePts
{
public:
  MePolyRedistributePtsImpl()
  : m_constSize(XM_NONE)
  , m_minLength(XM_DBL_HIGHEST)
  , m_maxLength(XM_DBL_LOWEST) // changed to be XM_DBL_LOWEST since it is an arbitrarily low value
  , m_sizeBias(1.0)
  , m_polyPts(new VecPt3d())
  , m_intersectWithTris(false)
  , m_distSqTol(0)
  , m_biasConstSize(false)
  , m_polyOffsetIter(0)
  , m_featureSizeCurvature(0)
  , m_meanSpacingCurvature(0)
  , m_minimumCurvature(0.001)
  , m_smoothCurvature(false)
  , m_curvatureRedist()
  {
  }

  virtual void SetSizeFunc(BSHP<InterpBase> a_interp) override;
  virtual void SetSizeFuncFromPoly(const VecPt3d& a_outPoly,
                                   const VecPt3d2d& a_inPolys,
                                   double a_sizeBias) override;

  //------------------------------------------------------------------------------
  /// \brief Sets the size function to a constant value
  /// \param a_size The element edge size.
  //------------------------------------------------------------------------------
  virtual void SetConstantSizeFunc(double a_size) override
  {
    m_constSize = a_size;
    m_curvatureRedist.reset(); // remove curvature redistribution
  }                            // SetConstantSizeFunc

  //------------------------------------------------------------------------------
  /// \brief Sets the bias for constant value size function
  /// \param a_sizeBias Transition rate for size function
  //------------------------------------------------------------------------------
  virtual void SetConstantSizeBias(double a_sizeBias) override
  {
    m_sizeBias = a_sizeBias;
    m_biasConstSize = true;
    m_curvatureRedist.reset(); // remove curvature redistribution
  }                            // SetConstantSizeBias

  void SetUseCurvatureRedistribution(double a_featureSize,
                                     double a_meanSpacing,
                                     double a_minimumCurvature,
                                     bool a_smooth) override;

  void Redistribute(const MePolyOffsetterOutput& a_input,
                    MePolyOffsetterOutput& a_out,
                    int a_polyOffsetIter);
  virtual VecPt3d Redistribute(const VecPt3d& a_polyLine) override;
  virtual double SizeFromLocation(const Pt3d& a_location) override;

  VecPt3d LoopToVecPt3d(const VecSizet& a_idx, const VecPt3d& a_pts);
  void IntersectWithTris(VecPt3d& a_pts);
  void InterpEdgeLengths(const VecPt3d& a_pts, VecDbl& lengths);
  void InterpToPoint(size_t a_idx,
                     const VecPt3d& a_pts,
                     VecDbl& a_lengths,
                     VecDbl& a_wt,
                     VecDbl& a_d2,
                     double& a_bias);
  void CalcInterpWeights(const Pt3d& a_pt, VecDbl& a_wt, VecDbl& a_d2, double& a_sumWt);
  VecPt3d RedistPts(const VecPt3d& a_pts, const VecDbl& lengths);
  VecPt3d RedistPts2(const VecPt3d& a_pts, const VecDbl& lengths);
  void RedistPtsToOutput(const VecPt3d& a_pts, int a_polyType, MePolyOffsetterOutput& a_out);

  void TvaluesForSeg(size_t a_idx,
                     size_t a_nextIdx,
                     const VecDbl& a_segLengths,
                     const VecDbl& a_segTvalues,
                     const VecDbl& a_interpLengths,
                     double& a_pcntin,
                     VecDbl& a_tVals);
  void SizeFromPolyCalcAveEdgeLengthAtPts(const VecPt3d& a_outPoly, const VecPt3d2d& a_inPolys);
  void SizeFromPolyAddEdgeLengths(const VecPt3d& a_pts);
  void CalcSegLengths(const VecPt3d& a_pts, VecDbl& a_segLength, VecDbl& _segTvalues);
  void GetSegmentFromTval(const VecDbl& a_segTvalues,
                          double a_tVal,
                          size_t& a_segIdx,
                          double& a_segT0,
                          double& a_segT1);
  void CreatePolyIntersector();

  double m_constSize,        ///< constant size function
    m_minLength,             ///< min segment length in polygon
    m_maxLength,             ///< max segment length in polygon
    m_sizeBias;              ///< transition factor for size function
  BSHP<VecPt3d> m_polyPts;   ///< polygon point locations
  VecDbl m_polyEdgeLengths;  ///< edge lengths of the polygon
  BSHP<InterpBase> m_interp; ///< interpolation class given to this class
  // BSHP<InterpIdw>              m_idw; ///< interpolation class used with larger numbers of
  // polygon points
  BSHP<VecPt3d> m_sizePts;                        ///< points from the m_interp class
  BSHP<VecInt> m_sizeTris;                        ///< triangles from the m_interp class
  BSHP<GmMultiPolyIntersector> m_polyIntersector; ///< used to intersect the polygon with m_sizeTris
  VecInt2d m_polys; ///< polygon definition of triangles for m_polyIntersector
  bool
    m_intersectWithTris; ///< flag to indicate that polygon should be intersected with the triangles
  double m_distSqTol;    ///< tolerance used to speed up interpolation
  bool m_biasConstSize;  ///< flag to indicate transitioning to constant size function
  int m_polyOffsetIter;  ///< number of iterations from the polygon boundary
  /// Used by curvature redistribution. The size of the smallest feature in the polyline to be
  /// detected. Large values will generate point distributions that follow coarser curvatures.
  double m_featureSizeCurvature;
  /// Used by curvature redistribution. The mean spacing between the distributed points.
  double m_meanSpacingCurvature;
  /// Used by curvature redistribution. The value of the curvature to be used instead of 0 in
  /// staight lines. It limits the maximum spacing between points. If not included, the default
  /// is 0.001.
  double m_minimumCurvature;
  /// Used by curvature redistribution. Detemines if the curvatures are to be averaged by a
  /// rolling 0.25-0.5-0.25 weighted rolling average.
  bool m_smoothCurvature;
  /// Point redistributor that uses curvature
  BSHP<MePolyRedistributePtsCurvature> m_curvatureRedist;
};

//------------------------------------------------------------------------------
/// \brief Creates an instance of this class
/// \return MePolyRedistributePts
//------------------------------------------------------------------------------
BSHP<MePolyRedistributePts> MePolyRedistributePts::New()
{
  BSHP<MePolyRedistributePts> ret(new MePolyRedistributePtsImpl);
  return ret;
} // MePolyRedistributePts::New
//------------------------------------------------------------------------------
/// \brief destructor
//------------------------------------------------------------------------------
MePolyRedistributePts::~MePolyRedistributePts()
{
} // MePolyRedistributePts::~MePolyRedistributePts
//------------------------------------------------------------------------------
/// \brief constructor
//------------------------------------------------------------------------------
MePolyRedistributePts::MePolyRedistributePts()
{
} // MePolyRedistributePts::MePolyRedistributePts

////////////////////////////////////////////////////////////////////////////////
/// \class MePolyRedistributePtsImpl
/// \brief Redistributes points along polylines using a size function.
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief Sets the size function interpolator
/// \param a_interp Size function interpolator class
//------------------------------------------------------------------------------
void MePolyRedistributePtsImpl::SetSizeFunc(BSHP<InterpBase> a_interp)
{
  m_curvatureRedist.reset(); // remove curvature redistribution
  m_interp = a_interp;
  m_intersectWithTris = true;
  // make sure we have triangles
  m_sizePts = m_interp->GetPts();
  m_sizeTris = m_interp->GetTris();
  if (m_intersectWithTris && m_sizeTris->empty())
  {
    VecInt vTris;
    TrTriangulatorPoints tri(*m_sizePts, vTris);
    tri.Triangulate();
    if (!vTris.empty())
    {
      m_sizeTris = BSHP<VecInt>(new VecInt());
      m_sizeTris->reserve(vTris.size());
      for (size_t i = 0; i < vTris.size(); ++i)
        m_sizeTris->push_back((int)vTris[i]);
    }
  }
  if (m_intersectWithTris)
    CreatePolyIntersector();
} // SetSizeFunc
//------------------------------------------------------------------------------
/// \brief Creates an interpolator that uses the spacing on the input polygon
/// as its scalar
/// \param a_outPoly The outside polygon
/// \param a_inPolys Inside polygons that are inside of a_outPoly
/// \param a_sizeBias A factor used in transitioning the size
//------------------------------------------------------------------------------
void MePolyRedistributePtsImpl::SetSizeFuncFromPoly(const VecPt3d& a_outPoly,
                                                    const VecPt3d2d& a_inPolys,
                                                    double a_sizeBias)
{
  m_curvatureRedist.reset(); // remove curvature redistribution
  m_sizeBias = a_sizeBias;
  // calculate the lengths of the segments
  SizeFromPolyCalcAveEdgeLengthAtPts(a_outPoly, a_inPolys);
  if (m_polyPts->empty())
    return;
  double diff = m_maxLength - m_minLength;
  if (diff / m_minLength < .05)
  {
    m_constSize = (m_maxLength + m_minLength) / 2;
  }
  // looking for a way to speed up interpolation of edge lengths
  // doesn't work well
  // else if (m_polyPts->size() > 100)
  //{
  //  m_idw = InterpIdw::New();
  //  m_idw->SetPtsTris(m_polyPts, m_sizeTris);
  //  BSHP<VecFlt> s(new VecFlt(m_polyEdgeLengths.size()));
  //  for (size_t i=0; i<m_polyEdgeLengths.size(); ++i)
  //    (*s)[i] = (float)m_polyEdgeLengths[i];
  //  m_idw->SetScalars(s);
  //  m_idw->SetSearchOpts(32, false);
  //  m_idw->SetWeightCalcMethod(InterpIdw::CLASSIC);
  //}
  else
  {
    size_t nPts = m_polyPts->size();
    Pt3d* p(&(*m_polyPts)[0]);
    Pt3d pMin(p[0]), pMax(p[0]);
    for (size_t i = 0; i < nPts; ++i)
    {
      if (p[i].x < pMin.x)
        pMin.x = p[i].x;
      if (p[i].y < pMin.y)
        pMin.y = p[i].y;
      if (p[i].x > pMax.x)
        pMax.x = p[i].x;
      if (p[i].y > pMax.y)
        pMax.y = p[i].y;
    }
    // compute distance to avoid interpolating too much
    double boundsDistSq = MdistSq(pMin.x, pMin.y, pMax.x, pMax.y);
    m_distSqTol = 1e-4 * boundsDistSq;
  }
} // MePolyRedistributePtsImpl::SetSizeFuncFromPoly
//------------------------------------------------------------------------------
/// \brief Specifies that curvature redistribution will be used
/// \param[in] a_featureSize: The size of the smallest feature in the polyline to be detected.
///   Large values will generate point distributions that follow coarser curvatures.
/// \param[in] a_meanSpacing: The mean spacing between the distributed points.
/// \param[in] a_minimumCurvature: The value of the curvature to be used instead of 0 in staight
///   lines. It limits the maximum spacing between points. If not included, the default is 0.001.
/// \param[in] a_smooth: Detemines if the curvatures are to be averaged by a rolling 0.25-0.5-0.25
///   weighted rolling average.
//------------------------------------------------------------------------------
void MePolyRedistributePtsImpl::SetUseCurvatureRedistribution(double a_featureSize,
                                                              double a_meanSpacing,
                                                              double a_minimumCurvature,
                                                              bool a_smooth)
{
  m_featureSizeCurvature = a_featureSize;
  m_meanSpacingCurvature = a_meanSpacing;
  m_minimumCurvature = a_minimumCurvature;
  m_smoothCurvature = a_smooth;
  m_curvatureRedist = MePolyRedistributePtsCurvature::New();
} // MePolyRedistributePtsImpl::SetUseCurvatureRedistribution
//------------------------------------------------------------------------------
/// \brief Redistributes points on closed loop polylines. The length of edges
/// in the redistribution comes from a size function that is interpolated to
/// the points that make up the polylines. By default this size function comes
/// from the edge lengths in the original polygon.
/// \param a_input Input closed loop polylines
/// \param a_out Redistributed closed loop polylines
/// \param a_polyOffsetIter Number of iterations from the polygon boundary
//------------------------------------------------------------------------------
void MePolyRedistributePtsImpl::Redistribute(const MePolyOffsetterOutput& a_input,
                                             MePolyOffsetterOutput& a_out,
                                             int a_polyOffsetIter)
{
  VecDbl lengths;
  m_polyOffsetIter = a_polyOffsetIter;
  a_out.m_loops.resize(0);
  a_out.m_pts.resize(0);
  a_out.m_loopTypes.resize(0);
  for (size_t i = 0; i < a_input.m_loops.size(); ++i)
  {
    const VecSizet& loop(a_input.m_loops[i]);
    const int& lType(a_input.m_loopTypes[i]);
    VecPt3d pts = LoopToVecPt3d(loop, a_input.m_pts);
    if (m_intersectWithTris)
    {
      IntersectWithTris(pts);
    }
    pts.push_back(pts.front());
    // interpolate edge lengths
    InterpEdgeLengths(pts, lengths);
    // redistribute the points
    VecPt3d redistPts = RedistPts(pts, lengths);
    if (!redistPts.empty())
      redistPts.pop_back();
    RedistPtsToOutput(redistPts, lType, a_out);
  }
  m_polyOffsetIter = 0;
} // MePolyPaverToMeshPtsImpl::Redistribute
//------------------------------------------------------------------------------
/// \brief Redistributes points on a polyline.
/// \param[in] a_polyLine Input polyline
/// \return Redistributed polyline
//------------------------------------------------------------------------------
VecPt3d MePolyRedistributePtsImpl::Redistribute(const VecPt3d& a_polyLine)
{
  if (m_curvatureRedist)
  {
    return m_curvatureRedist->Redistribute(a_polyLine, m_featureSizeCurvature,
                                           m_meanSpacingCurvature, m_minimumCurvature,
                                           m_smoothCurvature);
  }
  VecPt3d pts(a_polyLine), ret;
  VecDbl lengths;
  if (m_intersectWithTris)
  {
    IntersectWithTris(pts);
  }
  // interpolate edge lengths
  InterpEdgeLengths(pts, lengths);
  // redistribute the points
  ret = RedistPts(pts, lengths);
  return ret;
} // MePolyRedistributePtsImpl::Redistribute
//------------------------------------------------------------------------------
/// \brief returns a size based on the xy location
/// \param[in] a_location: the location
/// \return size function value for the desired edge length at the location
//------------------------------------------------------------------------------
double MePolyRedistributePtsImpl::SizeFromLocation(const Pt3d& a_location)
{
  if (m_curvatureRedist)
  {
    std::string msg = "MePolyRedistributePts set to use curvature redistribution; "
      "MePolyRedistributePtsImpl::SizeFromLocation can not be call with these "
      "settings. XM_NODATA will be returned.";
    boost::mutex mtx;
    mtx.lock();
    XM_LOG(xmlog::error, msg);
    mtx.unlock();
    return XM_NODATA;
  }
  VecPt3d pts(1, a_location);
  VecDbl lengths;
  InterpEdgeLengths(pts, lengths);
  return lengths.front();
} // MePolyRedistributePtsImpl::SizeFromLocation
//------------------------------------------------------------------------------
/// \brief Creates a vector of pts from indices into another vector of points
/// \param a_idx Vector of indices into the a_pts vector. Defines the closed
/// loop polyline.
/// \param a_pts Vector of locations.
/// \return Vector of locations defining the loop.
//------------------------------------------------------------------------------
VecPt3d MePolyRedistributePtsImpl::LoopToVecPt3d(const VecSizet& a_idx, const VecPt3d& a_pts)
{
  VecPt3d ret;
  ret.reserve(a_idx.size());
  for (size_t i = 0; i < a_idx.size(); ++i)
    ret.push_back(a_pts[a_idx[i]]);
  return ret;
} // MePolyRedistributePtsImpl::LoopToVecPt3d
//------------------------------------------------------------------------------
/// \brief Creates a vector of pts from indices into another vector of points
/// \param a_pts Vector of locations.
//------------------------------------------------------------------------------
void MePolyRedistributePtsImpl::IntersectWithTris(VecPt3d& a_pts)
{
  VecPt3d newPts, pts;
  VecInt polys;
  for (size_t i = 0; i < a_pts.size(); ++i)
  {
    Pt3d p0(a_pts[i]), p1(a_pts[0]);
    if (i < a_pts.size() - 1)
      p1 = a_pts[i + 1];
    m_polyIntersector->TraverseLineSegment(p0.x, p0.y, p1.x, p1.y, polys, pts);
    auto start = pts.begin();
    if (i > 0)
      ++start;
    newPts.insert(newPts.end(), start, pts.end());
  }
  a_pts.swap(newPts);
  a_pts.pop_back();
} // MePolyRedistributePtsImpl::IntersectWithTris
//------------------------------------------------------------------------------
/// \brief Interpolates edge lengths to each of the points that make up a
/// polyline
/// \param a_pts Vector of locations.
/// \param a_lengths Vector of interpolated lengths filled in by the method.
/// Will be the same size as a_idx.
//------------------------------------------------------------------------------
void MePolyRedistributePtsImpl::InterpEdgeLengths(const VecPt3d& a_pts, VecDbl& a_lengths)
{
  a_lengths.assign(a_pts.size(), 0.0);
  if (XM_NONE != m_constSize && !m_biasConstSize)
  {
    a_lengths.assign(a_pts.size(), m_constSize);
    return;
  }
  double bias(XM_NONE);
  if (XM_NONE != m_constSize && m_biasConstSize)
  {
    bias = m_sizeBias;
    if (bias > .99)
      bias = .99;
    else if (bias < .01)
      bias = .01;
    bias = 1 - bias;
    bias = pow(bias, (double)m_polyOffsetIter);
  }

  if (m_interp)
  { // size function interpolation
    VecFlt s;
    m_interp->InterpToPts(a_pts, s);
    for (size_t i = 0; i < s.size(); ++i)
      a_lengths[i] = (double)s[i];
  }
  // tried to use to speed up but was not successful
  // else if (m_idw)
  //{
  //  VecFlt s;
  //  m_idw->InterpToPts(a_pts, s);
  //  for (size_t i=0; i<s.size(); ++i) a_lengths[i] = (double)s[i];
  //}
  else if (!m_polyEdgeLengths.empty())
  {
    double lastInterp(-1), distSq(0);
    VecDbl wt, d2;
    Pt3d aPt;
    InterpToPoint(0, a_pts, a_lengths, wt, d2, bias);
    lastInterp = a_lengths[0];
    aPt = a_pts[0];

    // do the rest of the points
    for (size_t i = 1; i < a_pts.size(); ++i)
    {
      distSq = MdistSq(aPt.x, aPt.y, a_pts[i].x, a_pts[i].y);
      if (distSq < m_distSqTol)
      {
        a_lengths[i] = lastInterp;
        continue;
      }
      InterpToPoint(i, a_pts, a_lengths, wt, d2, bias);
      lastInterp = a_lengths[i];
      aPt = a_pts[i];
    }
  }
  else
  {
    boost::mutex mtx;
    mtx.lock();
    VecDbl lengths, tvals;
    CalcSegLengths(a_pts, lengths, tvals);
    double sum(0);
    for (size_t i = 0; i < lengths.size(); ++i)
      sum += lengths[i];
    sum = sum / lengths.size();
    m_constSize = sum;
    std::stringstream ss;
    ss << "Interpolator not defined in MePolyRedistributePts. Size function "
          "set to constant value: "
       << sum << ".";
    XM_LOG(xmlog::debug, ss.str());
    InterpEdgeLengths(a_pts, a_lengths);
    mtx.unlock();
  }
} // MePolyRedistributePtsImpl::InterpEdgeLengths
//------------------------------------------------------------------------------
/// \brief Interpolates to a point from the boundary of the polygon
/// \param[in] a_idx Index of point
/// \param[in] a_pts Vector of point locations
/// \param[out] a_lengths Vector of lengths. Value at a_idx computed.
/// \param[out] a_wt Vector of weights from polygon points filled in
/// CalcInterpWeights.
/// \param[out] a_d2 Vector of distance squared between polygon points and the
/// point at a_idx that is filled by CalcInterpWeights.
/// \param[in] a_bias A factor for scaling the length when doing a smooth
/// transition to a constant size function.
//------------------------------------------------------------------------------
void MePolyRedistributePtsImpl::InterpToPoint(size_t a_idx,
                                              const VecPt3d& a_pts,
                                              VecDbl& a_lengths,
                                              VecDbl& a_wt,
                                              VecDbl& a_d2,
                                              double& a_bias)
{
  double sumWt(0);
  CalcInterpWeights(a_pts[a_idx], a_wt, a_d2, sumWt);
  for (size_t j = 0; j < m_polyEdgeLengths.size(); ++j)
  {
    a_lengths[a_idx] += m_polyEdgeLengths[j] * (a_wt[j] / sumWt);
  }
  if (XM_NONE != a_bias)
  {
    double& d(a_lengths[a_idx]);
    if (d > m_constSize)
    {
      d = d * a_bias;
      if (d < m_constSize)
        d = m_constSize;
    }
    else
    {
      d = d / a_bias;
      if (d > m_constSize)
        d = m_constSize;
    }
  }
  else
  {
    if (a_lengths[a_idx] < m_minLength)
      a_lengths[a_idx] = m_minLength;
    if (a_lengths[a_idx] > m_maxLength)
      a_lengths[a_idx] = m_maxLength;
  }
} // MePolyRedistributePtsImpl::InterpToPoint
//------------------------------------------------------------------------------
/// \brief Calculates the interpolation weights of the points used to do the
/// interpolation
/// \param a_pt The location being interpolated to.
/// \param a_wt Vector of weights calculated by this function
/// \param a_d2 Vector of distance squared
/// \param a_sumWt Sum of the weights in a_wt
//------------------------------------------------------------------------------
void MePolyRedistributePtsImpl::CalcInterpWeights(const Pt3d& a_pt,
                                                  VecDbl& a_wt,
                                                  VecDbl& a_d2,
                                                  double& a_sumWt)
{
  a_sumWt = 0;
  double diff(0), factor(0);
  size_t nPts = m_polyPts->size();
  const Pt3d* p(&(*m_polyPts)[0]);
  a_wt.resize(nPts);
  a_d2.resize(nPts);
  for (size_t i = 0; i < nPts; ++i)
  {
    a_d2[i] = MdistSq(a_pt.x, a_pt.y, p[i].x, p[i].y);
    if (a_d2[i] < 10e-8)
      a_wt[i] = 10e10;
    else
      a_wt[i] = 1 / a_d2[i];

    diff = m_polyEdgeLengths[i] - m_minLength;
    factor = m_minLength + (m_sizeBias * diff);
    a_wt[i] *= factor;

    a_sumWt += a_wt[i];
  }
} // MePolyRedistributePtsImpl::CalcInterpWeights
//------------------------------------------------------------------------------
/// \brief Uses interpolated lengths to redistribute points on a polyline
/// \param a_pts Vector of locations.
/// \param a_interpLengths Vector of interpolated lengths.
/// \return Vector of locations of redistributed points
//------------------------------------------------------------------------------
VecPt3d MePolyRedistributePtsImpl::RedistPts(const VecPt3d& a_pts, const VecDbl& a_interpLengths)
{
  static bool debug(false);
  if (debug)
    return RedistPts2(a_pts, a_interpLengths);

  // Get lengths of all segments and the total length
  VecDbl segLength, segTvalues;
  CalcSegLengths(a_pts, segLength, segTvalues);

  // create array of T values for locations of new points along this loop
  double pcntin(0);
  size_t nextIdx;
  VecDbl tVals;
  tVals.reserve(a_pts.size());
  tVals.push_back(0.0); // this is the first and last point in our loop
  for (size_t i = 0; i < a_pts.size() - 1; ++i)
  {
    nextIdx = 0;
    if (i + 1 < a_pts.size() - 1)
      nextIdx = i + 1;
    TvaluesForSeg(i, nextIdx, segLength, segTvalues, a_interpLengths, pcntin, tVals);
  }
  double aveTincrement(0);
  for (size_t i = 1; i < tVals.size(); ++i)
  {
    aveTincrement += tVals[i] - tVals[i - 1];
  }
  aveTincrement /= (tVals.size() - 1);

  VecPt3d ret;
  // adjust the tvals based on what was left over when processing the last
  // segment. If the left over tval is less than .5 of the last segment then
  // we don't keep that last tval
  double leftOverTval = 1 - tVals.back();
  if (tVals.size() < 2)
    return ret;
  if (pcntin < .5)
  {
    tVals.pop_back();
    if (leftOverTval >= 0)
      leftOverTval = 1 - tVals.back();
  }
  if (tVals.size() > 2)
  { // redistribute the leftover tvalue to all the tVals
    double redistT = (leftOverTval - aveTincrement) / tVals.size();
    for (size_t i = 1; i < tVals.size(); ++i)
      tVals[i] += (i * redistT);
  }
  else if (segLength.size() > 3)
  { // we initially had at least 3 points
    tVals.assign(3, 0);
    tVals[1] = 1 / 3.0;
    tVals[2] = 2 / 3.0;
  }

  if (tVals.size() > 2)
  {
    tVals.push_back(1.0);
    ret.reserve(tVals.size());
    // create the points
    // put in the first point
    ret.push_back(a_pts[0]);
    size_t segIdx(0);
    for (size_t i = 1; i < tVals.size(); ++i)
    { // get the segment that the tval is on
      double segT0, segT1, tVal = tVals[i];
      GetSegmentFromTval(segTvalues, tVals[i], segIdx, segT0, segT1);
      // convert polyline tval to segment tval
      tVal = (tVal - segT0) / (segT1 - segT0);
      // calculate the point
      size_t nextIdx = 0;
      if (segIdx + 1 < a_pts.size())
        nextIdx = segIdx + 1;
      const Pt3d &p0(a_pts[segIdx]), &p1(a_pts[nextIdx]);
      Pt3d pt;
      pt.x = p0.x + tVal * (p1.x - p0.x);
      pt.y = p0.y + tVal * (p1.y - p0.y);
      ret.push_back(pt);
    }
  }
  return ret;
} // MePolyRedistributePtsImpl::RedistPts
//------------------------------------------------------------------------------
/// \brief Uses interpolated lengths to redistribute points on a polyline
/// \param a_pts Vector of locations.
/// \param a_interpLengths Vector of interpolated lengths.
/// \return Vector of locations of redistributed points
//------------------------------------------------------------------------------
VecPt3d MePolyRedistributePtsImpl::RedistPts2(const VecPt3d& a_pts, const VecDbl& a_interpLengths)
{
  VecPt4d finalpts;
  finalpts.reserve(a_pts.size());
  int numint((int)a_pts.size());
  double pcntin(0), pcntleft(0);
  int nextid(1);
  Pt3d prevpt(a_pts[0]);
  prevpt.z = a_interpLengths[0];
  Pt3d nextpt(a_pts[1]);
  nextpt.z = a_interpLengths[1];
  double len_segs(0), len_arc(0), length(0);
  len_arc = length = Mdist(prevpt.x, prevpt.y, nextpt.x, nextpt.y);
  int numfinal(0);
  double totdist(0), w1(0), w2(0);
  bool done(false);
  while (!done)
  {
    w1 = prevpt.z; // weight at first point
    w2 = nextpt.z; // weight at second point
    pcntleft = pcntin + 2 * length / (w1 + w2);
    // check if we should add a point
    if (pcntleft >= 1)
    { // calc new pt location
      double scale = 1.0 - pcntin;
      double ptlength = 2 * scale * w1 * length / (2 * length + scale * (w1 - w2));
      totdist += ptlength;
      // weight at new point
      double pcnt = ptlength / length;
      double newweight = w1 * (1 - pcnt) + w2 * pcnt;
      // add to final points
      ++numfinal;
      finalpts.push_back(Pt4d());
      finalpts[numfinal - 1].x = prevpt.x + pcnt * (nextpt.x - prevpt.x);
      finalpts[numfinal - 1].y = prevpt.y + pcnt * (nextpt.y - prevpt.y);
      finalpts[numfinal - 1].z = newweight;
      finalpts[numfinal - 1].w = totdist;
      len_segs += totdist;
      prevpt = finalpts[numfinal - 1];
      totdist = 0.0;
      // get length of segment left after inserting pt
      length -= ptlength;
      pcntin = 0.0;
    }
    else
    { // go on to next segment
      pcntin = pcntleft;
      if (++nextid < numint)
      {
        prevpt = nextpt;

        nextpt = a_pts[nextid];
        nextpt.z = a_interpLengths[nextid];
        totdist += length;
        length = Mdist(prevpt.x, prevpt.y, nextpt.x, nextpt.y);
        len_arc += length;
      }
      else
        done = true;
    }
  }

  // make sure we still have an arc
  if (numfinal > 0)
  {
    double len_segs_to_adjust = 0.0;
    // init length of the arc to adjust to remaining portion
    double len_arc_to_adjust = len_arc - len_segs;
    // init first point to move as the last in the list
    int first_pt_to_move = numfinal - 1;
    // eliminate last part_seg if less than half
    if (pcntin < 0.5)
      --numfinal;
    // keep the length of the last seg if more than half
    else
      len_segs_to_adjust = (length + totdist) / pcntleft;
    // find segs with length close to final length to find first pt to move
    while (first_pt_to_move >= 0 && finalpts[first_pt_to_move].w < w1 * 2.0 &&
           finalpts[first_pt_to_move].w > w1 * 0.5)
    {
      len_segs_to_adjust += finalpts[first_pt_to_move].w;
      len_arc_to_adjust += finalpts[first_pt_to_move].w;
      --first_pt_to_move;
    }
    ++first_pt_to_move;
    // find error adjustment scale
    double scale = len_arc_to_adjust / len_segs_to_adjust;
    /* reset the locations using scaled distance */
    nextid = 1;
    prevpt = a_pts[0];
    prevpt.z = a_interpLengths[0];
    nextpt = a_pts[1];
    nextpt.z = a_interpLengths[1];
    length = Mdist(prevpt.x, prevpt.y, nextpt.x, nextpt.y);
    for (int i = 0; i < numfinal; i++)
    {
      if (i < first_pt_to_move)
      {
        totdist = finalpts[i].w;
        // move along int points until reach finalpt
        while (totdist > length)
        {
          totdist -= length;
          // go onto next segment
          prevpt = nextpt;
          ++nextid;
          nextpt = a_pts[nextid];
          nextpt.z = a_interpLengths[nextid];
          length = Mdist(prevpt.x, prevpt.y, nextpt.x, nextpt.y);
        }
      }
      else
      {
        totdist = finalpts[i].w * scale;
        // move along int points until reach finalpt
        while (totdist > length)
        {
          totdist -= length;
          // go onto next segment
          prevpt = nextpt;
          ++nextid;
          nextpt = a_pts[nextid];
          nextpt.z = a_interpLengths[nextid];
          length = Mdist(prevpt.x, prevpt.y, nextpt.x, nextpt.y);
        }
        // reached point to move
        double pcnt = totdist / length;
        finalpts[i].x = prevpt.x + pcnt * (nextpt.x - prevpt.x);
        finalpts[i].y = prevpt.y + pcnt * (nextpt.y - prevpt.y);
      }
      prevpt = finalpts[i];
      length -= totdist;
    }
  }

  VecPt3d ret;
  // if arc closes, make sure there are at least 3 points on arc
  if (a_pts.front() == a_pts.back() && numfinal < 2)
  {
    if (a_pts.size() < 3)
      return ret;

    finalpts.resize(2);
    finalpts[0] = a_pts[a_pts.size() / 3];
    finalpts[1] = a_pts[a_pts.size() * 2 / 3];
  }
  ret.assign(finalpts.size() + 2, Pt3d());
  ret.front() = a_pts.front();
  ret.back() = a_pts.back();
  for (size_t i = 1; i < finalpts.size(); ++i)
    ret[i] = finalpts[i];
  return ret;
} // MePolyRedistributePtsImpl::RedistPts2
//------------------------------------------------------------------------------
/// \brief Move the redistributed points to the output class
/// \param a_pts Vector of locations of the redistributed points.
/// \param a_polyType The type of polygon that is represented by a_idx
/// (MePolyOffsetter::polytype)
/// \param a_out The output filled by this method. Creates a new polygon with
/// redistributed points.
//------------------------------------------------------------------------------
void MePolyRedistributePtsImpl::RedistPtsToOutput(const VecPt3d& a_pts,
                                                  int a_polyType,
                                                  MePolyOffsetterOutput& a_out)
{
  if (a_pts.empty())
    return;
  a_out.m_loopTypes.push_back(a_polyType);
  a_out.m_loops.push_back(VecSizet());
  for (size_t i = 0; i < a_pts.size(); ++i)
  {
    a_out.m_loops.back().push_back(a_out.m_pts.size());
    a_out.m_pts.push_back(a_pts[i]);
  }
} // MePolyRedistributePtsImpl::RedistPtsToOutput
//------------------------------------------------------------------------------
/// \brief A single segment of a polygon is analyzed to see where points are
/// redistributed on that segment.
/// \param a_idx Index of the segment
/// \param a_nextIdx The next index in the array. Since this is a closed loop
/// this value can wrap back to 0
/// \param a_segLengths Vector of segment lengths
/// \param a_segTvalues Vector of parametric values for the end points of
/// segments relative to the total polyline length
/// \param a_interpLengths Vector of interpolated edgelengths at the end points
/// of the segments
/// \param a_pcntin A running percentage that is added to by this method.
/// \param a_tVals Vector of parametric values for redistributed points. Filled
/// by this method.
//------------------------------------------------------------------------------
void MePolyRedistributePtsImpl::TvaluesForSeg(size_t a_idx,
                                              size_t a_nextIdx,
                                              const VecDbl& a_segLengths,
                                              const VecDbl& a_segTvalues,
                                              const VecDbl& a_interpLengths,
                                              double& a_pcntin,
                                              VecDbl& a_tVals)
{
  // current segment length
  double length(a_segLengths[a_idx]);
  // tvalues at the begin and end of current segment
  double t1(a_segTvalues[a_idx]), t2(a_segTvalues[a_nextIdx]);
  if (0 == a_nextIdx)
    t2 = 1;
  // target segment length at begin and end of segment
  double w1(a_interpLengths[a_idx]), w2(a_interpLengths[a_nextIdx]);

  bool done(false);
  while (!done)
  {
    double pcntleft = a_pcntin + (2 * length) / (w1 + w2);
    if (pcntleft >= 1.0)
    {
      double scale = 1.0 - a_pcntin;
      double ptlength = 2 * scale * w1 * length / (2 * length + scale * (w1 - w2));
      double pcnt = ptlength / length;
      double t = (pcnt * (t2 - t1)) + t1; // convert to polyline tval
      a_tVals.push_back(t);

      w1 = w1 * (1 - pcnt) + w2 * pcnt;
      t1 = t;
      length -= ptlength;
      a_pcntin = 0;
    }
    else
    {
      a_pcntin = pcntleft;
      done = true;
    }
  }
} // MePolyRedistributePtsImpl::TvaluesForSeg
//------------------------------------------------------------------------------
/// \brief Calculates average edge lengths at each point in a loop defining
/// a polygon.
/// \param a_outPoly Outside polygon
/// \param a_inPolys Inside polygons that are in a_outPoly
//------------------------------------------------------------------------------
void MePolyRedistributePtsImpl::SizeFromPolyCalcAveEdgeLengthAtPts(const VecPt3d& a_outPoly,
                                                                   const VecPt3d2d& a_inPolys)
{
  m_minLength = XM_DBL_HIGHEST;
  m_maxLength = XM_DBL_LOWEST;
  SizeFromPolyAddEdgeLengths(a_outPoly);
  for (size_t i = 0; i < a_inPolys.size(); ++i)
    SizeFromPolyAddEdgeLengths(a_inPolys[i]);
} // MePolyRedistributePtsImpl::SizeFromPolyCalcAveEdgeLengthAtPts
//------------------------------------------------------------------------------
/// \brief Calculates average edge lengths at each point in a loop defining
/// a polygon.
/// \param a_pts Locations that define a polygon
//------------------------------------------------------------------------------
void MePolyRedistributePtsImpl::SizeFromPolyAddEdgeLengths(const VecPt3d& a_pts)
{
  double l(0), dx, dy;
  size_t start = m_polyEdgeLengths.size();
  m_polyEdgeLengths.insert(m_polyEdgeLengths.end(), a_pts.size(), 0);
  for (size_t i = 0; i < a_pts.size(); ++i)
  {
    size_t nextIdx(0);
    Pt3d nextPt(a_pts[0]);
    if (i + 1 < a_pts.size())
    {
      nextIdx = i + 1;
      nextPt = a_pts[i + 1];
    }
    m_polyPts->push_back(a_pts[i]);

    dx = a_pts[i].x - nextPt.x;
    dy = a_pts[i].y - nextPt.y;
    l = sqrt(dx * dx + dy * dy);
    if (l < m_minLength)
      m_minLength = l;
    if (l > m_maxLength)
      m_maxLength = l;
    l = l / 2;
    m_polyEdgeLengths[i + start] += l;
    m_polyEdgeLengths[nextIdx + start] += l;
  }
} // MePolyRedistributePtsImpl::SizeFromPolyAddEdgeLengths
//------------------------------------------------------------------------------
/// \brief Calculates the lengths of segments and parametric values of segment
/// endpoints.
/// \param a_pts Vector of locations.
/// \param a_segLength Vector of segment lengths. Filled by this method.
/// \param a_segTvalues Vector of parametric values of segment endpoints
/// relative to the perimeter of the closed loop. Filled by this method.
//------------------------------------------------------------------------------
void MePolyRedistributePtsImpl::CalcSegLengths(const VecPt3d& a_pts,
                                               VecDbl& a_segLength,
                                               VecDbl& a_segTvalues)
{
  a_segLength.assign(a_pts.size() - 1, 0);
  a_segTvalues.assign(a_pts.size() - 1, 0);
  double totalLength(0);
  for (size_t i = 1; i < a_pts.size(); ++i)
  {
    const Pt3d &p0(a_pts[i - 1]), &p1(a_pts[i]);
    a_segLength[i - 1] = Mdist(p0.x, p0.y, p1.x, p1.y);
    totalLength += a_segLength[i - 1];
  }
  double sumLen(0);
  for (size_t i = 1; i < a_pts.size() - 1; ++i)
  {
    sumLen += a_segLength[i - 1];
    a_segTvalues[i] = sumLen / totalLength;
  }
} // MePolyRedistributePtsImpl::CalcSegLengths
//------------------------------------------------------------------------------
/// \brief Finds the segment that contains a_tVal
/// \param a_segTvalues Vector of parametric values of segment endpoints
/// relative to the perimeter of the closed loop.
/// \param a_tVal Parametric value along closed loop.
/// \param a_segIdx Index for starting the search of a_segTvalues and can be
/// modified by this method when the segment is found.
/// \param a_segT0 The parametric t value for the start of the found segment
/// \param a_segT1 The parametric t value for the end of the found segment
//------------------------------------------------------------------------------
void MePolyRedistributePtsImpl::GetSegmentFromTval(const VecDbl& a_segTvalues,
                                                   double a_tVal,
                                                   size_t& a_segIdx,
                                                   double& a_segT0,
                                                   double& a_segT1)
{
  if (a_tVal > 1)
  {
    XM_ASSERT(0);
    a_tVal = 1;
  }
  if (a_segIdx >= a_segTvalues.size())
    a_segIdx = 0;

  a_segT0 = a_segTvalues[a_segIdx];
  a_segT1 = 1;
  if (a_segIdx + 1 < a_segTvalues.size())
    a_segT1 = a_segTvalues[a_segIdx + 1];

  bool done = false;
  while (!done)
  {
    if (a_tVal >= a_segT0 && a_tVal <= a_segT1)
      done = true;
    else
    {
      a_segIdx++;
      a_segT0 = a_segTvalues[a_segIdx];
      a_segT1 = 1;
      if (a_segIdx + 1 < a_segTvalues.size())
        a_segT1 = a_segTvalues[a_segIdx + 1];
    }
  }
} // MePolyRedistributePtsImpl::GetSegmentFromTval
//------------------------------------------------------------------------------
/// \brief Creates the poly intersector class member variable
//------------------------------------------------------------------------------
void MePolyRedistributePtsImpl::CreatePolyIntersector()
{
  m_polys.assign(m_sizeTris->size() / 3, VecInt(3, 0));
  int* iPtr(&(*m_sizeTris)[0]);
  size_t idx;
  for (size_t i = 0; i < m_polys.size(); ++i)
  {
    idx = i * 3;
    m_polys[i][0] = iPtr[idx + 0];
    m_polys[i][1] = iPtr[idx + 1];
    m_polys[i][2] = iPtr[idx + 2];
  }
  BSHP<GmMultiPolyIntersectionSorterTerse> sorterTerse =
    boost::make_shared<GmMultiPolyIntersectionSorterTerse>();
  BSHP<GmMultiPolyIntersectionSorter> sorter = BDPC<GmMultiPolyIntersectionSorter>(sorterTerse);
  m_polyIntersector = GmMultiPolyIntersector::New(*m_sizePts, m_polys, sorter);

} // MePolyRedistributePtsImpl::CreatePolyIntersector
//------------------------------------------------------------------------------
/// \brief Free function access an implement method that needs to be hidden
/// from the public interface
/// \param[in] a_redist: a MePolyRedistributePts class
/// \param a_input Input closed loop polylines
/// \param a_out Redistributed closed loop polylines
/// \param a_polyOffsetIter Number of iterations from the polygon boundary
//------------------------------------------------------------------------------
void mePolyPaverRedistribute(BSHP<MePolyRedistributePts> a_redist,
                             const MePolyOffsetterOutput& a_input,
                             MePolyOffsetterOutput& a_out,
                             int a_polyOffsetIter)
{
  BSHP<MePolyRedistributePtsImpl> r = BDPC<MePolyRedistributePtsImpl>(a_redist);
  XM_ENSURE_TRUE(r);
  r->Redistribute(a_input, a_out, a_polyOffsetIter);
} // mePolyPaverRedistribute

} // namespace xms

#ifdef CXX_TEST
////////////////////////////////////////////////////////////////////////////////

#include <xmsmesh/meshing/MePolyRedistributePts.t.h>

#include <xmscore/testing/TestTools.h>
#include <xmsinterp/geometry/geoms.h>

// namespace xms {
using namespace xms;

////////////////////////////////////////////////////////////////////////////////
/// \class MePolyRedistributePtsUnitTests
/// \brief tester for the MePolyRedistributePts class
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief tests creating the class
//------------------------------------------------------------------------------
void MePolyRedistributePtsUnitTests::testCreateClass()
{
  BSHP<MePolyRedistributePts> b = MePolyRedistributePts::New();
  TS_ASSERT(b);
} // MePolyRedistributePtsUnitTests::testCreateClass
//------------------------------------------------------------------------------
/// \brief tests interpolating the edge length
//------------------------------------------------------------------------------
void MePolyRedistributePtsUnitTests::testInterpEdgeLengths()
{
  MePolyRedistributePtsImpl r;
  // input polygon
  VecPt3d outPoly = {{0, 0, 0}, {0, 10, 0}, {10, 10, 0}, {10, 0, 0}};
  VecPt3d2d inPolys;
  r.SetSizeFuncFromPoly(outPoly, inPolys, 1);
  VecPt3d pts = {{1, 1, 0}, {1, 9, 0}, {9, 9, 0}, {9, 1, 0}};
  VecDbl lengths;
  r.InterpEdgeLengths(pts, lengths);
  VecDbl baseLengths = {10.0, 10.0, 10.0, 10.0};
  TS_ASSERT_EQUALS_VEC(baseLengths, lengths);
} // MePolyRedistributePtsUnitTests::testInterpEdgeLengths
//------------------------------------------------------------------------------
/// \brief tests interpolating the edge length
//------------------------------------------------------------------------------
void MePolyRedistributePtsUnitTests::testInterpEdgeLengths1()
{
  MePolyRedistributePtsImpl r;
  // input polygon
  VecPt3d outPoly = {{0, 0, 0}, {0, 2, 0},  {0, 4, 0},   {0, 6, 0},
                     {0, 8, 0}, {0, 10, 0}, {10, 10, 0}, {10, 0, 0}};
  VecPt3d2d inPolys;
  r.SetSizeFuncFromPoly(outPoly, inPolys, 1);
  VecPt3d pts = {{1, 1, 0}, {1, 9, 0}, {9, 9, 0}, {9, 1, 0}};
  VecDbl lengths;
  r.InterpEdgeLengths(pts, lengths);
  VecDbl baseLengths = {5.00892, 5.00892, 9.79527, 9.79527};
  TS_ASSERT_DELTA_VEC(baseLengths, lengths, 1e-5);
} // MePolyRedistributePtsUnitTests::testInterpEdgeLengths1
//------------------------------------------------------------------------------
/// \brief tests interpolating the edge length
//------------------------------------------------------------------------------
void MePolyRedistributePtsUnitTests::testInterpEdgeLengths2()
{
  MePolyRedistributePtsImpl r;
  // input polygon
  VecPt3d outPoly = {{0, 0, 0}, {0, 2, 0},  {0, 4, 0},   {0, 6, 0},
                     {0, 8, 0}, {0, 10, 0}, {10, 10, 0}, {10, 0, 0}};
  VecPt3d2d inPolys;
  VecPt3d in = {{4.9, 5, 0}, {4.95, 4.9, 0}, {5.05, 4.9, 0},
                {5.1, 5, 0}, {5.50, 5.1, 0}, {4.95, 5.1, 0}};
  inPolys.push_back(in);

  r.SetSizeFuncFromPoly(outPoly, inPolys, 1);
  VecPt3d pts = {{1, 1, 0}, {1, 9, 0}, {9, 9, 0}, {9, 1, 0}};
  VecDbl lengths;
  r.InterpEdgeLengths(pts, lengths);
  VecDbl baseLengths = {4.96661, 4.96586, 9.71356, 9.71552};
  TS_ASSERT_DELTA_VEC(baseLengths, lengths, 1e-5);
  r.m_sizeBias = 0;
  r.InterpEdgeLengths(pts, lengths);
  baseLengths = {3.36034, 3.36167, 7.02831, 7.03054};
  TS_ASSERT_DELTA_VEC(baseLengths, lengths, 1e-5);
} // MePolyRedistributePtsUnitTests::testInterpEdgeLengths2
//------------------------------------------------------------------------------
/// \brief tests interpolating the edge length
//------------------------------------------------------------------------------
void MePolyRedistributePtsUnitTests::testInterpEdgeLengths3()
{
  MePolyRedistributePtsImpl r;
  r.SetConstantSizeFunc(5.0);
  VecPt3d pts = {{0, 0, 0}, {0, 100, 0}, {100, 100, 0}, {100, 0, 0}};
  VecDbl lengths;
  r.InterpEdgeLengths(pts, lengths);
  VecDbl baseLengths(4, 5.0);
  TS_ASSERT_EQUALS_VEC(baseLengths, lengths);
} // MePolyRedistributePtsUnitTests::testInterpEdgeLengths3
//------------------------------------------------------------------------------
/// \brief tests interpolating the edge length
//------------------------------------------------------------------------------
void MePolyRedistributePtsUnitTests::testInterpEdgeLengths4()
{
  MePolyRedistributePtsImpl r;
  // r.SetConstantSizeFunc(5.0);
  VecPt3d pts = {{0, 0, 0}, {0, 100, 0}, {100, 100, 0}, {100, 0, 0}};
  VecDbl lengths;
  r.InterpEdgeLengths(pts, lengths);
  VecDbl baseLengths(4, 100.0);
  TS_ASSERT_EQUALS_VEC(baseLengths, lengths);
} // MePolyRedistributePtsUnitTests::testInterpEdgeLengths4
//------------------------------------------------------------------------------
/// \brief test redistributing the points on the polygon boundary
//------------------------------------------------------------------------------
void MePolyRedistributePtsUnitTests::testRedistPts()
{
  VecPt3d pts = {{0, 0, 0}, {0, 10, 0}, {10, 10, 0}, {10, 0, 0}, {0, 0, 0}};
  VecDbl lengths = {5, 5, 5, 5, 5};

  MePolyRedistributePtsImpl r;
  VecPt3d outPts = r.RedistPts(pts, lengths);
  VecPt3d basePts = {{0, 0, 0},  {0, 5, 0},  {0, 10, 0}, {5, 10, 0}, {10, 10, 0},
                     {10, 5, 0}, {10, 0, 0}, {5, 0, 0},  {0, 0, 0}};
  TS_ASSERT_EQUALS_VEC(basePts, outPts);
} // MePolyRedistributePtsUnitTests::testRedistPts
//------------------------------------------------------------------------------
/// \brief test redistributing the points on the polygon boundary
//------------------------------------------------------------------------------
void MePolyRedistributePtsUnitTests::testRedistPts1()
{
  VecPt3d pts = {{0, 0, 0}, {0, 10, 0}, {10, 10, 0}, {10, 0, 0}, {0, 0, 0}};
  VecDbl lengths = {4, 4, 4, 4, 4};

  MePolyRedistributePtsImpl r;
  VecPt3d outPts = r.RedistPts(pts, lengths);
  VecPt3d basePts = {{0, 0, 0},  {0, 4, 0},  {0, 8, 0}, {2, 10, 0}, {6, 10, 0}, {10, 10, 0},
                     {10, 6, 0}, {10, 2, 0}, {8, 0, 0}, {4, 0, 0},  {0, 0, 0}};
  TS_ASSERT_DELTA_VECPT3D(basePts, outPts, 1e-9);
} // MePolyRedistributePtsUnitTests::testRedistPts1
//------------------------------------------------------------------------------
/// \brief test redistributing the points on the polygon boundary
//------------------------------------------------------------------------------
void MePolyRedistributePtsUnitTests::testRedistPts2()
{
  VecPt3d pts = {{0, 0, 0}, {0, 10, 0}, {10, 10, 0}, {10, 0, 0}, {0, 0, 0}};
  VecDbl lengths = {3, 3, 3, 3, 3};

  MePolyRedistributePtsImpl r;
  VecPt3d outPts = r.RedistPts(pts, lengths);
  VecPt3d basePts = {{0, 0, 0},     {0, 3.08, 0},  {0, 6.15, 0},  {0, 9.23, 0},  {2.31, 10, 0},
                     {5.38, 10, 0}, {8.46, 10, 0}, {10, 8.46, 0}, {10, 5.38, 0}, {10, 2.31, 0},
                     {9.23, 0, 0},  {6.15, 0, 0},  {3.08, 0, 0},  {0, 0, 0}};
  TS_ASSERT_DELTA_VECPT3D(basePts, outPts, 0.01);
} // MePolyRedistributePtsUnitTests::testRedistPts2
//------------------------------------------------------------------------------
/// \brief test redistributing the points on the polygon boundary
//------------------------------------------------------------------------------
void MePolyRedistributePtsUnitTests::testRedistPts3()
{
  VecPt3d pts = {{0, 0, 0}, {0, 10, 0}, {10, 10, 0}, {10, 0, 0}, {0, 0, 0}};
  VecDbl lengths = {3, 5, 3, 5, 3};

  MePolyRedistributePtsImpl r;
  VecPt3d outPts = r.RedistPts(pts, lengths);
  VecPt3d basePts = {{0, 0, 0},     {0, 3.34, 0},  {0, 7.43, 0},  {2.2, 10, 0},
                     {6.36, 10, 0}, {9.76, 10, 0}, {10, 6.95, 0}, {10, 2.95, 0},
                     {8.19, 0, 0},  {3.95, 0, 0},  {0, 0, 0}};
  TS_ASSERT_DELTA_VECPT3D(basePts, outPts, 0.01);
} // MePolyRedistributePtsUnitTests::testRedistPts3
//------------------------------------------------------------------------------
/// \brief test redistributing the points on the polygon boundary
//------------------------------------------------------------------------------
void MePolyRedistributePtsUnitTests::testRedistPts4()
{
  VecPt3d pts = {{0, 0, 0}, {0, 10, 0}, {10, 10, 0}, {10, 0, 0}, {0, 0, 0}};
  VecDbl lengths = {15, 15, 15, 15, 15};

  MePolyRedistributePtsImpl r;
  VecPt3d outPts = r.RedistPts(pts, lengths);
  VecPt3d basePts = {{0, 0, 0}, {3.33, 10, 0}, {10, 3.33, 0}, {0, 0, 0}};
  TS_ASSERT_DELTA_VECPT3D(basePts, outPts, 0.01);
} // MePolyRedistributePtsUnitTests::testRedistPts4
//------------------------------------------------------------------------------
/// \brief test redistributing the points on the polygon boundary
//------------------------------------------------------------------------------
void MePolyRedistributePtsUnitTests::testRedistPts5()
{
  VecPt3d pts = {{0, 0, 0}, {0, 10, 0}, {10, 10, 0}, {10, 0, 0}, {0, 0, 0}};
  VecDbl lengths = {20, 20, 20, 20, 20};

  MePolyRedistributePtsImpl r;
  VecPt3d outPts = r.RedistPts(pts, lengths);
  VecPt3d basePts = {{0, 0, 0}, {3.33, 10, 0}, {10, 3.33, 0}, {0, 0, 0}};
  TS_ASSERT_DELTA_VECPT3D(basePts, outPts, 0.01);
} // MePolyRedistributePtsUnitTests::testRedistPts4
//------------------------------------------------------------------------------
/// \brief tests intersecting the polygon with the triangles from the interp
/// size function
//------------------------------------------------------------------------------
void MePolyRedistributePtsUnitTests::testIntersectWithTris()
{
  VecPt3d loop = {{0, 0, 0}, {0, 10, 0}, {10, 10, 0}, {10, 0, 0}};
  BSHP<VecPt3d> triPts(new VecPt3d());
  *triPts = {{5, 5, 0},   {-5, -5, 0}, {5, -5, 0},  {15, -5, 0}, {15, 5, 0},
             {15, 15, 0}, {5, 15, 0},  {-5, 15, 0}, {-5, 5, 0}};
  BSHP<VecInt> triTris(new VecInt());
  *triTris = {0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 5, 6, 0, 6, 7, 0, 7, 8, 0, 8, 1};

  MePolyRedistributePtsImpl r;
  r.m_sizePts = triPts;
  r.m_sizeTris = triTris;
  r.CreatePolyIntersector();
  r.IntersectWithTris(loop);
  VecPt3d baseLoop = {{0, 0, 0},   {0, 5, 0},  {0, 10, 0}, {5, 10, 0},
                      {10, 10, 0}, {10, 5, 0}, {10, 0, 0}, {5, 0, 0}};
  TS_ASSERT_DELTA_VECPT3D(baseLoop, loop, FLT_EPSILON);
} // MePolyRedistributePtsUnitTests::testIntersectWithTris
//------------------------------------------------------------------------------
/// \brief test redistributing the points on a polyline
//------------------------------------------------------------------------------
void MePolyRedistributePtsUnitTests::testRedistPolyLine()
{
  VecPt3d pts = {{0, 0, 0}, {0, 10, 0}, {10, 10, 0}, {10, 0, 0}};
  VecDbl lengths = {5, 5, 5, 5};

  MePolyRedistributePtsImpl r;
  VecPt3d outPts = r.RedistPts(pts, lengths);
  VecPt3d basePts = {{0, 0, 0},   {0, 5, 0},  {0, 10, 0}, {5, 10, 0},
                     {10, 10, 0}, {10, 5, 0}, {10, 0, 0}};
  TS_ASSERT_DELTA_VECPT3D(basePts, outPts, 0.01);
} // MePolyRedistributePtsUnitTests::testRedistPolyLine
//------------------------------------------------------------------------------
/// \brief test redistributing the points on a polyline
//------------------------------------------------------------------------------
void MePolyRedistributePtsUnitTests::testRedistPolyLine1()
{
  VecPt3d pts = {{0, 0, 0}, {0, 10, 0}, {10, 10, 0}, {10, 0, 0}};
  VecDbl lengths = {6, 6, 6, 6};

  MePolyRedistributePtsImpl r;
  VecPt3d outPts = r.RedistPts(pts, lengths);
  VecPt3d basePts = {{0, 0, 0}, {0, 6, 0}, {2, 10, 0}, {8, 10, 0}, {10, 6, 0}, {10, 0, 0}};
  TS_ASSERT_DELTA_VECPT3D(basePts, outPts, 0.01);
} // MePolyRedistributePtsUnitTests::testRedistPolyLine1
//------------------------------------------------------------------------------
/// \brief test redistributing the points on a polyline
//------------------------------------------------------------------------------
void MePolyRedistributePtsUnitTests::testRedistPolyLine2()
{
  VecPt3d pts = {{0, 0, 0}, {0, 10, 0}, {10, 10, 0}, {10, 0, 0}};
  VecDbl lengths = {3, 5, 3, 5};

  MePolyRedistributePtsImpl r;
  VecPt3d outPts = r.RedistPts(pts, lengths);
  VecPt3d basePts = {{0, 0, 0},     {0, 3.42, 0},  {0, 7.58, 0},  {2.43, 10, 0}, {6.66, 10, 0},
                     {10, 9.85, 0}, {10, 6.76, 0}, {10, 3.68, 0}, {10, 0, 0}};
  TS_ASSERT_DELTA_VECPT3D(basePts, outPts, 0.01);
} // MePolyRedistributePtsUnitTests::testRedistPolyLine2

//} // namespace xms
#endif
