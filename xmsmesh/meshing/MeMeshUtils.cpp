//------------------------------------------------------------------------------
/// \file
/// \brief Utilities meshing algorithms
/// \ingroup meshing
/// \copyright (C) Copyright Aquaveo 2018. Distributed under FreeBSD License
/// (See accompanying file LICENSE or https://aqaveo.com/bsd/license.txt)
//------------------------------------------------------------------------------

//----- Included files ---------------------------------------------------------

// 1. Precompiled header

// 2. My own header
#include <xmsmesh/meshing/MeMeshUtils.h>

// 3. Standard library headers
#include <map>
#include <sstream>

// 4. External library headers

// 5. Shared code headers
#include <xmscore/math/math.h>
#include <xmscore/misc/XmError.h>
#include <xmscore/misc/DynBitset.h>
#include <xmsinterp/triangulate/TrTin.h>

// 6. Non-shared code headers

//----- Forward declarations ---------------------------------------------------

//----- External globals -------------------------------------------------------

//----- Namespace declaration --------------------------------------------------
namespace xms
{
//----- Constants / Enumerations -----------------------------------------------

//----- Classes / Structs ------------------------------------------------------

//------------------------------------------------------------------------------
/// \brief Helper class for size function smoothing
//------------------------------------------------------------------------------
class SmoothIo
{
public:
  SmoothIo()
  : m_tin()
  , m_sizes(nullptr)
  , m_anchorType(0)
  , m_ptsFlag()
  , m_smoothSize(nullptr)
  , m_checkMinSize(false)
  , m_sizeRatio(0)
  , m_minSize(0)
  , m_scaleFactor(0)
  , m_percentGrowth(0)
  , m_logPercentGrowth(0)
  , m_maxSize(0)
  {
  }

  bool CalcMaxSize(double a_length, float a_smoothVal, double& a_maxSize);
  double CalcMinSize(double a_length, float a_smoothVal, double a_calcMaxSize);

  BSHP<TrTin> m_tin;     ///< geometry defining connections between points
  const VecFlt* m_sizes; ///< array of size values
  int m_anchorType;      ///< anchor to min or max value
  DynBitset m_ptsFlag;   ///< flags indicating whether a point should be processed
  VecFlt* m_smoothSize;  ///< the output smoothed sizes

  // used with size function smoothing
  bool m_checkMinSize;       ///< flag to indicate that the min size should be checked
  double m_sizeRatio;        ///< size ratio value
  double m_minSize;          ///< min size
  double m_scaleFactor;      ///< scaling factor
  double m_percentGrowth;    ///< growth factor
  double m_logPercentGrowth; ///< growth factor

  double m_maxSize; ///< max size used with elevation smoothing
};

//----- Internal functions -----------------------------------------------------
//------------------------------------------------------------------------------
/// \brief Calculates a max size for use in meiDoSmooth
/// \param[in] a_length Length between points being considered
/// \param[in] a_smoothVal Current size function value at point being evaluated
/// \param[out] a_maxSize The max size at the point being evaluated
/// \return True if no problems encountered
//------------------------------------------------------------------------------
bool SmoothIo::CalcMaxSize(double a_length, float a_smoothVal, double& a_maxSize)
{
  if (m_checkMinSize)
  {
    double sizeOneAway = a_length / (m_scaleFactor * a_smoothVal);
    double negLogCheck = sizeOneAway * (m_percentGrowth - 1.0) + 1.0;
    XM_ENSURE_TRUE(negLogCheck > 0.0, false);
    double numElems = log(negLogCheck) / m_logPercentGrowth;
    a_maxSize = pow(m_percentGrowth, numElems) * a_smoothVal;
  }
  else
  {
    a_maxSize = (double)a_smoothVal + a_length * m_maxSize;
  }
  return true;
} // SmoothIo::CalcMaxSize
//------------------------------------------------------------------------------
/// \brief Calculates a max size for use in meiDoSmooth
/// \param[in] a_length Length between points being considered
/// \param[in] a_smoothVal Current size function value at point being evaluated
/// \param[in] a_calcMaxSize The max size calculated at the point being evaluated
/// \return the min size
//------------------------------------------------------------------------------
double SmoothIo::CalcMinSize(double a_length, float a_smoothVal, double a_calcMaxSize)
{
  double minSize(0);
  if (m_checkMinSize)
  {
    minSize = a_smoothVal - (a_calcMaxSize - a_smoothVal);
  }
  else
  {
    minSize = (double)a_smoothVal - a_length * m_maxSize;
  }
  return minSize;
} // SmoothIo::CalcMaxSize
//------------------------------------------------------------------------------
/// \brief Smooths a size function. Ensures that the size function transitions
/// over a sufficient distance so that the area change of adjacent elements
/// does not exceed the size maximum passed in.
/// \param[in] a_ SmoothIo class with variables for the operation
//------------------------------------------------------------------------------
static void meiDoSmooth(SmoothIo& a_)
{
  if (a_.m_ptsFlag.empty())
    a_.m_ptsFlag.resize(a_.m_tin->NumPoints(), true);
  XM_ENSURE_TRUE(a_.m_ptsFlag.size() == (size_t)a_.m_tin->NumPoints());

  DynBitset ptsFlag(a_.m_tin->NumPoints(), false);

  VecFlt& smoothSize(*a_.m_smoothSize);
  const VecFlt& sz(*a_.m_sizes);
  // set the smoothed size equal to the incoming size
  smoothSize = sz;
  // sort points based on size
  std::multimap<float, int> mapSizeIdx;
  std::vector<std::multimap<float, int>::iterator> vecIt;
  float val(1.0);
  // if we are anchoring to the max size then make the map key the negative
  // of the size
  if (1 == a_.m_anchorType)
    val = -1.0;
  for (int i = 0; i < (int)sz.size(); ++i)
  {
    vecIt.push_back(mapSizeIdx.insert(std::make_pair(val * sz[i], i)));
  }

  // get references to the Tin
  const VecPt3d& pts(a_.m_tin->Points());
  const VecInt2d& trisAdjToPts(a_.m_tin->TrisAdjToPts());
  const VecInt& tris(a_.m_tin->Triangles());

  // iterate through the points
  auto it = mapSizeIdx.begin();
  auto endIt = mapSizeIdx.end();
  for (; it != endIt; ++it)
  {
    int i = it->second; // index of the point
    if (!a_.m_ptsFlag[i])
      continue;
    if (a_.m_checkMinSize && (double)smoothSize[i] < a_.m_minSize)
    {
      smoothSize[i] = (float)a_.m_minSize;
    }
    const VecInt adjTris(trisAdjToPts[i]); // triangles adjacent to the point
    for (const auto& t : adjTris)
    {
      int tIdx = t * 3;
      int triPts[3];
      // points in the adjacent triangle
      triPts[0] = tris[tIdx];
      triPts[1] = tris[tIdx + 1];
      triPts[2] = tris[tIdx + 2];
      for (int t1 = 0; t1 < 3; ++t1)
      {
        int ix = triPts[t1]; // index of triangle point
        // make sure this is not the vertex in the "i" loop
        if (ix == i)
          continue;
        // make sure we have not already adjusted this point
        if (ptsFlag[ix])
          continue;

        const Pt3d &p0(pts[i]), &p1(pts[ix]);
        // calculate what the min or max size can be based on how close
        // the elements are
        double length = Mdist(p0.x, p0.y, p1.x, p1.y);
        double maxSize(0);
        XM_ENSURE_TRUE(a_.CalcMaxSize(length, smoothSize[i], maxSize));

        switch (a_.m_anchorType)
        {
        case 0: // anchor to the min size
        {
          if (a_.m_checkMinSize)
            XM_ENSURE_TRUE(maxSize > a_.m_minSize);
          XM_ENSURE_TRUE(maxSize >= (double)smoothSize[i]);
          if (maxSize < (double)smoothSize[ix])
          {
            smoothSize[ix] = (float)maxSize;
            ptsFlag.set(ix, true);
          }
          else if (a_.m_checkMinSize && (double)smoothSize[ix] < a_.m_minSize)
          {
            smoothSize[ix] = (float)a_.m_minSize;
            ptsFlag.set(ix, true);
          }
        }
        break;
        case 1: // anchor to the max size
        {
          double minSize = a_.CalcMinSize(length, smoothSize[i], maxSize);
          // double minSize = smoothSize[i] - (maxSize - smoothSize[i]);
          XM_ENSURE_TRUE(minSize < smoothSize[i]);
          if (minSize > (double)smoothSize[ix])
          {
            smoothSize[ix] = (float)minSize;
            ptsFlag.set(ix, true);
          }
          else if (a_.m_checkMinSize && (double)smoothSize[ix] < a_.m_minSize)
          {
            smoothSize[ix] = (float)a_.m_minSize;
            ptsFlag.set(ix, true);
          }
        }
        break;
        default:
          XM_ASSERT(0);
          break;
        }
        // resort the point if it changed size
        if (ptsFlag[ix])
        {
          ptsFlag.set(ix, false);
          mapSizeIdx.erase(vecIt[ix]);
          vecIt[ix] = mapSizeIdx.insert(std::make_pair(val * smoothSize[ix], ix));
        }
      }
    }
  }
} // meiDoSmooth

//------------------------------------------------------------------------------
/// \brief Creates a size at each point based on the depth at the point and the
/// min and max sizes
/// the equation is  min_depth +
/// ( (depth - min_depth) / (max_depth - min_depth) ) * (max_size - min_size)
/// \param[in] a_depths The measured depths at point locations
/// \param[out] a_size A target size value based on the depth and the min/max
/// specified size.
/// \param[in] a_minSize The minimum element edge size
/// \param[in] a_maxSize The maximum element edge size
//------------------------------------------------------------------------------
void meSizeFunctionFromDepth(const VecDbl& a_depths,
                             VecDbl& a_size,
                             double a_minSize,
                             double a_maxSize)
{
  XM_ENSURE_TRUE(a_minSize > 0);
  XM_ENSURE_TRUE(a_maxSize > a_minSize);
  double minD = *std::min_element(a_depths.begin(), a_depths.end());
  double maxD = *std::max_element(a_depths.begin(), a_depths.end());
  double dDiff = maxD - minD;
  XM_ENSURE_TRUE(dDiff > 0);
  double sDiff = a_maxSize - a_minSize;
  XM_ENSURE_TRUE(sDiff > 0);

  a_size.assign(a_depths.size(), a_minSize);
  int i(0);
  for (auto& d : a_size)
  {
    d += ((a_depths[i++] - minD) / dDiff) * sDiff;
  }
} // meSizeFunctionFromDepth
//------------------------------------------------------------------------------
/// \brief Smooths a size function. Ensures that the size function transitions
/// over a sufficient distance so that the area change of adjacent elements
/// meets the size ratio passed in.
/// \param[in] a_tin Points and triangles defining the connectivity of the
/// size function.
/// \param[in] a_sizes Array of the current sizes
/// \param[in] a_sizeRatio Allowable size difference between adjacent elements
/// \param[in] a_minSize Minimum specified element size
/// \param[in] a_anchorType The minimum element edge size
/// \param[in] a_ptsFlag Flag to indicate if the value at the point should be
/// adjusted (a value of true will skip the point).
/// Leave the bitset empty to process all points.
/// \param[out] a_smoothSize Array of smoothed sizes
//------------------------------------------------------------------------------
void meSmoothSizeFunction(BSHP<TrTin> a_tin,
                          const VecFlt& a_sizes,
                          double a_sizeRatio,
                          double a_minSize,
                          int a_anchorType,
                          const DynBitset& a_ptsFlag,
                          VecFlt& a_smoothSize)
{
  a_smoothSize.resize(0);

  // error checking
  XM_ENSURE_TRUE(!a_sizes.empty());
  XM_ENSURE_TRUE(a_sizeRatio > 0.0);
  XM_ENSURE_TRUE(a_minSize > 0.0);
  double percentGrowth = sqrt(1.0 / a_sizeRatio);
  XM_ENSURE_TRUE(percentGrowth > 0.0);
  double logPercentGrowth = log(percentGrowth);
  double scaleFactor = (percentGrowth - 1.0) / logPercentGrowth;
  XM_ENSURE_TRUE(scaleFactor != 0.0);

  SmoothIo io;
  io.m_tin = a_tin;
  io.m_sizes = &a_sizes;
  io.m_anchorType = a_anchorType;
  io.m_ptsFlag = a_ptsFlag;
  io.m_smoothSize = &a_smoothSize;

  io.m_checkMinSize = true;
  io.m_sizeRatio = a_sizeRatio;
  io.m_minSize = a_minSize;
  io.m_percentGrowth = percentGrowth;
  io.m_logPercentGrowth = logPercentGrowth;
  io.m_scaleFactor = scaleFactor;
  meiDoSmooth(io);
} // meSmoothSizeFunction
//------------------------------------------------------------------------------
/// \brief Smooths a elevations based on max specified slope (a_maxSlope)
/// preserving either the min or max based on a_anchorType
/// \param[in] a_tin Points and triangles defining the connectivity of the
/// size function.
/// \param[in] a_elevs Array of the current elevations
/// \param[in] a_maxSlope Maximum allowable slope
/// \param[in] a_anchorType Indicates if you are anchoring to the top or bottom
/// of the slope
/// \param[in] a_ptsFlag Flag to indicate if the value at the point should be
/// adjusted (a value of true will skip the point).
/// Leave the bitset empty to process all points.
/// \param[out] a_smoothElevs Array of smoothed elevations
//------------------------------------------------------------------------------
void meSmoothElevBySlope(BSHP<TrTin> a_tin,
                         const VecFlt& a_elevs,
                         double a_maxSlope,
                         int a_anchorType,
                         const DynBitset& a_ptsFlag,
                         VecFlt& a_smoothElevs)
{
  a_smoothElevs.resize(0);

  // error checking
  XM_ENSURE_TRUE(!a_elevs.empty());
  XM_ENSURE_TRUE(a_maxSlope > 0.0);

  SmoothIo io;
  io.m_tin = a_tin;
  io.m_sizes = &a_elevs;
  io.m_anchorType = a_anchorType;
  io.m_ptsFlag = a_ptsFlag;
  io.m_smoothSize = &a_smoothElevs;

  io.m_maxSize = a_maxSlope;

  meiDoSmooth(io);
} // meSmoothSizeFunction
//------------------------------------------------------------------------------
/// \brief Prepends the polygon id as part of an error messge
/// \param[in] a_polyId The id of the polygon
/// \param[in,out] a_msg The error message
//------------------------------------------------------------------------------
void meModifyMessageWithPolygonId(int a_polyId, std::string& a_msg)
{
  if (a_polyId > -1)
  {
    std::stringstream ss;
    ss << "Error meshing polygon id: " << a_polyId << ". ";
    a_msg = ss.str() + a_msg;
  }
} // meModifyMessageWithPolygonId

} // namespace xms

#if CXX_TEST

#include <xmsmesh/meshing/MeMeshUtils.t.h>

#include <xmscore/testing/TestTools.h>
#include <xmsinterp/triangulate/TrTriangulatorPoints.h>

////////////////////////////////////////////////////////////////////////////////
/// \class MeMeshUtilsUnitTests
/// \brief Unit tests for MeMeshUtils.
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief Tests creating a size function from depth
//------------------------------------------------------------------------------
void MeMeshUtilsUnitTests::testSizeFuncFromDepth()
{
  // array of depths
  xms::VecDbl depths = {0, 5, 10, 20, 25, 5, 0};
  // array for the computed sizes
  xms::VecDbl elemSize;
  // set the value of the min and max element size
  double minElem(2), maxElem(102);
  // generate the size array
  xms::meSizeFunctionFromDepth(depths, elemSize, minElem, maxElem);
  // verify that the sizes are as expected
  xms::VecDbl baseElemSize = {2, 22, 42, 82, 102, 22, 2};
  TS_ASSERT_DELTA_VEC(baseElemSize, elemSize, 1e-9);
} // MeMeshUtilsUnitTests::testSizeFuncFromDepth
//------------------------------------------------------------------------------
/// \brief Tests size function smoothing
/// \verbatim
/// Input
///
///       100     100     100      100
///  20- *[8]    *[9]     *[10]    *[11]
///    |
///    |
///    |  1       100      100      100
///  10- *[4]     *[5]     *[6]     *[7]
///    |
///    |
///    |  100     100      100      100
///   0- *[0]     *[1]     *[2]     *[3]
///      0-------10-------20-------30
///
/// Output
///
///       4.4     7.9     11.3      14.8
///  20- *[8]    *[9]     *[10]    *[11]
///    |
///    |
///    |  1       4.4      7.9      11.3
///  10- *[4]     *[5]     *[6]     *[7]
///    |
///    |
///    |  4.4     5.9      9.3       12.8
///   0- *[0]     *[1]     *[2]     *[3]
///      0-------10-------20-------30
/// \endverbatim
//------------------------------------------------------------------------------
//! [snip_MeMeshUtilsTests::testSmoothSizeFunc]
void MeMeshUtilsUnitTests::testSmoothSizeFunc()
{
  BSHP<xms::VecPt3d> pts(new xms::VecPt3d());
  *pts = {{0, 0},   {10, 0},  {20, 0}, {30, 0},  {0, 10},  {10, 10},
          {20, 10}, {30, 10}, {0, 20}, {10, 20}, {20, 20}, {30, 20}};
  xms::VecFlt sizes(pts->size(), 100);
  sizes[4] = 1;
  BSHP<xms::VecInt> tris(new xms::VecInt());
  BSHP<xms::VecInt2d> adjTris(new xms::VecInt2d());
  xms::TrTriangulatorPoints tr(*pts, *tris, &*adjTris);
  tr.Triangulate();
  BSHP<xms::TrTin> tin = xms::TrTin::New();
  tin->SetGeometry(pts, tris, adjTris);

  xms::VecFlt vSmooth;
  xms::DynBitset ptFlags;
  xms::meSmoothSizeFunction(tin, sizes, 0.5, 1.0, 0, ptFlags, vSmooth);
  xms::VecFlt baseSmooth = {4.46f, 5.90f,  9.36f, 12.83f, 1.0f,   4.46f,
                            7.93f, 11.39f, 4.46f, 7.93f,  11.39f, 14.86f};
  TS_ASSERT_DELTA_VEC(baseSmooth, vSmooth, .1);
} // MeMeshUtilsUnitTests::testSmoothSizeFunc
//! [snip_MeMeshUtilsTests::testSmoothSizeFunc]
//------------------------------------------------------------------------------
/// \brief Tests size function smoothing
/// \verbatim
/// Input
///
///       1        1        1        1
///  20- *[8]     *[9]     *[10]    *[11]
///    |
///    |
///    |  100      1        1        1
///  10- *[4]     *[5]     *[6]     *[7]
///    |
///    |
///    |  1        1        1        1
///   0- *[0]     *[1]     *[2]     *[3]
///      0-------10-------20-------30
///
/// Output
///
///       95.5     96.5     89.6     86.1
///  20- *[8]     *[9]     *[10]    *[11]
///    |
///    |
///    |  100      96.5     93.0     89.6
///  10- *[4]     *[5]     *[6]     *[7]
///    |
///    |
///    |  95.5    95.0     91.6      88.1
///   0- *[0]     *[1]     *[2]     *[3]
///      0-------10-------20-------30
/// \endverbatim
//------------------------------------------------------------------------------
//! [snip_MeMeshUtilsTests::testSmoothSizeFunc1]
void MeMeshUtilsUnitTests::testSmoothSizeFunc1()
{
  BSHP<xms::VecPt3d> pts(new xms::VecPt3d());
  *pts = {{0, 0},   {10, 0},  {20, 0}, {30, 0},  {0, 10},  {10, 10},
          {20, 10}, {30, 10}, {0, 20}, {10, 20}, {20, 20}, {30, 20}};
  xms::VecFlt sizes(pts->size(), 1);
  sizes[4] = 100;
  BSHP<xms::VecInt> tris(new xms::VecInt());
  BSHP<xms::VecInt2d> adjTris(new xms::VecInt2d());
  xms::TrTriangulatorPoints tr(*pts, *tris, &*adjTris);
  tr.Triangulate();
  BSHP<xms::TrTin> tin = xms::TrTin::New();
  tin->SetGeometry(pts, tris, adjTris);

  xms::VecFlt vSmooth;
  xms::DynBitset ptFlags;
  xms::meSmoothSizeFunction(tin, sizes, 0.5, 1.0, 1, ptFlags, vSmooth);
  xms::VecFlt baseSmooth = {96.53f, 95.10f, 91.63f, 88.17f, 100.0f, 96.53f,
                            93.07f, 89.60f, 96.53f, 93.07f, 89.60f, 86.14f};
  TS_ASSERT_DELTA_VEC(baseSmooth, vSmooth, .1);
} // MeMeshUtilsUnitTests::testSmoothSizeFunc1
//! [snip_MeMeshUtilsTests::testSmoothSizeFunc1]

//------------------------------------------------------------------------------
/// \brief Tests size function smoothing
/// \verbatim
/// Input
///
///       100     100     100      100
///  20- *[8]    *[9]     *[10]    *[11]
///    |
///    |
///    |  1       100      100      100
///  10- *[4]     *[5]     *[6]     *[7]
///    |
///    |
///    |  100     100      100      100
///   0- *[0]     *[1]     *[2]     *[3]
///      0-------10-------20-------30
///
/// Output
///
///       6.0    11.0     16.0      21.0
///  20- *[8]    *[9]     *[10]    *[11]
///    |
///    |
///    |  1       6.0     11.0      16.0
///  10- *[4]     *[5]     *[6]     *[7]
///    |
///    |
///    |  6.0     8.0     13.0       18.0
///   0- *[0]     *[1]     *[2]     *[3]
///      0-------10-------20-------30
/// \endverbatim
//------------------------------------------------------------------------------
//! [snip_MeMeshUtilsTests::testSmoothSizeFunc2]
void MeMeshUtilsUnitTests::testSmoothSizeFunc2()
{
  BSHP<xms::VecPt3d> pts(new xms::VecPt3d());
  *pts = {{0, 0},   {10, 0},  {20, 0}, {30, 0},  {0, 10},  {10, 10},
          {20, 10}, {30, 10}, {0, 20}, {10, 20}, {20, 20}, {30, 20}};
  xms::VecFlt sizes(pts->size(), 100);
  sizes[4] = 1;
  BSHP<xms::VecInt> tris(new xms::VecInt());
  BSHP<xms::VecInt2d> adjTris(new xms::VecInt2d());
  xms::TrTriangulatorPoints tr(*pts, *tris, &*adjTris);
  tr.Triangulate();
  BSHP<xms::TrTin> tin = xms::TrTin::New();
  tin->SetGeometry(pts, tris, adjTris);

  xms::VecFlt vSmooth;
  xms::DynBitset ptsFlag;
  xms::meSmoothElevBySlope(tin, sizes, 0.5, 0, ptsFlag, vSmooth);
  xms::VecFlt baseSmooth = {6.00f,  8.07f,  13.07f, 18.07f, 1.0f,   6.00f,
                            11.00f, 16.00f, 6.00f,  11.00f, 16.00f, 21.00f};
  TS_ASSERT_DELTA_VEC(baseSmooth, vSmooth, .1);
} // MeMeshUtilsUnitTests::testSmoothSizeFunc2
//! [snip_MeMeshUtilsTests::testSmoothSizeFunc2]
//------------------------------------------------------------------------------
/// \brief Tests size function smoothing
/// \verbatim
/// Input
///
///       1        1        1        1
///  20- *[8]     *[9]     *[10]    *[11]
///    |
///    |
///    |  100      1        1        1
///  10- *[4]     *[5]     *[6]     *[7]
///    |
///    |
///    |  1        1        1        1
///   0- *[0]     *[1]     *[2]     *[3]
///      0-------10-------20-------30
///
/// Output
///
///       95.0     90.0     85.0     80.0
///  20- *[8]     *[9]     *[10]    *[11]
///    |
///    |
///    |  100      95.0     90.0     85.0
///  10- *[4]     *[5]     *[6]     *[7]
///    |
///    |
///    |  95.0    92.9     87.9      82.9
///   0- *[0]     *[1]     *[2]     *[3]
///      0-------10-------20-------30
/// \endverbatim
//------------------------------------------------------------------------------
//! [snip_MeMeshUtilsTests::testSmoothSizeFunc3]
void MeMeshUtilsUnitTests::testSmoothSizeFunc3()
{
  BSHP<xms::VecPt3d> pts(new xms::VecPt3d());
  *pts = {{0, 0},   {10, 0},  {20, 0}, {30, 0},  {0, 10},  {10, 10},
          {20, 10}, {30, 10}, {0, 20}, {10, 20}, {20, 20}, {30, 20}};
  xms::VecFlt sizes(pts->size(), 1);
  sizes[4] = 100;
  BSHP<xms::VecInt> tris(new xms::VecInt());
  BSHP<xms::VecInt2d> adjTris(new xms::VecInt2d());
  xms::TrTriangulatorPoints tr(*pts, *tris, &*adjTris);
  tr.Triangulate();
  BSHP<xms::TrTin> tin = xms::TrTin::New();
  tin->SetGeometry(pts, tris, adjTris);

  xms::VecFlt vSmooth;
  xms::DynBitset ptsFlag;
  xms::meSmoothElevBySlope(tin, sizes, 0.5, 1, ptsFlag, vSmooth);
  xms::VecFlt baseSmooth = {95.00f, 92.92f, 87.92f, 82.92f, 100.0f, 95.00f,
                            90.00f, 85.00f, 95.00f, 90.00f, 85.00f, 80.00f};
  TS_ASSERT_DELTA_VEC(baseSmooth, vSmooth, .1);
} // MeMeshUtilsUnitTests::testSmoothSizeFunc3
  //! [snip_MeMeshUtilsTests::testSmoothSizeFunc3]

#endif