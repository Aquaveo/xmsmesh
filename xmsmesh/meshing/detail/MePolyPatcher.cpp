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
#include <xmsmesh/meshing/detail/MePolyPatcher.h>

// 3. Standard library headers
#include <set>

// 4. External library headers

// 5. Shared code headers
#include <xmscore/math/math.h>
#include <xmscore/misc/DynBitset.h>
#include <xmscore/misc/XmError.h>
#include <xmscore/misc/XmLog.h>
#include <xmscore/misc/xmstype.h>
#include <xmsinterp/geometry/geoms.h>
#include <xmsinterp/geometry/GmPtSearch.h>
#include <xmsmesh/meshing/MeMeshUtils.h>

// 6. Non-shared code headers

//----- Forward declarations ---------------------------------------------------
class MePolyPatcherUnitTests;

//----- External globals -------------------------------------------------------

//----- Namespace declaration --------------------------------------------------
namespace xms
{
//----- Constants / Enumerations -----------------------------------------------

//----- Classes / Structs ------------------------------------------------------

//----- Internal functions -----------------------------------------------------

//----- Class / Function definitions -------------------------------------------
////////////////////////////////////////////////////////////////////////////////
/// \class MePolyPatcherImpl
/// \brief Generates an adaptive patch mesh from a polygon
////////////////////////////////////////////////////////////////////////////////
class MePolyPatcherImpl : public MePolyPatcher
{
  friend MePolyPatcherUnitTests; ///< tests the class
public:
  MePolyPatcherImpl()
  : m_prog()
  , m_n0(0)
  , m_n1(0)
  , m_dn(0)
  , m_np_end(0)
  , m_m0(0)
  , m_m1(0)
  , m_dm(0)
  , m_np_side(0)
  , m_ptSearch(GmPtSearch::New(true))
  , m_meshPts(new VecPt3d())
  , m_polyId(-1)
  {
    m_ptSearch->VectorThatGrowsToSearch(m_meshPts);
  }

  virtual bool MeshIt(int a_polyId,
                      const VecPt3d& a_outPoly,
                      const VecInt& a_polyCorners,
                      double a_xytol,
                      VecPt3d& a_points,
                      VecInt& a_cells) override;
  //------------------------------------------------------------------------------
  /// \brief sets the observer to report progress
  /// \param a_: The observer.
  //------------------------------------------------------------------------------
  void SetObserver(BSHP<Observer> a_) override { m_prog = a_; }

  size_t PolyCornerClosestToTopLeft(const VecPt3d& a_oPoly,
                                    const VecInt& a_idx,
                                    const Pt3d& a_topLeft);
  void PolyPtsToSides(const VecPt3d& a_outPoly, const VecInt& a_polyCorners);
  void OrderPoints();
  bool QuadPatchErrors(const VecPt3d& a_outPoly, const VecInt& a_polyCorners);
  void QuadPatch(const VecPt3d& a_outPoly, const VecInt& a_polyCorners);
  bool SetupSideInfo();
  bool AdjustNodesInCol3a(VecInt& a_nodesincol, int flag);
  void GetPercentages(int a_np, const VecPt3d& a_pts, VecDbl& pcnts);
  void PointsFromPercentages(const VecPt3d& a_pts,
                             const VecDbl& a_pcnts,
                             int a_np_new,
                             VecPt3d& a_ptsNew);
  void InterpolatePercentages(const VecDbl& a_pcnt_in, VecDbl& a_pcnt_out);
  void CalcPointPercentages();
  void GenerateMeshPts();
  void SetUpColumnForRectPatch3a(int a_j, VecPt3d& a_newcol);
  void GenerateMeshCells();
  void ValidateMeshCells();
  bool CellOverlapsAdj(int a_cellIdx, DynBitset& a_cellFlags);

  void CreateCellAdjacencyInfo();
  void GetCellPtInfo(int a_cellIdx, VecInt& a_ptIdxs, VecInt& a_adjCellIdx, VecInt2d& a_adjPts);

  BSHP<Observer> m_prog;   ///< observer to report mesh generation progress
  VecPt3d2d m_pts,         ///< polygon points divided into 4 sides
    m_side1pts,            ///< 1 side of polygon points
    m_side2pts,            ///< 1 side of polygon points
    m_end1pts,             ///< 1 side of polygon points
    m_end2pts;             ///< 1 side of polygon points
  int m_n0,                ///< variable refactored out of patch code
    m_n1,                  ///< variable refactored out of patch code
    m_dn,                  ///< variable refactored out of patch code
    m_np_end,              ///< variable refactored out of patch code
    m_m0,                  ///< variable refactored out of patch code
    m_m1,                  ///< variable refactored out of patch code
    m_dm,                  ///< variable refactored out of patch code
    m_np_side;             ///< variable refactored out of patch code
  VecInt m_nodesinrow,     ///< number of nodes in a row
    m_nodesincol;          ///< number of nodes in a column
  VecDbl m_side1pcnt,      ///< percentage of node location along side
    m_side2pcnt,           ///< percentage of node location along side
    m_end1pcnt,            ///< percentage of node location along side
    m_end2pcnt;            ///< percentage of node location along side
  VecDbl2d m_newside1pcnt, ///< percentage of node location along side
    m_newside2pcnt,        ///< percentage of node location along side
    m_newend1pcnt,         ///< percentage of node location along side
    m_newend2pcnt;         ///< percentage of node location along side

  BSHP<GmPtSearch> m_ptSearch; ///< spatial index for searching points
  BSHP<VecPt3d> m_meshPts;     ///< generated mesh nodes
  VecInt m_ptIdxToRemove;      ///< mesh nodes that will be removed
  VecInt m_meshCells;          ///< generated mesh cells
  VecInt2d m_nodeIdx;          ///< index to m_meshPts identifying mesh nodes
  double m_xytol;              ///< tolerance for geometry comparison
  VecInt m_ptCellCnt,          ///< number of cells connected to each point
    m_ptCellIdx,               ///< index into cell adjacency array for each point
    m_ptAdjCells,              ///< adjacent cells for each point
    m_cellIdx;                 ///< location of each cell in m_meshCells
  int m_polyId;                ///< id of the polygon
};

//------------------------------------------------------------------------------
/// \brief Creates a new instance of this class
/// \return Shared pointer to an instance of MePolyPatcher
//------------------------------------------------------------------------------
BSHP<MePolyPatcher> MePolyPatcher::New()
{
  BSHP<MePolyPatcher> ret(new MePolyPatcherImpl());
  return ret;
} // MePolyPaverToMeshPts::New
//------------------------------------------------------------------------------
/// \brief Perform MESH_PATCH meshing on a polygon. Return the
///        mesh points through a_points and a_cells (mix of quad/tri).
/// \param[in]  a_polyId:      Id of the polygon or -1
/// \param[in]  a_outPoly:     Array of points defining the polygon
/// \param[in]  a_polyCorners: Indices of corner points in a_outPoly
/// \param[in]  a_xytol:       Tolerance for geometric comparisons.
/// \param[out] a_points:      Computed mesh points.
/// \param[out] a_cells:       Computed mesh cells, a mix of quads and tris.
/// \return true on success.
//------------------------------------------------------------------------------
bool MePolyPatcherImpl::MeshIt(int a_polyId,
                               const VecPt3d& a_outPoly,
                               const VecInt& a_polyCorners,
                               double a_xytol,
                               VecPt3d& a_points,
                               VecInt& a_cells)
{
  m_polyId = a_polyId;
  m_xytol = a_xytol;
  size_t nCorner = a_polyCorners.size();
  // if all corner indices are the same we have a problem
  bool same(true);
  for (size_t i = 1; same && i < a_polyCorners.size(); ++i)
  {
    if (a_polyCorners[i] != a_polyCorners[0])
      same = false;
  }
  if (same)
    nCorner = 1;

  if (2 == nCorner || 3 == nCorner)
    QuadPatch(a_outPoly, a_polyCorners);
  else
  {
    std::string msg = "Polygon corners incorrectly specified. Aborting patch mesh generation.";
    meModifyMessageWithPolygonId(m_polyId, msg);
    XM_LOG(xmlog::error, msg);
    return false;
  }

  a_points.swap(*m_meshPts);
  a_cells.swap(m_meshCells);
  return true;
} // MePolyPatcherImpl::MeshIt
//------------------------------------------------------------------------------
/// \brief puts the polygon points in the
/// \param[in]  a_oPoly:      Array of points defining the polygon
/// \param[in]  a_idx:        Indices of corner points in a_outPoly
/// \param[out] a_topLeft:    The location of point closest to top left
/// \return index to the corner closest to the top left
//------------------------------------------------------------------------------
size_t MePolyPatcherImpl::PolyCornerClosestToTopLeft(const VecPt3d& a_oPoly,
                                                     const VecInt& a_idx,
                                                     const Pt3d& a_topLeft)
{
  int polyCornerIdx(0);
  double dist = MdistSq(a_oPoly[0].x, a_oPoly[0].y, a_topLeft.x, a_topLeft.y);
  for (size_t i = 0; i < a_idx.size(); ++i)
  {
    int q = a_idx[i];
    double d1 = MdistSq(a_oPoly[q].x, a_oPoly[q].y, a_topLeft.x, a_topLeft.y);
    if (d1 < dist)
    {
      dist = d1;
      polyCornerIdx = (int)i;
    }
  }
  return polyCornerIdx;
} // MePolyPatcherImpl::PolyCornerClosestToTopLeft
//------------------------------------------------------------------------------
/// \brief puts the polygon points in the
/// \param[in]  a_outPoly:      Array of points defining the polygon
/// \param[in]  a_polyCorners: Indices of corner points in a_outPoly
//------------------------------------------------------------------------------
void MePolyPatcherImpl::PolyPtsToSides(const VecPt3d& a_outPoly, const VecInt& a_polyCorners)
{
  XM_ENSURE_TRUE(!a_outPoly.empty());
  // create pt vectors (like the figure) from a_outPoly
  //                pts[3]
  //       |  ----------------->  |
  //       |                      |
  //       |                      |
  //       |                      |
  //       |                      |
  // pts[0] |                      | pts[2]
  //       |                      |
  //       |                      |
  //       |                      |
  //       |                      |
  //       |                      |
  //       v  ----------------->  v
  //                pts[1]

  Pt3d pMin, pMax;
  gmEnvelopeOfPts(a_outPoly, pMin, pMax);
  // find corner that is closest to the upper left of the poly envelope
  VecInt idx(1, 0);
  idx.insert(idx.end(), a_polyCorners.begin(), a_polyCorners.end());
  Pt3d topLeft(pMin);
  topLeft.y = pMax.y;
  size_t polyCornerIdx = PolyCornerClosestToTopLeft(a_outPoly, idx, topLeft);

  // fill in the arrays of points from a_outPoly
  VecPt3d pts(a_outPoly);
  idx.push_back((int)pts.size());
  pts.push_back(pts.front());
  VecPt3d2d tmpPts;
  tmpPts.assign(idx.size() - 1, VecPt3d());
  for (size_t i = 1; i < idx.size(); ++i)
  {
    for (int j = idx[i - 1]; j <= idx[i]; ++j)
    {
      tmpPts[i - 1].push_back(pts[j]);
    }
  }
  // this array tells the order of tmpPts to m_pts shown in the figure above
  int ix[4][4] = {{3, 2, 1, 0}, {0, 3, 2, 1}, {1, 0, 3, 2}, {2, 1, 0, 3}};
  m_pts.assign(4, VecPt3d());
  m_pts[0].swap(tmpPts[ix[polyCornerIdx][0]]);
  m_pts[1].swap(tmpPts[ix[polyCornerIdx][1]]);
  m_pts[2].swap(tmpPts[ix[polyCornerIdx][2]]);
  m_pts[3].swap(tmpPts[ix[polyCornerIdx][3]]);
  std::reverse(m_pts[0].begin(), m_pts[0].end());
  std::reverse(m_pts[1].begin(), m_pts[1].end());
} // MePolyPatcherImpl::PolyPtsToSides
//------------------------------------------------------------------------------
/// \brief checks that the points are ordered in the expected order
/// \verbatim
///                           pts[3]
///                  |  ----------------->  |
///                  |                      |
///                  |                      |
///                  |                      |
///                  |                      |
///           pts[0] |                      | pts[2]
///                  |                      |
///                  |                      |
///                  |                      |
///                  |                      |
///                  |                      |
///                  v  ----------------->  v
///                           pts[1]
///
///             However, this function makes sure that the following
///             conditions hold true:
///               1. np[0] <= np[2]
///               2. np[3] <= np[1]
///               3. if (np[2] - np[0]) == 0 then (np[3] - np[1]) == 0 also.
///                  Otherwise, edges & sides will be swapped to get
///                  pts[1] == pts[3].
///               4. (np[2] - np[0]) >= (np[1] - np[3]).
/// \endverbatim
//------------------------------------------------------------------------------
void MePolyPatcherImpl::OrderPoints()
{
  VecInt np(m_pts.size(), 0);
  for (size_t i = 0; i < m_pts.size(); ++i)
    np[i] = (int)m_pts[i].size();

  int tmpnp;
  // see if we need to switch ends and sides
  if ((np[2] - np[0] == 0) && (np[3] - np[1] != 0))
  { // swap end 1 with side 1
    tmpnp = np[0];
    np[0] = np[3];
    np[3] = tmpnp;
    m_pts[0].swap(m_pts[3]);
    // swap end 2 with side 2
    tmpnp = np[1];
    np[1] = np[2];
    np[2] = tmpnp;
    m_pts[1].swap(m_pts[2]);
  }
  // make the end with less nodes the right end
  if (np[0] > np[2])
  { // swap two ends
    tmpnp = np[0];
    np[0] = np[2];
    np[2] = tmpnp;
    m_pts[0].swap(m_pts[2]);
    // reorder two sides
    std::reverse(m_pts[1].begin(), m_pts[1].end());
    std::reverse(m_pts[3].begin(), m_pts[3].end());
  }
  // make side with less nodes the top side
  if (np[3] > np[1])
  { // swap two sides
    tmpnp = np[1];
    np[1] = np[3];
    np[3] = tmpnp;
    m_pts[1].swap(m_pts[3]);
    // reorder two ends
    std::reverse(m_pts[0].begin(), m_pts[0].end());
    std::reverse(m_pts[2].begin(), m_pts[2].end());
  }
  // make difference between columns more than that between rows
  if ((np[1] - np[3]) > (np[2] - np[0]))
  { // swap end 1 with side 1
    tmpnp = np[0];
    np[0] = np[3];
    np[3] = tmpnp;
    m_pts[0].swap(m_pts[3]);
    // swap end 2 with side 2
    tmpnp = np[1];
    np[1] = np[2];
    np[2] = tmpnp;
    m_pts[1].swap(m_pts[2]);
  }
} // MePolyPatcherImpl::OrderPoints
//------------------------------------------------------------------------------
/// \brief Checks the polygon to make sure that each side of the polygon has
/// the same number of segments.
/// \param[in]  a_outPoly:      Array of points defining the polygon
/// \param[in]  a_polyCorners: Indices of corner points in a_outPoly
/// \return true if errors were found.
//------------------------------------------------------------------------------
bool MePolyPatcherImpl::QuadPatchErrors(const VecPt3d& a_outPoly, const VecInt& a_polyCorners)
{
  bool err(false);
  XM_ENSURE_TRUE(a_outPoly.size() > 1, true);

  // the angles at the corners must be convex
  VecInt c;
  c.reserve(a_polyCorners.size() + 1);
  c.push_back(0);
  for (size_t i = 0; i < a_polyCorners.size(); ++i)
    c.push_back(a_polyCorners[i]);

  Pt3d p0, p1, p2;
  p0 = a_outPoly.back();
  p1 = a_outPoly.front();
  p2 = a_outPoly[1];
  double angle = gmAngleBetweenEdges(p0, p1, p2);
  if (angle >= XM_PI)
    err = true;
  for (size_t i = 1; i < c.size() && !err; ++i)
  {
    if (c[i] == c[i - 1])
      continue; // skip when we have a merged node
    if (c[i] - 1 < 0)
    {
      std::string msg = "Invalid corner detected in polygon. Aborting patch.";
      meModifyMessageWithPolygonId(m_polyId, msg);
      XM_LOG(xmlog::error, msg);
      return true;
    }
    p0 = a_outPoly[c[i] - 1];
    p1 = a_outPoly[c[i]];
    if (c[i] + 1 < (int)a_outPoly.size())
      p2 = a_outPoly[c[i] + 1];
    else
      p2 = a_outPoly[0];
    angle = gmAngleBetweenEdges(p0, p1, p2);
    if (angle >= XM_PI)
      err = true;
  }

  if (err)
  {
    std::string msg = "Non-convex corner detected in polygon. Aborting patch.";
    meModifyMessageWithPolygonId(m_polyId, msg);
    XM_LOG(xmlog::error, msg);
    return true;
  }
  // Check number of sides. There should be 3 or 4. If there are 3 we treat
  // the corner with the smallest angle as the "merged node."
  return false;
} // MePolyPatcherImpl::QuadPatchErrors
//------------------------------------------------------------------------------
/// \brief Creates mesh cells using the patch algorithm for triangles.
/// \param[in]  a_outPoly:      Array of points defining the polygon
/// \param[in]  a_polyCorners: Indices of corner points in a_outPoly
//------------------------------------------------------------------------------
void MePolyPatcherImpl::QuadPatch(const VecPt3d& a_outPoly, const VecInt& a_polyCorners)
{
  if (QuadPatchErrors(a_outPoly, a_polyCorners))
    return;
  PolyPtsToSides(a_outPoly, a_polyCorners);
  OrderPoints();
  if (!SetupSideInfo())
    return;
  GenerateMeshPts();
  GenerateMeshCells();
  ValidateMeshCells();
} // MePolyPatcherImpl::QuadPatch
//------------------------------------------------------------------------------
/// \brief Set up some variables based on the points being setup.
/// \return True on success.
//------------------------------------------------------------------------------
bool MePolyPatcherImpl::SetupSideInfo()
{
  // get info about sides 0, 2
  m_n0 = (int)m_pts[0].size();
  m_n1 = (int)m_pts[2].size();
  m_dn = m_n1 - m_n0;
  m_np_end = Mmax(m_n0, m_n1);
  // get info about sides
  m_m0 = (int)m_pts[3].size();
  m_m1 = (int)m_pts[1].size();
  m_np_side = Mmax(m_m0, m_m1);
  m_dm = m_m1 - m_m0;
  // set up these arrays
  // nodes  = (noderec***)malloc(np_end  * sizeof(noderec**));
  m_nodesinrow.assign(m_np_end, 0);
  m_nodesincol.assign(m_np_side, 0);
  // determine number of nodes in each column
  double tmpd;
  for (int i = 0; i < m_np_side; i++)
  {
    tmpd = ((double)i / (double)(m_np_side - 1) * m_dn);
    if (m_dn < 0)
      m_nodesincol[i] = m_n0 + (int)(tmpd - .5);
    else
      m_nodesincol[i] = m_n0 + (int)(tmpd + .5);
  }
  // determine number of nodes in each row
  for (int i = 0; i < m_np_end; i++)
  {
    tmpd = ((double)i / (double)(m_np_end - 1) * m_dm);
    if (m_dm < 0)
      m_nodesinrow[i] = m_m0 + (int)(tmpd - .5);
    else
      m_nodesinrow[i] = m_m0 + (int)(tmpd + .5);
  }

  if (!AdjustNodesInCol3a(m_nodesincol, 1))
    return false;
  if (!AdjustNodesInCol3a(m_nodesinrow, 0))
    return false;
  return true;
} // MePolyPatcherImpl::SetupSideInfo
//------------------------------------------------------------------------------
/// \brief performs adjustments to the nodes to make the correct desired
/// transitions from column to column. Both odd and even numbers of nodes
/// transition better to an even number of nodes.
/// \param[in]  a_nodesincol: Array of the number of nodes in each column
/// \param[in]  a_flag:       To adjust just one side, do this.
/// \return true if no errors encountered
//------------------------------------------------------------------------------
bool MePolyPatcherImpl::AdjustNodesInCol3a(VecInt& a_nodesincol, int a_flag)
{
  XM_ENSURE_TRUE(!a_nodesincol.empty(), false);

  int changes = a_nodesincol.back() - a_nodesincol[0];

  if (a_flag)
  {
    if (changes > 1)
    {
      for (size_t i = 1; i < a_nodesincol.size() - 1; ++i)
      {
        if (a_nodesincol[i] > a_nodesincol[i - 1])
        {
          if (a_nodesincol[i] & 1)
          {
            if (!(a_nodesincol[i - 1] & 1))
            {
              (a_nodesincol[i])--;
            }
          }
        }
      }
      // if both boundaries odd, make all columns odd */
      if (a_nodesincol[0] % 2 && a_nodesincol.back() % 2)
      {
        for (size_t i = 1; i < a_nodesincol.size() - 1; ++i)
        {
          if (a_nodesincol[i] < a_nodesincol.back() && a_nodesincol[i] % 2 == 0)
            (a_nodesincol[i])++;
        }
      }
    }
    else if (changes < -1)
    {
      std::string msg = "Error in polygon to patch algorithm.";
      meModifyMessageWithPolygonId(m_polyId, msg);
      XM_LOG(xmlog::debug, msg);
      return false;
    }
  }
  // handle degenerate case
  for (size_t i = 1; i < a_nodesincol.size() - 1; ++i)
    a_nodesincol[i] = Mmax(a_nodesincol[i], 2);
  return true;
} // MePolyPatcherImpl::AdjustNodesInCol3a
//------------------------------------------------------------------------------
/// \brief returns the real percentage for each point across the line.
/// \param[in]  a_np:     The number of points.
/// \param[in]  a_pts:    The point locations
/// \param[out] a_pcnts:  The percentage of each point
//------------------------------------------------------------------------------
void MePolyPatcherImpl::GetPercentages(int a_np, const VecPt3d& a_pts, VecDbl& a_pcnts)
{
  XM_ENSURE_TRUE(a_np > 0);
  XM_ENSURE_TRUE(a_np <= (int)a_pts.size() && a_np <= (int)a_pcnts.size());
  double len(0), delta(0);

  // get total length of points
  for (int i = 1; i < a_np; ++i)
    len += Mdist(a_pts[i].x, a_pts[i].y, a_pts[i - 1].x, a_pts[i - 1].y);
  a_pcnts[0] = 0;
  a_pcnts[a_np - 1] = 1.0;
  // loop through middle points
  for (int i = 1; i < a_np - 1; ++i)
  {
    delta += Mdist(a_pts[i].x, a_pts[i].y, a_pts[i - 1].x, a_pts[i - 1].y);
    // AKZ
    if (len < 0.00001)
      a_pcnts[i] = (double)i / (double)(a_pts.size() - 1);
    else
      a_pcnts[i] = delta / len;
  }
} // MePolyPatcherImpl::GetPercentages
//------------------------------------------------------------------------------
/// \brief reinterpolates the array of points a_pts to a new array of points
/// a_ptsNew. The interpolation is performed based on the spacing of the
/// original a_pts.
/// \param[in]  a_pts:    The point locations
/// \param[in]  a_pcnt:   The percentage of each point in a_pts
/// \param[in]  a_np_new: Number of new points
/// \param[out] a_ptsNew: The interpolated point locations
//------------------------------------------------------------------------------
void MePolyPatcherImpl::PointsFromPercentages(const VecPt3d& a_pts,
                                              const VecDbl& a_pcnt,
                                              int a_np_new,
                                              VecPt3d& a_ptsNew)
{
  // get percentages for new column/row
  VecDbl pcnt_new(a_np_new, 0);
  // degenerate case
  if (1 == a_pts.size())
  {
    for (size_t i = 0; i < pcnt_new.size() && i < a_ptsNew.size(); ++i)
    {
      a_ptsNew[i] = a_pts[0];
      pcnt_new[i] = 1.0;
    }
  }
  else
  {
    InterpolatePercentages(a_pcnt, pcnt_new);
    for (size_t i = 0, j = 0; i < pcnt_new.size(); i++)
    { // get the index before or on current percentage
      while (j < (a_pts.size() - 1) && a_pcnt[j + 1] <= pcnt_new[i])
        j++;

      if (j == (a_pts.size() - 1))
        a_ptsNew[i] = a_pts[j];
      else
      { // get values of current percentage
        a_ptsNew[i].x = a_pts[j].x + (pcnt_new[i] - a_pcnt[j]) * (a_pts[j + 1].x - a_pts[j].x) /
                                       (a_pcnt[j + 1] - a_pcnt[j]);
        a_ptsNew[i].y = a_pts[j].y + (pcnt_new[i] - a_pcnt[j]) * (a_pts[j + 1].y - a_pts[j].y) /
                                       (a_pcnt[j + 1] - a_pcnt[j]);
        a_ptsNew[i].z = a_pts[j].z + (pcnt_new[i] - a_pcnt[j]) * (a_pts[j + 1].z - a_pts[j].z) /
                                       (a_pcnt[j + 1] - a_pcnt[j]);
      }
    }
  }
} // MePolyPatcherImpl::PointsFromPercentages
//------------------------------------------------------------------------------
/// \brief gets the new percentages, interpolated from old percentages
///
///                             |
///    y-axis = % along string  |                          ___+
///      (from the function     |                    _+___/
///       paiGetPercentages())  |               ____/
///                             |           __+/
///                             |        __/
///                             |     _+/
///                             |    /
///                             |  _/
///                             | /
///                             |/_____________________________
///                                    |      |       |       |
///                           node1  node2   node3   node4  node5
///
///                      x-axis = % of points (always an even distribution)
///
/// \param[in]  a_pcnt_in:    The input percentages
/// \param[in]  a_pcnt_out:   The input percentages
//------------------------------------------------------------------------------
void MePolyPatcherImpl::InterpolatePercentages(const VecDbl& a_pcnt_in, VecDbl& a_pcnt_out)
{
  // fill in the x-axis arrays
  VecDbl x_in(a_pcnt_in.size(), 0);
  for (size_t i = 1; i < x_in.size() - 1; ++i)
  {
    x_in[i] = (double)i / (x_in.size() - 1);
  }
  x_in.back() = 1.0;

  VecDbl x_out(a_pcnt_out.size(), 0);
  for (size_t i = 1; i < x_out.size() - 1; ++i)
  {
    x_out[i] = (double)i / (x_out.size() - 1);
  }
  x_out.back() = 1.0;
  // find the y-axis array to fill in
  // DSG: reset end percentages to be whatever is passed in
  // DSG   pcnt_out[0]        = 0.0;
  // DSG   pcnt_out[np_out-1] = 1.0;
  a_pcnt_out.front() = x_out.front();
  a_pcnt_out.back() = x_out.back();

  for (size_t i = 1, j = 0; i < a_pcnt_out.size() - 1; ++i)
  { // get the index before or on current x-loc
    while (j < a_pcnt_in.size() - 1 && x_in[j + 1] <= x_out[i])
      j++;
    // get percentage of current x-loc
    a_pcnt_out[i] = a_pcnt_in[j] + (x_out[i] - x_in[j]) * (a_pcnt_in[j + 1] - a_pcnt_in[j]) /
                                     (x_in[j + 1] - x_in[j]);
  }
} // MePolyPatcherImpl::InterpolatePercentages
//------------------------------------------------------------------------------
/// \brief Calculates percentages that are used to determine the locations of
/// the mesh points
//------------------------------------------------------------------------------
void MePolyPatcherImpl::CalcPointPercentages()
{
  m_side1pcnt.assign(m_m0, 0);
  m_side2pcnt.assign(m_m1, 0);
  m_end1pcnt.assign(m_n0, 0);
  m_end2pcnt.assign(m_n1, 0);
  GetPercentages(m_m0, m_pts[3], m_side1pcnt);
  GetPercentages(m_m1, m_pts[1], m_side2pcnt);
  GetPercentages(m_n0, m_pts[0], m_end1pcnt);
  GetPercentages(m_n1, m_pts[2], m_end2pcnt);

  VecDbl vd(m_np_side, 0);
  VecPt3d tmp(m_np_side, Pt3d());
  m_side1pts.assign(m_np_end, tmp);
  m_side2pts.assign(m_np_end, tmp);
  m_newside1pcnt.assign(m_np_end, vd);
  m_newside2pcnt.assign(m_np_end, vd);

  vd.resize(m_np_end, 0);
  tmp.resize(m_np_end, Pt3d());
  m_end1pts.assign(m_np_side, tmp);
  m_end2pts.assign(m_np_side, tmp);
  m_newend1pcnt.assign(m_np_side, vd);
  m_newend2pcnt.assign(m_np_side, vd);

  // get points for each column percentages
  for (int i = 0; i < m_np_end; ++i)
  {
    PointsFromPercentages(m_pts[3], m_side1pcnt, m_nodesinrow[i], m_side1pts[i]);
    GetPercentages(m_nodesinrow[i], m_side1pts[i], m_newside1pcnt[i]);
    PointsFromPercentages(m_pts[1], m_side2pcnt, m_nodesinrow[i], m_side2pts[i]);
    GetPercentages(m_nodesinrow[i], m_side2pts[i], m_newside2pcnt[i]);
  }
  // get points for each row percentages
  for (int j = 0; j < m_np_side; ++j)
  {
    PointsFromPercentages(m_pts[0], m_end1pcnt, m_nodesincol[j], m_end1pts[j]);
    GetPercentages(m_nodesincol[j], m_end1pts[j], m_newend1pcnt[j]);
    PointsFromPercentages(m_pts[2], m_end2pcnt, m_nodesincol[j], m_end2pts[j]);
    GetPercentages(m_nodesincol[j], m_end2pts[j], m_newend2pcnt[j]);
  }
} // MePolyPatcherImpl::CalcPointPercentages
//------------------------------------------------------------------------------
/// \brief Creates the mesh pts.
//------------------------------------------------------------------------------
void MePolyPatcherImpl::GenerateMeshPts()
{
  CalcPointPercentages();
  // init m_nodeIdx to -1
  VecInt vn(m_np_side, -1);
  m_nodeIdx.assign(m_np_end, vn);
  // create the first column
  int idx;
  for (int i = 0, j = 0; i < m_nodesincol[j]; ++i)
  {
    m_ptSearch->AddPtToVectorIfUnique(m_pts[0][i], XM_ZERO_TOL, idx);
    m_nodeIdx[i][j] = (int)idx;
  }
  // create the last column
  for (int i = 0, j = (int)m_nodesincol.size() - 1; i < m_nodesincol[j]; ++i)
  {
    m_ptSearch->AddPtToVectorIfUnique(m_pts[2][i], XM_ZERO_TOL, idx);
    m_nodeIdx[i][j] = (int)idx;
  }
  // loop through columns and create mesh points
  size_t m0(m_pts[3].size());
  VecPt3d newcol(m_np_end, Pt3d());
  for (int j = 1; j < m_np_side - 1; ++j)
  {
    SetUpColumnForRectPatch3a(j, newcol);
    // top node
    newcol[0] = m_pts[3][(int)((double)j / (m_np_side - 1) * (m0 - 1) + 0.5)];
    // bottom node
    newcol[m_nodesincol[j] - 1] = m_pts[1][j];
    // create the nodes
    for (int i = 0; i < m_nodesincol[j]; i++)
    {
      Pt3d pt1, pt2;
      int idx1 = m_nodeIdx[i][j - 1];
      if (idx1 > 0 && idx1 < (int)m_meshPts->size())
        pt1 = (*m_meshPts)[idx1];
      int idx2 = m_nodeIdx[m_nodesincol[j - 1] - 1][j - 1];
      if (idx2 > 0 && idx2 < (int)m_meshPts->size())
        pt2 = (*m_meshPts)[idx2];

      if (fabs(newcol[i].x + 1234.5) < XM_ZERO_TOL && fabs(newcol[i].y + 1234.5) < XM_ZERO_TOL &&
          fabs(newcol[i].z + 1234.5) < XM_ZERO_TOL)
      {
        m_nodeIdx[i][j] = m_nodeIdx[0][0];
      }
      else if (i > 0 &&
               Mdist(newcol[i].x, newcol[i].y, newcol[i - 1].x, newcol[i - 1].y) < XM_ZERO_TOL)
      {
        m_nodeIdx[i][j] = m_nodeIdx[i - 1][j];
      }
      else if (i == 0 && Mdist(newcol[i].x, newcol[i].y, pt1.x, pt1.y) < XM_ZERO_TOL)
      {
        m_nodeIdx[i][j] = m_nodeIdx[i][j - 1];
      }
      else if (i == m_nodesincol[j] && Mdist(newcol[i].x, newcol[i].y, pt2.x, pt2.y) < XM_ZERO_TOL)
      {
        m_nodeIdx[i][j] = m_nodeIdx[m_nodesincol[j - 1] - 1][j - 1];
      }
      else
      {
        m_ptSearch->AddPtToVectorIfUnique(newcol[i], XM_ZERO_TOL, idx);
        m_nodeIdx[i][j] = (int)idx;
      }
    }
  }
} // MePolyPatcherImpl::GenerateMeshPts
//------------------------------------------------------------------------------
/// \brief initialize an array of nodes for a column of patch nodes
/// \param[in]  a_j:         The column index
/// \param[out] a_newcol:    The mesh point locations in a_j column
//------------------------------------------------------------------------------
void MePolyPatcherImpl::SetUpColumnForRectPatch3a(int a_j, VecPt3d& a_newcol)
{
  a_newcol.resize(m_np_end, Pt3d());
  int np = m_nodesincol[a_j];
  int np_side = m_nodesinrow.back();
  // get corners
  Pt3d c[4];
  c[0] = m_end1pts[a_j][0];
  c[1] = m_end1pts[a_j][np - 1];
  c[2] = m_end2pts[a_j][0];
  c[3] = m_end2pts[a_j][np - 1];
  // find points on new column
  for (int i = 1; i < np - 1; ++i)
  {
    double u_pcnt = (double)a_j / (np_side - 1);
    int row = (int)((double)i / (np - 1) * (m_np_end - 1) + 0.5);
    int col = (int)(u_pcnt * (m_nodesinrow[row] - 1) + 0.5);
    // avoid messing up the mesh, find least row that gives same row
    for (int i_test = i; row > 0 && i_test == i;)
    {
      i_test = (int)((double)(row - 1) / (m_np_end - 1) * (np - 1) + 0.5);
      if (i_test == i)
      {
        row--;
        col = (int)(u_pcnt * (m_nodesinrow[row] - 1) + 0.5);
      }
    }
    // don't create extra nodes on boundary
    if (col == m_nodesinrow[row] - 1 && col == 1)
    {
      a_newcol[i].x = a_newcol[i].y = a_newcol[i].z = -1234.5;
    }
    else
    {
      // don't create extra nodes on boundary
      if (col == 0)
        col++;
      else if (col == m_nodesinrow[row] - 1)
        col--;
      a_newcol[0] = m_side1pts[row][col];
      a_newcol[np - 1] = m_side2pts[row][col];
      // avoid messing up the mesh, find least j that gives same column
      int j_test, col_test;
      for (j_test = a_j, col_test = col; j_test > 1 && col_test == col;)
      {
        col_test =
          (int)((double)(j_test - 1) / (m_nodesinrow[m_np_end - 1] - 1) * (m_nodesinrow[row] - 1) +
                0.5);
        if (col_test == 0)
          col_test++;
        else if (col_test == m_nodesinrow[row] - 1)
          col_test--;
        if (col_test == col)
          j_test--;
      }
      // find correct i value
      int i_test = (int)((double)i / (m_nodesincol[a_j] - 1) * (m_nodesincol[j_test] - 1) + 0.5);
      // handle degenerate case with this if statement
      if (m_nodesincol[j_test] > 2)
      {
        if (i_test == 0)
          i_test++;
        else if (i_test == m_nodesincol[j_test] - 1)
          i_test--;
      }
      // compute these
      int j_min = Mmin(j_test, m_nodesinrow[i_test] - 1);
      u_pcnt = (double)j_min / (double)(m_nodesinrow[i_test] - 1);         // j_min
      double v_pcnt = (double)i_test / (double)(m_nodesincol[j_test] - 1); // j_test
      // now get the correct percentages
      double u1 =
        m_newside1pcnt[i_test][j_min] * (1.0 - v_pcnt) + m_newside2pcnt[i_test][j_min] * v_pcnt;
      double u2 = 1.0 - u1;
      double v1 =
        m_newend1pcnt[j_test][i_test] * (1.0 - u_pcnt) + m_newend2pcnt[j_test][i_test] * u_pcnt;
      double v2 = 1.0 - v1;
      // finally, compute the x, y, z components from the Coons equation
      a_newcol[i].x = u2 * m_end1pts[j_test][i_test].x + u1 * m_end2pts[j_test][i_test].x +
                      v2 * a_newcol[0].x + v1 * a_newcol[np - 1].x -
                      u2 * (v2 * c[0].x + v1 * c[1].x) - u1 * (v2 * c[2].x + v1 * c[3].x);
      a_newcol[i].y = u2 * m_end1pts[j_test][i_test].y + u1 * m_end2pts[j_test][i_test].y +
                      v2 * a_newcol[0].y + v1 * a_newcol[np - 1].y -
                      u2 * (v2 * c[0].y + v1 * c[1].y) - u1 * (v2 * c[2].y + v1 * c[3].y);
      a_newcol[i].z = u2 * m_end1pts[j_test][i_test].z + u1 * m_end2pts[j_test][i_test].z +
                      v2 * a_newcol[0].z + v1 * a_newcol[np - 1].z -
                      u2 * (v2 * c[0].z + v1 * c[1].z) - u1 * (v2 * c[2].z + v1 * c[3].z);
    }
  }
} // MePolyPatcherImpl::SetUpColumnForRectPatch3a
//------------------------------------------------------------------------------
/// \brief Creates the mesh cells and puts the mesh pts into 1 vector.
//------------------------------------------------------------------------------
void MePolyPatcherImpl::GenerateMeshCells()
{
  // These correspond to defines in vtkCellType. There are classes downstream from
  // meshing that convert a cell stream into a vtkUnstructured grid that use these
  // magic numbers.
  const int VTK_QUAD(9);
  const int VTK_TRI(5);

  // loop through columns
  for (int j = 1; j < m_np_side; ++j)
  {
    // find nodes for next element as follows
    //------------------------------------------------------------------------------
    //          node1a     node2a
    //              *-------*
    //              |       |
    //              |       |
    //              |       |
    //              *-------*
    //          node1b     node2b
    //------------------------------------------------------------------------------
    int i(0), index1(0), index2(0);
    for (i = 1, index1 = 0; i < m_nodesincol[j]; ++i)
    { // find index2
      if (m_nodesincol[j] == m_nodesincol[j - 1])
        index2 = i;
      else
      {
        double fval = (double)i / (m_nodesincol[j] - 1);
        if (fval >= 0.50)
          index2 = (int)(fval * (m_nodesincol[j - 1] - 1) + 0.51);
        else
          index2 = (int)(fval * (m_nodesincol[j - 1] - 1) + 0.49);
        // make sure index2 is okay
        index2 = Mmin(index2, index1 + 1);
      }
      // get nodes for new cell
      int node1a(m_nodeIdx[index1][j - 1]), node1b(m_nodeIdx[index2][j - 1]),
        node2a(m_nodeIdx[i - 1][j]), node2b(m_nodeIdx[i][j]);
      index1 = index2;
      // see if the element got pinched to a line
      if ((node1a == node1b && node2a == node2b) || (node1a == node2a && node1b == node2b) ||
          (node1a == node1b && node1a == node2a) || (node1b == node1a && node1b == node2b) ||
          (node2b == node1b && node2b == node2a) || (node2a == node1a && node2a == node2b))
      { // do not create element in this case
      }
      else
      { // make sure nodes are ccw
        int x[4] = {node1a, node1b, node2b, node2a};
        VecInt n_ccw(&x[0], &x[4]);
        VecInt idxs(4);
        VecPt3d pts(4, Pt3d());
        for (int i = 0; i < 4; ++i)
        {
          idxs[i] = i;
          pts[i] = (*m_meshPts)[n_ccw[i]];
        }
        gmOrderPointsCounterclockwise(pts, idxs);
        node1a = n_ccw[idxs[0]];
        node1b = n_ccw[idxs[1]];
        node2b = n_ccw[idxs[2]];
        node2a = n_ccw[idxs[3]];
        // see if new quad will be built
        if (node1a != node1b && node2a != node2b && node1a != node2a && node1b != node2b)
        {
          m_meshCells.push_back(VTK_QUAD); // cell type
          m_meshCells.push_back(4);        // number of points
          m_meshCells.push_back(node1a);   // the points that make up the
          m_meshCells.push_back(node1b);   // cell
          m_meshCells.push_back(node2b);
          m_meshCells.push_back(node2a);
        }
        else
        { // triangle will be built -- find which one
          int ix[3];
          ix[0] = node1a;
          if (node1a == node1b)
          {
            ix[1] = node2b;
            ix[2] = node2a;
          }
          else if (node2a == node2b)
          {
            ix[1] = node1b;
            ix[2] = node2a;
          }
          else if (node1a == node2a)
          {
            ix[1] = node1b;
            ix[2] = node2b;
          }
          else
          {
            ix[1] = node1b;
            ix[2] = node2a;
          }
          m_meshCells.push_back(VTK_TRI);
          m_meshCells.push_back(3);
          for (int i = 0; i < 3; ++i)
            m_meshCells.push_back(ix[i]);
        }
      }
    }
  }
} // MePolyPatcherImpl::GenerateMeshCells
//------------------------------------------------------------------------------
/// \brief Ensures that the generated mesh cells are valid (no ill formed cells)
//------------------------------------------------------------------------------
void MePolyPatcherImpl::ValidateMeshCells()
{
  if (m_meshCells.size() < 7)
    return;

  CreateCellAdjacencyInfo();
  bool err(false);
  int nCell(0);
  size_t i = 0;
  while (i < m_meshCells.size())
  {
    ++i; // celltype
    int numPts = m_meshCells[i++];
    i += numPts;
    ++nCell;
  }

  DynBitset cellFlags;
  cellFlags.resize(nCell, false);
  // loop through cells
  for (int i = 0; !err && i < nCell; ++i)
  {
    err = CellOverlapsAdj(i, cellFlags);
  }
  if (err)
  {
    m_meshCells.clear();
    m_meshPts->clear();
  }
} // MePolyPatcherImpl::ValidateMeshCells
//------------------------------------------------------------------------------
/// \brief Checks to see if the cells overlaps with any adjacent cells
/// \param[in] a_cellIdx Cell index
/// \param[in] a_cellFlags Flag to tell if cell has already been checked
/// \return True if the cell overlaps any adjacent cell
//------------------------------------------------------------------------------
bool MePolyPatcherImpl::CellOverlapsAdj(int a_cellIdx, DynBitset& a_cellFlags)
{
  const VecPt3d& mp(*m_meshPts);
  VecInt cellPtIdx, adjCellIdx;
  VecInt2d adjPtIdx;
  GetCellPtInfo(a_cellIdx, cellPtIdx, adjCellIdx, adjPtIdx);
  VecPt3d cellPoints;
  Pt3d cellCenter;
  for (size_t i = 0; i < cellPtIdx.size(); ++i)
  {
    cellPoints.push_back(mp[cellPtIdx[i]]);
    cellCenter += cellPoints.back();
  }
  cellCenter /= cellPoints.size();

  // add the first idx to the end of the vector so that we can test the edge defined by
  // the last point to the first point
  cellPtIdx.push_back(cellPtIdx.front());
  // check each edge to see if it overlaps the edges of any neighbor cell
  for (size_t i = 1; i < cellPtIdx.size(); ++i)
  {
    int i0(cellPtIdx[i - 1]), i1(cellPtIdx[i]);
    const Pt3d &p0(mp[i0]), &p1(mp[i1]);
    int cnt(-1);
    for (auto v : adjPtIdx)
    {
      cnt++;
      if (a_cellFlags[adjCellIdx[cnt]])
        continue;
      Pt3d adjCellCenter;
      VecPt3d adjacentPoints;
      for (size_t j = 0; j < v.size(); ++j)
      {
        adjacentPoints.push_back(mp[v[j]]);
        adjCellCenter += adjacentPoints.back();
      }
      adjCellCenter /= adjacentPoints.size();

      bool overlap(false);
      // see if centroid of a_cellIdx is inside of neigh
      if (gmPointInPolygon2D(&adjacentPoints[0], adjacentPoints.size(), cellCenter) > -1)
      {
        std::string msg = "Invalid patch. Centroid of base cell inside of adjacent cell.";
        meModifyMessageWithPolygonId(m_polyId, msg);
        XM_LOG(xmlog::error, msg);
        return true;
      }

      // see if centroid of neigh is inside of a_cellIdx
      if (gmPointInPolygon2D(&cellPoints[0], cellPoints.size(), adjCellCenter) > -1)
      {
        std::string msg = "Invalid patch. Centroid of adjacent cell inside of base cell.";
        meModifyMessageWithPolygonId(m_polyId, msg);
        XM_LOG(xmlog::error, msg);
        return true;
      }

      // add the first idx to the end of the vector so that we can test the edge defined by
      // the last point to the first point
      v.push_back(v.front());
      for (size_t j = 1; !overlap && j < v.size(); ++j)
      {
        int j0(v[j - 1]), j1(v[j]);
        // edges share a point in the grid so don't check if they intersect
        if (i0 == j0 || i0 == j1 || i1 == j0 || i1 == j1)
          continue;

        const Pt3d &p2(mp[j0]), &p3(mp[j1]);
        if (gmLinesIntersect(p0, p1, p2, p3))
        {
          std::string msg = "Invalid patch. Edges of adjacent cells overlap.";
          meModifyMessageWithPolygonId(m_polyId, msg);
          XM_LOG(xmlog::error, msg);
          return true;
        }
      }
    }
  }
  a_cellFlags[a_cellIdx] = true;
  return false;
} // MePolyPatcherImpl::CellOverlapsAdj
//------------------------------------------------------------------------------
/// \brief Creates arrays with cell adjacency information
//------------------------------------------------------------------------------
void MePolyPatcherImpl::CreateCellAdjacencyInfo()
{
  m_ptCellCnt.assign(m_meshPts->size(), 0);
  m_cellIdx.reserve(m_meshCells.size() / 3);
  // for each point count the number of attached cells
  for (size_t i = 0; i < m_meshCells.size(); ++i)
  {
    m_cellIdx.push_back(static_cast<int>(i));
    for (int j = 0; j < m_meshCells[i + 1]; ++j)
    {
      m_ptCellCnt[m_meshCells[i + 2 + j]]++;
    }
    i += (m_meshCells[i + 1] + 1);
  }
  // for each point create an index into an array with adjacent cells listed
  m_ptCellIdx.assign(m_ptCellCnt.size(), 0);
  int last(0);
  for (size_t i = 0; i < m_ptCellCnt.size(); ++i)
  {
    m_ptCellIdx[i] = last;
    last += m_ptCellCnt[i];
  }
  // fill an array with point/cell adjacency info
  m_ptAdjCells.assign(last, -1);
  for (size_t i = 0, c = 0; i < m_meshCells.size(); ++i, ++c)
  {
    for (int j = 0; j < m_meshCells[i + 1]; ++j)
    {
      int ptIdx = m_meshCells[i + 2 + j];
      int ix = m_ptCellIdx[ptIdx];
      while (m_ptAdjCells[ix] != -1)
        ix++;
      m_ptAdjCells[ix] = (int)c;
    }
    i += (m_meshCells[i + 1] + 1);
  }
} // MePolyPatcherImpl::CreateCellAdjacencyInfo
//------------------------------------------------------------------------------
/// \brief Gets the points that make up the cell and points from adjacent cells
/// \param[in] a_cellIdx index to the cell
/// \param[out] a_ptIdxs the point indexes for the points that make up this cell
/// \param[out] a_adjCellIdx the indexes to the adjacent cells
/// \param[out] a_adjPts the point indexes from adjacent cells
//------------------------------------------------------------------------------
void MePolyPatcherImpl::GetCellPtInfo(int a_cellIdx,
                                      VecInt& a_ptIdxs,
                                      VecInt& a_adjCellIdx,
                                      VecInt2d& a_adjPts)
{
  // get points for this cell
  a_ptIdxs.resize(0);
  int ix = m_cellIdx[a_cellIdx];
  for (int i = 0; i < m_meshCells[ix + 1]; ++i)
  {
    a_ptIdxs.push_back(m_meshCells[ix + 2 + i]);
  }
  // get adjacent cell points
  std::set<int> adjCells;
  adjCells.insert(a_cellIdx);
  a_adjPts.resize(0);
  a_adjCellIdx.resize(0);
  for (size_t i = 0; i < a_ptIdxs.size(); ++i)
  {
    const int& px = a_ptIdxs[i];
    ix = m_ptCellIdx[px];
    // adjacent cells
    for (int j = 0; j < m_ptCellCnt[px]; ++j)
    {
      const int& cx1 = m_ptAdjCells[ix + j]; // adjacent cell
      if (adjCells.find(cx1) != adjCells.end())
        continue;
      adjCells.insert(cx1);
      a_adjPts.push_back(VecInt());
      a_adjCellIdx.push_back(cx1);
      const int& cx = m_cellIdx[cx1];
      for (int k = 0; k < m_meshCells[cx + 1]; ++k)
      { // points on this cell
        const int& pt(m_meshCells[cx + 2 + k]);
        a_adjPts.back().push_back(pt);
      }
    }
  }
} // MePolyPatcherImpl::GetCellPtInfo

} // namespace xms

#ifdef CXX_TEST
////////////////////////////////////////////////////////////////////////////////

#include <xmsmesh/meshing/detail/MePolyPatcher.t.h>

#include <xmscore/testing/TestTools.h>

////////////////////////////////////////////////////////////////////////////////
/// \class MePolyPatcherUnitTests
/// \brief tester for the MePolyPaverToMeshPts class
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief test creating the class
//------------------------------------------------------------------------------
void MePolyPatcherUnitTests::testCreateClass()
{
  BSHP<xms::MePolyPatcher> p = xms::MePolyPatcher::New();
  TS_ASSERT(p);
} // MeIntersectPolysTest::testCreateClass
//------------------------------------------------------------------------------
/// \brief tests converting the polygon definition into the "side" definition
//------------------------------------------------------------------------------
void MePolyPatcherUnitTests::testPolyPtsToSides()
{
  xms::VecPt3d pts;
  {
    using namespace xms;
    pts = {{0, 0, 0},   {0, 10, 0},  {0, 20, 0},  {0, 30, 0}, {10, 30, 0}, {20, 30, 0},
           {30, 30, 0}, {30, 20, 0}, {30, 10, 0}, {30, 0, 0}, {20, 0, 0},  {10, 0, 0}};
  }
  xms::VecInt corner = {3, 6, 9};
  xms::MePolyPatcherImpl p;
  p.PolyPtsToSides(pts, corner);
  TS_ASSERT_EQUALS(4, p.m_pts.size());
  if (4 != p.m_pts.size())
    return;
  xms::VecPt3d2d basePts(4);
  {
    using namespace xms;
    basePts[0] = {{0, 30, 0}, {0, 20, 0}, {0, 10, 0}, {0, 0, 0}};
    basePts[1] = {{0, 0, 0}, {10, 0, 0}, {20, 0, 0}, {30, 0, 0}};
    basePts[2] = {{30, 30, 0}, {30, 20, 0}, {30, 10, 0}, {30, 0, 0}};
    basePts[3] = {{0, 30, 0}, {10, 30, 0}, {20, 30, 0}, {30, 30, 0}};
  }
  TS_ASSERT_DELTA_VECPT3D(basePts[0], p.m_pts[0], 1e-9);
  TS_ASSERT_DELTA_VECPT3D(basePts[1], p.m_pts[1], 1e-9);
  TS_ASSERT_DELTA_VECPT3D(basePts[2], p.m_pts[2], 1e-9);
  TS_ASSERT_DELTA_VECPT3D(basePts[3], p.m_pts[3], 1e-9);

  // change the starting point "[0]" of input polygon
  {
    using namespace xms;
    pts = {{0, 30, 0}, {10, 30, 0}, {20, 30, 0}, {30, 30, 0}, {30, 20, 0}, {30, 10, 0},
           {30, 0, 0}, {20, 0, 0},  {10, 0, 0},  {0, 0, 0},   {0, 10, 0},  {0, 20, 0}};
  }
  p.PolyPtsToSides(pts, corner);
  TS_ASSERT_EQUALS(4, p.m_pts.size());
  if (4 != p.m_pts.size())
    return;
  {
    using namespace xms;
    basePts[0] = {{0, 30, 0}, {0, 20, 0}, {0, 10, 0}, {0, 0, 0}};
    basePts[1] = {{0, 0, 0}, {10, 0, 0}, {20, 0, 0}, {30, 0, 0}};
    basePts[2] = {{30, 30, 0}, {30, 20, 0}, {30, 10, 0}, {30, 00, 0}};
    basePts[3] = {{0, 30, 0}, {10, 30, 0}, {20, 30, 0}, {30, 30, 0}};
  }
  TS_ASSERT_DELTA_VECPT3D(basePts[0], p.m_pts[0], 1e-9);
  TS_ASSERT_DELTA_VECPT3D(basePts[1], p.m_pts[1], 1e-9);
  TS_ASSERT_DELTA_VECPT3D(basePts[2], p.m_pts[2], 1e-9);
  TS_ASSERT_DELTA_VECPT3D(basePts[3], p.m_pts[3], 1e-9);
} // MePolyPatcherUnitTests::testPolyPtsToSides
//------------------------------------------------------------------------------
/// \brief tests converting the polygon definition into the "side" definition
//------------------------------------------------------------------------------
void MePolyPatcherUnitTests::testPolyPtsToSides1()
{
  xms::VecPt3d pts;
  {
    using namespace xms;
    pts = {{0, 0, 0},    {0, 70, 0},   {0, 100, 0}, {35, 100, 0}, {100, 100, 0},
           {100, 70, 0}, {100, 35, 0}, {100, 0, 0}, {70, 0, 0},   {35, 0, 0}};
  }
  xms::VecInt corner = {2, 4, 7};
  xms::MePolyPatcherImpl p;
  p.PolyPtsToSides(pts, corner);
  TS_ASSERT_EQUALS(4, p.m_pts.size());
  if (4 != p.m_pts.size())
    return;
  xms::VecPt3d2d basePts(4);
  {
    using namespace xms;
    basePts[0] = {{0, 100, 0}, {0, 70, 0}, {0, 0, 0}};
    basePts[1] = {{0, 0, 0}, {35, 0, 0}, {70, 0, 0}, {100, 0, 0}};
    basePts[2] = {{100, 100, 0}, {100, 70, 0}, {100, 35, 0}, {100, 0, 0}};
    basePts[3] = {{0, 100, 0}, {35, 100, 0}, {100, 100, 0}};
  }
  TS_ASSERT_DELTA_VECPT3D(basePts[0], p.m_pts[0], 1e-9);
  TS_ASSERT_DELTA_VECPT3D(basePts[1], p.m_pts[1], 1e-9);
  TS_ASSERT_DELTA_VECPT3D(basePts[2], p.m_pts[2], 1e-9);
  TS_ASSERT_DELTA_VECPT3D(basePts[3], p.m_pts[3], 1e-9);
} // MePolyPatcherUnitTests::testPolyPtsToSides1
//------------------------------------------------------------------------------
/// \brief test patch generation of a square with equal side spacing
//------------------------------------------------------------------------------
void MePolyPatcherUnitTests::testPatch00()
{
  xms::VecPt3d pts;
  {
    using namespace xms;
    pts = {{0, 0, 0},   {0, 10, 0},  {0, 20, 0},  {0, 30, 0}, {10, 30, 0}, {20, 30, 0},
           {30, 30, 0}, {30, 20, 0}, {30, 10, 0}, {30, 0, 0}, {20, 0, 0},  {10, 0, 0}};
  }
  xms::VecInt corner = {3, 6, 9};
  xms::MePolyPatcherImpl p;
  xms::VecPt3d mPts;
  xms::VecInt mCells;
  p.MeshIt(-1, pts, corner, 1e-9, mPts, mCells);
} // MePolyPatcherUnitTests::testPatch00
//------------------------------------------------------------------------------
/// \brief tests error detection in polygon for quad patch
//------------------------------------------------------------------------------
void MePolyPatcherUnitTests::testQuadPatchErrors()
{
  xms::VecPt3d pts;
  {
    using namespace xms;
    pts = {{1, 1, 0},   {0, 1, 0},   {0, 10, 0},  {0, 20, 0}, {0, 30, 0}, {10, 30, 0}, {20, 30, 0},
           {30, 30, 0}, {30, 20, 0}, {30, 10, 0}, {30, 0, 0}, {20, 0, 0}, {10, 0, 0},  {1, 0, 0}};
  }
  xms::VecInt corner = {4, 7, 10};
  xms::MePolyPatcherImpl p;
  TS_ASSERT_EQUALS(true, p.QuadPatchErrors(pts, corner));
  {
    TS_ASSERT_STACKED_ERRORS("---Non-convex corner detected in polygon. Aborting patch.\n\n");
  }

  {
    using namespace xms;
    pts = {{0, 0, 0},   {0, 10, 0},  {0, 20, 0},  {0, 29, 0}, {0, 30, 0}, {10, 30, 0}, {20, 30, 0},
           {30, 30, 0}, {30, 20, 0}, {30, 10, 0}, {30, 0, 0}, {20, 0, 0}, {10, 0, 0}};
  }
  corner = {3, 7, 10};
  TS_ASSERT_EQUALS(true, p.QuadPatchErrors(pts, corner));
  {
    TS_ASSERT_STACKED_ERRORS("---Non-convex corner detected in polygon. Aborting patch.\n\n");
  }
} // MePolyPatcherUnitTests::testQuadPatchErrors
//------------------------------------------------------------------------------
/// \brief tests for bug report 9226
//------------------------------------------------------------------------------
void MePolyPatcherUnitTests::testBug9226()
{
  xms::VecPt3d pts = {{-15, 38}, {17, 47}, {52, 50}, {95, 40},
                      {92, 27},  {51, 35}, {18, 33}, {-12, 26}};
  xms::VecInt corner = {3, 4, 7};
  xms::MePolyPatcherImpl p;
  xms::VecPt3d mpts;
  xms::VecInt mcells;
  TS_ASSERT_EQUALS(true, p.MeshIt(-1, pts, corner, 1e-9, mpts, mcells));
  xms::VecPt3d basePts = {{-15, 38}, {-12, 26}, {95, 40}, {92, 27},
                          {17, 47},  {18, 33},  {52, 50}, {51, 35}};
  TS_ASSERT_EQUALS_VEC(basePts, mpts);
  xms::VecInt baseCells = {9, 4, 0, 1, 5, 4, 9, 4, 4, 5, 7, 6, 9, 4, 6, 7, 3, 2};
  TS_ASSERT_EQUALS_VEC(baseCells, mcells);
} // MePolyPatcherUnitTests::testBug9226

#endif
