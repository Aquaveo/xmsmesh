//------------------------------------------------------------------------------
/// \file
/// \ingroup meshing
/// \copyright (C) Copyright Aquaveo 2018. Distributed under FreeBSD License
/// (See accompanying file LICENSE or https://aqaveo.com/bsd/license.txt)
//------------------------------------------------------------------------------

//----- Included files ---------------------------------------------------------

// 1. Precompiled header

// 2. My own header
#include <xmsmesh/meshing/detail/MeQuadBlossom.h>

// 3. Standard library headers
#include <cmath>
#include <numeric>

// 4. External library headers

// 5. Shared code headers
#include <xmscore/misc/XmError.h>
#include <xmscore/misc/xmstype.h>
#include <xmsgrid/ugrid/XmEdge.h>
#include <xmsgrid/ugrid/XmUGrid.h>
#include <xmsinterp/geometry/geoms.h>
#include <xmsmesh/meshing/detail/MeWeightMatcher.h>

// 6. Non-shared code headers

//----- Forward declarations ---------------------------------------------------

//----- External globals -------------------------------------------------------

//----- Namespace declaration --------------------------------------------------

namespace xms
{
//----- Constants / Enumerations -----------------------------------------------
const double PI_OVER_TWO = XM_PI / 2.0; ///< PI divided by 2
const double TWO_OVER_PI = 2.0 / XM_PI; ///< 2 divided by PI
const double TWO_PI = 2.0 * XM_PI;      ///< 2 times PI
const int DEFAULT_COST = -1000;         ///< The default weight for extra edges
//----- Classes / Structs ------------------------------------------------------

namespace
{
/// Edge information capturing points on edge, faces on left and right and
/// points to the left and right.
struct Edge
{
  /// Equals operator.
  bool operator==(const Edge& a_rhs) const
  {
    return m_p0 == a_rhs.m_p0 && m_p1 == a_rhs.m_p1 && m_fL == a_rhs.m_fL && m_fR == a_rhs.m_fR &&
           m_pL == a_rhs.m_pL && m_pR == a_rhs.m_pR;
  }

  int m_p0; ///< Point 0
  int m_p1; ///< Point 1
  int m_fL; ///< Face on the left looking from m_p0 to m_p1
  int m_fR; ///< Face on the right looking from m_p0 to m_p1
  int m_pL; ///< Point on the face (triangle) to the left opposite the edge
  int m_pR; ///< Point on the face (triangle) to the right opposite the edge
};

/// Vector of edge information.
typedef std::vector<Edge> VecEdge;

/// Identifies cell adjacent to a particular point p.
struct CellData
{
  /// Constructor.
  CellData(int a_priorPt, int a_nextPt, int a_cell)
  : m_priorPoint(a_priorPt)
  , m_nextPoint(a_nextPt)
  , m_cell(a_cell)
  {
  }

  // Equals operator.
  bool operator==(const CellData& a_rhs) const
  {
    return m_priorPoint == a_rhs.m_priorPoint && m_nextPoint == a_rhs.m_nextPoint &&
           m_cell == a_rhs.m_cell;
  }

  int m_priorPoint; ///< Index of the point coming into point p.
  int m_nextPoint;  ///< Index of the point going out from point p.
  int m_cell;       ///< Cell index.
};

typedef std::vector<CellData> VecCellData;             ///< Vector of CellData.
typedef std::vector<VecCellData> VecVecCellData;       ///< Vector of VecCellData.
typedef std::vector<VecVecCellData> VecVecVecCellData; ///< Vector of VecVecCellData.

/// Implementation of MeQuadBlossom.
class MeQuadBlossomImpl : public MeQuadBlossom
{
public:
  MeQuadBlossomImpl(const VecPt3d& a_points, const VecInt2d& a_triangles);

  virtual int PreMakeQuads() override;
  virtual BSHP<XmUGrid> MakeQuads(bool a_splitBoundaryPoints,
                                  bool a_useAngle) override;

  virtual BSHP<XmUGrid> _MakeQuads(bool a_splitBoundaryPoints = true,
                                   bool a_useAngle = false,
                                   int cost = -10);
  const VecPt3d& GetLocations() const;

  void PrepareForMatch(int a_cost,
                       bool a_useAngle);
  VecInt MatchTriangles(bool a_splitBoundaryPoints);
  void GetEdges(int a_cost);
  void ProcessPointChains(int a_p,
                          VecVecCellData& a_chains,
                          int extra_cost,
                          VecEdge& a_interiorEdges,
                          VecEdge& a_boundaryEdges,
                          VecMeEdge& a_extraEdges,
                          VecInt2d& a_extraPoints);
  void ExtractCellsAdjacentToPoint(VecVecCellData& a_pointFaceChains,
                                   VecVecVecCellData& a_sortedFaceChains);
  int FindFaceThatStartsWith(const VecCellData& facesData, int a_p);
  int FindChainThatStartsWith(const VecVecCellData& a_chains, int a_p);
  VecVecCellData SortIntoPointChains(VecCellData& faces);
  VecMeEdge CalculateEdgeCost(bool a_useAngle);
  VecInt2d EliminateEdges(const VecInt& eliminate);
  VecEdge GetInteriorEdges(int a_p, const VecCellData& a_chain);
  Edge GetBoundaryEdges(int a_p,
                        const VecCellData& a_chain,
                        int cost,
                        VecMeEdge& a_extraEdges,
                        VecInt2d& a_extraPoints);

//private: // all public for testing
  VecPt3d m_points;        ///< points of the triangular mesh
  VecInt2d m_faces;        ///< initially triangles becomes quads
  VecPt3d m_splitPoints;   ///< The potential new points for split points.
  VecEdge m_interiorEdges; ///< The interior edges of the triangular mesh.
  VecEdge m_boundaryEdges; ///< The boundary edges of the triangular mesh.
  VecMeEdge m_costs;       ///< The costs of each edge between designated faces.
  /// The adjacent boundary triangles that can be split. { face0, face1, cost }
  VecMeEdge m_extraEdges;
  /// For every extra edge the boundary point and prior and next point of each
  /// of the two extra edge faces.
  /// { p, face0_prior, face0_next, face1_prior, face1_next }
  VecInt2d m_extraPoints;
}; // class MeQuadBlossomImpl

/// An adjacent point and midpoint pair.
typedef std::pair<int, int> AdjPointMidpoint;
/// A vector of adjacent point and midpoints for a given point where the
/// adjacent point's index (AdjPointMidpoints.first) is greater than the
/// index of the base point.
typedef std::vector<AdjPointMidpoint> VecAdjPointMidpoint;
/// A vector of adjacent point midpoint vectors. Indexed point UGrid point
/// indices.
typedef std::vector<VecAdjPointMidpoint> VecVecAdjPointMidpoint;

//----- Internal functions -----------------------------------------------------

//----- Class / Function definitions -------------------------------------------

//------------------------------------------------------------------------------
/// \brief Converts the 3 edge 2D cells of a UGrid into a 2D vector.
/// \param[in] a_ugrid The UGrid to convert.
/// \return A vector of triangle indices for each cell.
//------------------------------------------------------------------------------
VecInt2d UGrid2dToVecInt2d(BSHP<XmUGrid> a_ugrid)
{
  int numCells = a_ugrid->GetCellCount();
  VecInt2d cells(numCells, VecInt());
  for (int cellIdx = 0; cellIdx < numCells; ++cellIdx)
  {
    VecInt& cell = cells[cellIdx];
    if (a_ugrid->GetCellDimension(cellIdx) == 2 && a_ugrid->GetCellEdgeCount(cellIdx) == 3)
    {
      VecInt cellPoints = a_ugrid->GetCellPoints(cellIdx);
      for (auto pointIdx : cellPoints)
      {
        cell.push_back(pointIdx);
      }
    }
  }

  return cells;
} // UGrid2dToVecInt2d
//------------------------------------------------------------------------------
/// \brief Build a UGrid based on quads and triangles given a vector of points
/// and a vector of cell triangles or quads.
/// \param[in] a_points The new UGrid's points.
/// \param[in] a_faces A vector of vectors containing 3 or 4 point indices.
/// \return The UGrid generated from the points and faces.
//------------------------------------------------------------------------------
BSHP<XmUGrid> BuildUGrid(const VecPt3d& a_points, const VecInt2d& a_faces)
{
  VecInt cells;
  for (auto& face : a_faces)
  {
    if (face.size() == 4)
    {
      cells.push_back(XMU_QUAD);
      cells.push_back(4);
      cells.insert(cells.end(), face.begin(), face.end());
    }
    else if (face.size() == 3)
    {
      cells.push_back(XMU_TRIANGLE);
      cells.push_back(3);
      cells.insert(cells.end(), face.begin(), face.end());
    }
    else
    {
      XM_ASSERT(0);
    }
  }

  return XmUGrid::New(a_points, cells);
} // BuildUGrid
//------------------------------------------------------------------------------
/// \brief Get the angle between 3 points.
/// \param[in] p0 The first end point.
/// \param[in] pc The center point.
/// \param[in] p1 The second end point.
/// \return The angle in radians (range of 0 -> 2 PI).
//------------------------------------------------------------------------------
double GetAngle(const Pt3d& p0, const Pt3d& pc, const Pt3d& p1)
{
  Pt3d d0c = gmCreateVector(p0, pc);
  Pt3d d1c = gmCreateVector(p1, pc);
  double angle_c0 = atan2(d0c[1], d0c[0]);
  double angle_c1 = atan2(d1c[1], d1c[0]);
  double angle_0c1 = angle_c0 - angle_c1;
  if (angle_0c1 < 0)
  {
    angle_0c1 += TWO_PI;
  }
  else if (angle_0c1 > TWO_PI)
  {
    angle_0c1 -= TWO_PI;
  }

  return angle_0c1;
} // GetAngle
//------------------------------------------------------------------------------
/// \brief Get the cost based on angle of two adjacent triangles {p0,p1,p2} and
/// {p3,p2,p1}.
/// \param[in] p0 The point opposite the shared edge of the triangle on the left.
/// \param[in] p1 The first point of the shared edge between the triangles.
/// \param[in] p2 The second point of the shared edge between the triangles.
/// \param[in] p3 The point opposite the shared edge of the triangle on the right.
/// \return A cost or weight integer value between 0 and 1,000.  If the two
/// triangles form a square, the cost will be the maximum (1,000).  If they
/// combine to form a more skewed quadrilateral, the cost tends toward 0.
//------------------------------------------------------------------------------
int GetEtaAngle(const Pt3d& p0, const Pt3d& p1, const Pt3d& p2, const Pt3d& p3)
{
  double a0 = GetAngle(p1, p0, p2);
  double a1 = GetAngle(p3, p1, p0);
  double a2 = GetAngle(p0, p2, p3);
  double a3 = GetAngle(p2, p3, p1);
  if (a0 > XM_PI)
  {
    a0 = a0 - XM_PI;
  }
  if (a1 > XM_PI)
  {
    a1 = a1 - XM_PI;
  }
  if (a2 > XM_PI)
  {
    a2 = a2 - XM_PI;
  }
  if (a3 > XM_PI)
  {
    a3 = a3 - XM_PI;
  }
  double etaMax1 = std::max(fabs(PI_OVER_TWO - a0), fabs(PI_OVER_TWO - a1));
  double etaMax2 = std::max(fabs(PI_OVER_TWO - a2), fabs(PI_OVER_TWO - a3));
  double etaMax = 1.0 - TWO_OVER_PI * std::max(etaMax1, etaMax2);
  double eta = std::max(0.0, etaMax);
  return int(1000.0 * eta + 0.5);
} // GetEtaAngle
//------------------------------------------------------------------------------
/// \brief Get the cost based on distance of two adjacent triangles {p0,p1,p2}
/// and {p3,p2,p1}.
/// \param[in] p0 The point opposite the shared edge of the triangle on the left.
/// \param[in] p1 The first point of the shared edge between the triangles.
/// \param[in] p2 The second point of the shared edge between the triangles.
/// \param[in] p3 The point opposite the shared edge of the triangle on the right.
/// \return A cost or weight integer value between 0 and 1,000.  If the two
/// triangles form a square, the cost will be the maximum (1,000).  If they
/// combine to form a more skewed quadrilateral, the cost tends toward 0.
//------------------------------------------------------------------------------
int GetEtaDistance(const Pt3d& p0, const Pt3d& p1, const Pt3d& p2, const Pt3d& p3)
{
  double a = gmXyDistanceSquared(p1, p0);
  double b = gmXyDistanceSquared(p0, p2);
  double c = gmXyDistanceSquared(p1, p3);
  double d = gmXyDistanceSquared(p3, p2);
  double e = gmXyDistanceSquared(p1, p2);
  double f = gmXyDistanceSquared(p0, p3);
  //XM_ASSERT(a + b >= e);
  //XM_ASSERT(c + d >= e);
  //XM_ASSERT(a + c >= f);
  //XM_ASSERT(b + d >= f);
  double inveta1 = std::max((a + b) / e, e / (a + b));
  double inveta2 = std::max((c + d) / e, e / (c + d));
  double inveta3 = std::max((a + c) / f, f / (a + c));
  double inveta4 = std::max((b + d) / f, f / (b + d));
  double inveta = std::max(inveta1 + inveta2, inveta3 + inveta4);
  //XM_ASSERT(inveta - (a + b + c + d)/std::min(e, f) < 1e-5);
  return int(1000.0 * 2.0 / inveta + 0.5);
} // GetEtaDistance
//------------------------------------------------------------------------------
/// \brief Gets the opposite point for each edge adjacent to a given point.
/// \param[in] a_ugrid The UGrid.
/// \param[in] a_pointIdx The point to get adjacent points from.
/// \param[out] a_edgePoints The adjacent points.
//------------------------------------------------------------------------------
void GetPointIdxsAttachedByEdge(BSHP<XmUGrid> a_ugrid, int a_pointIdx, VecInt& a_edgePoints)
{
  // TODO: Use new UGrid code in xmsgrid
  a_edgePoints.clear();
  VecInt associatedCells = a_ugrid->GetPointAdjacentCells(a_pointIdx);
  if (associatedCells.size() == 0)
  {
    return;
  }
  for (int i = 0; i < associatedCells.size(); ++i)
  {
    for (int j = 0; j < a_ugrid->GetCellEdgeCount(associatedCells[i]); ++j)
    {
      XmEdge temp = a_ugrid->GetCellEdge(associatedCells[i], j);
      if (temp.GetFirst() == a_pointIdx)
      {
        a_edgePoints.push_back(temp.GetSecond());
      }
      else if (temp.GetSecond() == a_pointIdx)
      {
        a_edgePoints.push_back(temp.GetFirst());
      }
    }
  }
  std::sort(a_edgePoints.begin(), a_edgePoints.end());
  auto it = std::unique(a_edgePoints.begin(), a_edgePoints.end());
  a_edgePoints.erase(it, a_edgePoints.end());
} // XmUGridImpl::GetPointIdxsAttachedByEdge
//------------------------------------------------------------------------------
/// \brief Finds the pre-calculated midpoint index a given a pair a UGrid
/// points that form an edge.
/// \param[in] a_edgeMidsides A vector indexed by the lesser of the point
/// indices for an edge that gives a precalculated vector of pairs for edges
/// emanating from that point.  The first element of each pair is the adjacent
/// point index on that edge (which is greater than the first point index).
/// The second element of the pair is the index of the midpoint for that edge.
/// \param[in] a_idx0 The first edge point.
/// \param[in] a_idx1 The second edge point.
/// \return the index of the midpoint.
//------------------------------------------------------------------------------
int FindMidpoint(const VecVecAdjPointMidpoint& a_edgeMidsides, int a_idx0, int a_idx1)
{
  if (a_idx0 > a_idx1)
    std::swap(a_idx0, a_idx1);

  auto& adjPoints = a_edgeMidsides[a_idx0];
  for (auto& adjPoint : adjPoints)
  {
    if (adjPoint.first == a_idx1)
    {
      return adjPoint.second;
    }
  }

  XM_ASSERT(0);
  return XM_NONE;
} // FindMidpoint

//------------------------------------------------------------------------------
/// \brief Constructor
//------------------------------------------------------------------------------
MeQuadBlossomImpl::MeQuadBlossomImpl(const VecPt3d& a_points, const VecInt2d& a_triangles)
: m_points(a_points)
, m_faces(a_triangles)
{
} // MeQuadBlossomImpl::MeQuadBlossomImpl
//------------------------------------------------------------------------------
/// \brief Estimate time to make quads and return whether triangles will be
/// expected.
/// \param[out] a_timeEstimateInMinutes The estimated time to generate quads
/// in minutes.
/// \return The number of boundary edges. If odd then at least one triangles
/// will be expected in the resulting mesh.
//------------------------------------------------------------------------------
int MeQuadBlossomImpl::PreMakeQuads()
{
  if (m_boundaryEdges.empty())
  {
    GetEdges(DEFAULT_COST);
  }
  return (int)m_boundaryEdges.size();
} // MeQuadBlossomImpl::PreMakeQuads
//------------------------------------------------------------------------------
/// \brief Turn faces from triangles into quads by using MeWeightMatcher
/// to identify edges to drop and then by calling EliminateEdges to remove them.
/// \param[in] a_splitBoundaryPoints If necessary, split boundary points to
/// create quads from "pseudo" edges between unmatched boundary triangles
/// separated by at least one other triangle.
/// \param[in] a_useAngle If true use GetEtaAngle else use GetEtaDistance. Used
/// to compute the interior edge cost.
/// \return An XmUGrid with the quads and any new points.
//------------------------------------------------------------------------------
BSHP<XmUGrid> MeQuadBlossomImpl::MakeQuads(bool a_splitBoundaryPoints,
                                           bool a_useAngle)
{
  return _MakeQuads(a_splitBoundaryPoints, a_useAngle);
} // MeQuadBlossomImpl::MakeQuads
//------------------------------------------------------------------------------
/// \brief Turn faces from triangles into quads by using MeWeightMatcher
/// to identify edges to drop and then by calling EliminateEdges to remove them.
/// \param[in] a_splitBoundaryPoints If necessary, split boundary points to
/// create quads from "pseudo" edges between unmatched boundary triangles
/// separated by at least one other triangle.
/// \param[in] a_cost The cost to associate with extra boundary edges.
/// \param[in] a_useAngle If true use GetEtaAngle else use GetEtaDistance. Used
/// to compute the interior edge cost.
/// \return An XmUGrid with the quads and any new points.
//------------------------------------------------------------------------------
BSHP<XmUGrid> MeQuadBlossomImpl::_MakeQuads(bool a_splitBoundaryPoints /*= true */,
                                            bool a_useAngle /*= false*/,
                                            int a_cost /*= -10*/)
{
  PrepareForMatch(a_cost, a_useAngle);
  VecInt eliminate = MatchTriangles(a_splitBoundaryPoints);
  EliminateEdges(eliminate);
  BSHP<XmUGrid> ugrid = BuildUGrid(m_points, m_faces);
  return ugrid;
} // MeQuadBlossomImpl::MakeQuads
//------------------------------------------------------------------------------
/// \brief Prepare for MakeQuads call by computing the edges and costs.
///
/// Calculates an array of costs for each edge (or pair of adjacent triangles).
/// The first elements of a_cost correspond to the interior edge returned. The
/// remaining elements will each identify two boundary triangles that share a
/// point but not an edge (i.e., there is at least one triangle between them
/// and the cost of treating that point like a pseudo edge that can be split to
/// form an edge between the triangles and turn them into quads.
/// \param[in] a_cost The cost to associate with extra pseudo edges.
/// \param[in] a_useAngle If true, use the angle method of computing edge costs;
/// otherwise, use the distance method.
//------------------------------------------------------------------------------
void MeQuadBlossomImpl::PrepareForMatch(int a_cost,
                                        bool a_useAngle)
{
  if (m_boundaryEdges.empty())
  {
    GetEdges(a_cost);
  }
  else if (a_cost != DEFAULT_COST)
  {
    for (auto& edge : m_extraEdges)
    {
      edge.m_weight = a_cost;
    }
  }
  m_costs = CalculateEdgeCost(a_useAngle);
  m_costs.insert(m_costs.end(), m_extraEdges.begin(), m_extraEdges.end());
} // MeQuadBlossomImpl::PrepareForMatch
//------------------------------------------------------------------------------
/// \brief Determine which edges to eliminate and which boundary points to split
/// to turn adjacent triangles into quads.
/// \param[in] a_splitBoundaryPoints if true, all matching the point between
/// adjacent boundary triangles that do not share an edge as if it were an edge.
/// If true you will get at most one triangle. You will get no triangles if
/// given an even number of boundary edges.
/// \return The edges to eliminate and boundary points to split to turn adjacent
/// triangles into quads. If elimnate[i] > the number of interiorEdges, then it
/// identifies a boundary point to be split.
//------------------------------------------------------------------------------
VecInt MeQuadBlossomImpl::MatchTriangles(bool a_splitBoundaryPoints)
{
  bool maxCardinality = a_splitBoundaryPoints;
  BSHP<MeWeightMatcher> matcher = MeWeightMatcher::New();
  VecInt result = matcher->MatchWeights(m_costs, maxCardinality);

  // VecInt unmatched;
  VecInt eliminate;
  for (int i = 0; i < (int)result.size(); ++i)
  {
    int j = result[i];
    // if (j != -1)
    //{
    // unmatched.push_back(i);
    //}
    // if j == -1, then j < i
    if (j > i)
    {
      for (int k = 0; k < m_costs.size(); ++k)
      {
        if ((i == m_costs[k].m_f0 && j == m_costs[k].m_f1) ||
            (i == m_costs[k].m_f1 && j == m_costs[k].m_f0))
        {
          eliminate.push_back(k);
        }
      }
    }
  }
  return eliminate;
} // MeQuadBlossomImpl::MatchTriangles
//------------------------------------------------------------------------------
/// \brief Returns interior, boundary, and extra edges.
/// \param[in] a_cost The cost to associate with psuedo-extra edges that are
/// conceptually between two boundary edges that share a point that has more
/// than 3 total edges.
//------------------------------------------------------------------------------
void MeQuadBlossomImpl::GetEdges(int a_cost)
{
  VecVecCellData dummy;
  VecVecVecCellData sorted;
  ExtractCellsAdjacentToPoint(dummy, sorted);

  int p = 0;
  m_interiorEdges.clear();
  m_boundaryEdges.clear();
  m_extraEdges.clear();
  m_extraPoints.clear();
  m_splitPoints.clear();
  for (auto& chains : sorted)
  {
    VecEdge interior, boundary;
    ProcessPointChains(p, chains, a_cost, interior, boundary, m_extraEdges, m_extraPoints);
    m_interiorEdges.insert(m_interiorEdges.begin(), interior.begin(), interior.end());
    m_boundaryEdges.insert(m_boundaryEdges.begin(), boundary.begin(), boundary.end());
    p += 1;
  }
} // MeQuadBlossomImpl::GetEdges
//------------------------------------------------------------------------------
/// \brief Process a chain of sorted groups of face data for a point, producing
/// interior, boundary, and extra edges.
/// \param[in] a_p The point edge index.
/// \param[in] a_chains A vector of vectors of sorted FaceData.
/// \param[in] extra_cost Cost to associate with extra edges.
/// \param[out] a_interiorEdges The interior edges.
/// \param[out] a_boundaryEdges The boundary edges (edge.fR and edge.pR are
/// XM_NONE).
/// \param[out] a_extraEdge An array of [face0, face1, and a cost]. Face0 and
/// Face1 are boundary faces that share a point with more than 3 edges.
//------------------------------------------------------------------------------
void MeQuadBlossomImpl::ProcessPointChains(int a_p,
                                           VecVecCellData& a_chains,
                                           int extra_cost,
                                           VecEdge& a_interiorEdges,
                                           VecEdge& a_boundaryEdges,
                                           VecMeEdge& a_extraEdges,
                                           VecInt2d& a_extraPoints)
{
  a_interiorEdges.clear();
  a_boundaryEdges.clear();

  if (a_chains.size() == 0)
  {
    return;
  }

  if (a_chains.size() == 1)
  {
    auto& chain = a_chains[0];
    auto& first = chain[0];
    auto& last = chain.back();
    if (last.m_nextPoint == first.m_priorPoint) // This is an interior point
    {
      chain.push_back(first);
      a_interiorEdges = GetInteriorEdges(a_p, chain);
      return;
    }
    a_interiorEdges = GetInteriorEdges(a_p, chain);
    Edge boundaryEdge = GetBoundaryEdges(a_p, chain, extra_cost, a_extraEdges, a_extraPoints);
    a_boundaryEdges.push_back(boundaryEdge);
    return;
  }
  // There is  more than one chain, so we have possibly several groups of faces sharing the point
  for (auto& chain : a_chains)
  {
    VecEdge edges = GetInteriorEdges(a_p, chain);
    a_interiorEdges.insert(a_interiorEdges.begin(), edges.begin(), edges.end());
    auto boundaryEdge = GetBoundaryEdges(a_p, chain, extra_cost, a_extraEdges, a_extraPoints);
    a_boundaryEdges.push_back(boundaryEdge);
  }
} // MeQuadBlossomImpl::ProcessVertexChains
//------------------------------------------------------------------------------
/// \brief Extract a sorted vector of triangles attached to each point.
/// \param[out]: a_pointFaceChains a vector of vectors of face data, but
/// unsorted.
/// \param[out]: a_sortedFaceChains vector of vector of vectors of face data,
/// where the faces adjacent to each point have be sorted such that the
/// element[i].m_nextPoint == elements[i+1].m_priorPoint.  There are separate
/// vectors if there is a break in that pattern.  If there are multiple chains
/// for any point, then that means there are multiple triangles separated by
/// boundary edges.  If there is only one chain for a point, and the last
/// element's m_nextPoint == first element's m_priorPoint, then this is the
/// chain for an interior point.
//------------------------------------------------------------------------------
void MeQuadBlossomImpl::ExtractCellsAdjacentToPoint(VecVecCellData& a_pointFaceChains,
                                                    VecVecVecCellData& a_sortedFaceChains)
{
  a_pointFaceChains.clear();
  a_pointFaceChains.resize(m_points.size(), VecCellData());
  for (int fi = 0; fi < (int)m_faces.size(); ++fi)
  {
    VecInt face = m_faces[fi];

    int p0 = face[0];
    int p1 = face[1];
    int p2 = face[2];
    a_pointFaceChains[p0].push_back(CellData(p2, p1, fi));
    a_pointFaceChains[p1].push_back(CellData(p0, p2, fi));
    a_pointFaceChains[p2].push_back(CellData(p1, p0, fi));
  }

  a_sortedFaceChains.clear();
  for (auto& pointFace : a_pointFaceChains)
  {
    VecVecCellData chains = SortIntoPointChains(pointFace);
    a_sortedFaceChains.push_back(chains);
  }
} // MeQuadBlossomImpl::ExtractCellsAdjacentToPoint
//------------------------------------------------------------------------------
/// \brief Finds a face in faceData whose m_priorPoint == a_p
/// \param[in] facesData A vector of face data for faces incident to point a_p.
/// \param[in] a_p A point index in the mesh.
/// \return Index of the element in faces that whose m_priorPoint == a_p.
//------------------------------------------------------------------------------
int MeQuadBlossomImpl::FindFaceThatStartsWith(const VecCellData& facesData, int a_p)
{
  for (size_t i = 0; i < facesData.size(); ++i)
  {
    if (a_p == facesData[i].m_priorPoint)
    {
      return (int)i;
    }
  }
  return XM_NONE;
} // MeQuadBlossomImpl::FindFaceThatStartsWith
//------------------------------------------------------------------------------
/// \brief Finds a vector of FaceData whose first element's m_priorPoint == a_p.
///
/// Note: After searching forward for a elements to find a chain, the last
/// element might chain together with another chain already formed.
///
/// \param[in] a_chains A chain of face data elements.
/// \param[in] a_p Index of the point we want to match.
/// \return Index of the matching chain in chains.
//------------------------------------------------------------------------------
int MeQuadBlossomImpl::FindChainThatStartsWith(const VecVecCellData& a_chains, int a_p)
{
  for (size_t i = 0; i < a_chains.size(); ++i)
  {
    if (a_p == a_chains[i][0].m_priorPoint)
    {
      return (int)i;
    }
  }
  return XM_NONE;
} // MeQuadBlossomImpl::FindChainThatStartsWith
//------------------------------------------------------------------------------
/// \brief Helper function to ExtractCellsAdjacentToPoint to sort the faces
/// incident to the point into adjacent lists.
/// \param[in] faces Faces incident to a point to be sorted into adjacency
/// lists.
/// \return A list of chains, where each chain is a list of face data where
/// chain[i][1] == chain[i+1][0].  In most cases, chains will have just one
/// chain, and the number of elements in each chain will be small. (In most
/// meshes, there are just 3-6 triangles incident on any one point.)
//------------------------------------------------------------------------------
VecVecCellData MeQuadBlossomImpl::SortIntoPointChains(VecCellData& facesData)
{
  VecVecCellData chains;
  VecCellData remaining = facesData;
  int face_i = 0;
  VecCellData chain;
  while (remaining.size() > 0)
  {
    auto iter = remaining.begin() + face_i;
    CellData face = remaining[face_i];
    remaining.erase(iter);
    int p = face.m_nextPoint;
    chain.push_back(face);
    face_i = FindFaceThatStartsWith(remaining, p);
    if (face_i == XM_NONE)
    {
      face_i = 0;
      // Now search chains for a chain that begins with p
      int chain_i = FindChainThatStartsWith(chains, p);
      if (chain_i != XM_NONE)
      {
        // insert chain at the front of chains[chain_i]
        chains[chain_i].insert(chains[chain_i].begin(), chain.begin(), chain.end());
      }
      else
      {
        chains.push_back(chain);
      }
      chain.clear();
    }
  }
  return chains;
} // MeQuadBlossomImpl::SortIntoVertexChains
//------------------------------------------------------------------------------
/// \brief Compute a cost associated with each edge in m_interiorEdges for use
/// in converting faces to quads.
/// \param[in] a_useAngle If true use GetEtaAngle else use GetEtaDistance.
/// \return An array of arrays for each edge of the form [left_face_index,
/// right_face_index, cost].
//------------------------------------------------------------------------------
VecMeEdge MeQuadBlossomImpl::CalculateEdgeCost(bool a_useAngle)
{
  auto& eta = a_useAngle ? GetEtaAngle : GetEtaDistance;
  VecMeEdge costs;
  // Note: m_interiorEdges contains only one of each edge and its twin, not both
  for (auto& edge : m_interiorEdges)
  {
    // p1 and p2 are the points of the edge
    const Pt3d& p1 = m_points[edge.m_p0];
    const Pt3d& p2 = m_points[edge.m_p1];

    // p0 and p3 are the opposite points from the edge in the respective faces of the edge
    // If those faces are faces, they are the triangle points not on the edge.
    const Pt3d& p0 = m_points[edge.m_pL];
    const Pt3d& p3 = m_points[edge.m_pR];

    int cost = int(eta(p0, p1, p2, p3));
    costs.push_back({edge.m_fL, edge.m_fR, cost});
  }
  return costs;
} // MeQuadBlossomImpl::CalculateEdgeCost
//------------------------------------------------------------------------------
/// \brief Eliminate the interior edges indexed in eliminate to turn self.face
/// from triangles to quad.
///
/// Primary helper to MakeQuads.  It makes compacts faces by dropping edge
/// between triangles as indicated by the results of maxWeightMatching().  For
/// each edge identified in 'eliminate', replace the triangle face in fL with a
/// quad with points [p0, pL, p1, pR] and the triangle face as fR with None.
/// Then compact self.faces to remove the empty elements.
///
/// \param[in] eliminate Array of indices into interiorEdges of edges that are
/// to be dopped. If value is -1 then it's not dropped. Otherwise, the edge
/// between point i and point eliminate[i] is to be dropped.
/// \return Array of resulting faces (also the same as member faces).
//------------------------------------------------------------------------------
VecInt2d MeQuadBlossomImpl::EliminateEdges(const VecInt& eliminate)
{
  int nInteriorEdges = (int)m_interiorEdges.size();

  VecInt splitPointsRedir(m_points.size(), 0);
  std::iota(splitPointsRedir.begin(), splitPointsRedir.end(), 0);
  for (auto i : eliminate)
  {
    if (i >= nInteriorEdges)
    {
      // split the point
      //
      // if this splits, it will produce quads:
      //   for first.m_cell: { first.m_nextPoint, first.m_priorPoint, p, px }
      //   for last.m_cell: { p, last.m_nextPoint, last.m_priorPoint, px }
      //
      //          first.m_nextPoint  \   px  / last.m_priorPoint
      //                              \     /
      //                 first.m_cell  \   /  last.m_cell
      //                                \ /
      //   first.m_priorPoint ---------- p -------- last.m_nextPoint
      //
      int j = i - nInteriorEdges;
      const auto& extra = m_extraPoints[j];
      const auto& edge = m_costs[i];
      int p = extra[0];
      int firstVertexOut = extra[1];
      int firstVertexIn = extra[2];
      int lastVertexOut = extra[3];
      int lastVertexIn = extra[4];
      int first = edge.m_f0;
      int last = edge.m_f1;
      int px = (int)m_points.size();
      m_points.push_back(m_splitPoints[j]);
      m_faces[first] = { firstVertexOut, firstVertexIn, p, px };
      m_faces[last] = { p, lastVertexOut, lastVertexIn, px };

      splitPointsRedir[p] = px;
    }
  }

  for (auto i : eliminate)
  {
    if (i < nInteriorEdges)
    {
      // eliminate the edge
      const Edge& edge = m_interiorEdges[i];
      int p0 = splitPointsRedir[edge.m_p0];
      int pR = splitPointsRedir[edge.m_pR];
      int p1 = splitPointsRedir[edge.m_p1];
      int pL = splitPointsRedir[edge.m_pL];
      m_faces[edge.m_fL] = {p0, pR, p1, pL};
      m_faces[edge.m_fR].clear();
    }
  }

  VecInt2d tempFaces;
  for (auto& face : m_faces)
  {
    if (!face.empty())
    {
      tempFaces.push_back(face);
    }
  }

  m_faces = tempFaces;
  return tempFaces;
} // MeQuadBlossomImpl::EliminateEdges
//------------------------------------------------------------------------------
/// \brief Return a list of interior edges for the edge between every face pair
/// in the chain if p < p0 of the 2nd.
/// \param[in] a_p The point this chain is adjacent to.  (Will be the leading
/// point of every edge created.)
/// \param[in] chain A sorted list of adjacent faces incident to p.
/// \return A vector of the interior edges.
//------------------------------------------------------------------------------
VecEdge MeQuadBlossomImpl::GetInteriorEdges(int a_p, const VecCellData& chain)
{
  VecEdge edges;
  XM_ASSERT(chain.size() > 0);
  for (size_t i = 0; i < chain.size() - 1; ++i)
  {
    const CellData& last = chain[i];
    const CellData& face = chain[i + 1];
    int p1 = face.m_priorPoint;
    // only create an edge if p1 > a_p.  Let the chain of p1 create it otherwise.
    if (p1 > a_p)
    {
      // last: [m_priorPoint=2, m_nextPoint=5, face=3],
      // current: [m_priorPoint=5, m_nextPoint=8, face=6],
      // a_p: 4=> [p0=4, p1=5, fL=3, fR=6, pL=2, pR=8]

      // last: [m_priorPoint=3, m_nextPoint=0, face=0],
      // current: [m_priorPoint=0, m_nextPoint=1, face=1], a_p: 4
      //  => [p0=4, p1=0, fL=0, fR=1, pL=3, pR=1]
      edges.push_back({a_p, p1, last.m_cell, face.m_cell, last.m_priorPoint, face.m_nextPoint});
    }
  }
  return edges;
} // MeQuadBlossomImpl::GetInteriorEdges
//------------------------------------------------------------------------------
/// \brief Return a boundary edge and optional extra_edge for point p from
/// chain.
/// \param[in] a_p The point this chain is adjacent to.
/// \param[in] chain A list of faces of the form [priorPt, nextPt, face]
/// \param[in] cost The cost to associated with an extra edge.
/// \param[in/out] a_extraEdges An extra edge of the form [f0, f1, cost] if
/// there are more than 2 faces in the chain.  f0 and f1 are boundary triangle
/// incident to a_p.
/// \return A vector of the boundary edges.
//------------------------------------------------------------------------------
Edge MeQuadBlossomImpl::GetBoundaryEdges(int a_p,
                                         const VecCellData& chain,
                                         int cost,
                                         VecMeEdge& a_extraEdges,
                                         VecInt2d& a_extraPoints)
{
  // a_p = 5
  const CellData& last = chain.back(); // [4, 3, 0] => [5, 3, 0, None, 4, None]
  Edge edge = {a_p, last.m_nextPoint, last.m_cell, XM_NONE, last.m_priorPoint, XM_NONE};

  if (chain.size() > 2)
  {
    const CellData& first = chain[0]; // [m_priorPoint=1, m_nextPoint=7, face=6]
    a_extraEdges.push_back({first.m_cell, last.m_cell, cost}); // [6, 0, cost]
    // if this splits, it will produce quads:
    //   for first.m_cell: { a_p, first.m_nextPoint, first.m_priorPoint, p' }
    //   for last.m_cell: { last.m_nextPoint, last.m_priorPoint, a_p, p' }
    //
    //          last.m_nextP \        / first.m_priorP
    //                        \  p'  /
    //           last.m_cell   \    /   first.m_cell
    //  last.m_priorP ---------- a_p -------- first.m_nextP
    a_extraPoints.push_back(
      {a_p, first.m_nextPoint, first.m_priorPoint, last.m_nextPoint, last.m_priorPoint});
    // compute p'
    Pt3d px;
    int nTris = (int)chain.size();
    const CellData& tri = chain[nTris >> 1];
    const Pt3d& p0 = m_points[a_p];
    const Pt3d& p1 = m_points[tri.m_priorPoint];
    if (nTris & 0x1)
    {
      // There are an odd number of faces, so use the centroid of the
      // middle triangle.
      const Pt3d& p2 = m_points[tri.m_nextPoint];
      px = (p0 + p1 + p2) / 3.0;
    }
    else
    {
      // There are an even number of faces, so use the middle triangle's
      // incoming edge.
      px = (p0 + p1) / 2.0;
    }
    m_splitPoints.push_back(px);
  }
  return edge;

} // MeQuadBlossomImpl::GetBoundaryEdges

} // namespace

////////////////////////////////////////////////////////////////////////////////
/// \class MeQuadBlossom
/// \brief Class to convert 2D grid of triangles to quads.
/// \see MeQuadBlossomImpl
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief Create new MeQuadBlossom.
/// \param[in] a_ugrid The UGrid to be converted from triangles to quads.
/// \return The new MeQuadBlossom.
//------------------------------------------------------------------------------
BSHP<MeQuadBlossom> MeQuadBlossom::New(BSHP<XmUGrid> a_ugrid)
{
  const VecPt3d& points = a_ugrid->GetLocations();
  VecInt2d triangles = UGrid2dToVecInt2d(a_ugrid);
  BSHP<MeQuadBlossom> wm(new MeQuadBlossomImpl(points, triangles));
  return wm;
} // MeQuadBlossom::New
//------------------------------------------------------------------------------
/// \brief Constructor.
//------------------------------------------------------------------------------
MeQuadBlossom::MeQuadBlossom()
{
} // MeQuadBlossom::MeQuadBlossom
//------------------------------------------------------------------------------
/// \brief Destructor.
//------------------------------------------------------------------------------
MeQuadBlossom::~MeQuadBlossom()
{
} // MeQuadBlossom::~MeQuadBlossom
//------------------------------------------------------------------------------
/// \brief Get the estimated time to run the Quad Blossom algorithm in minutes.
/// \param[in] a_numPoints The number of mesh points.
/// \return The estimated minutes to generate the quad mesh.
//------------------------------------------------------------------------------
double MeQuadBlossom::EstimatedRunTimeInMinutes(int a_numPoints)
{
  double minutes = 6.8E-11 * a_numPoints * a_numPoints * a_numPoints;
  return minutes;
} // MeQuadBlossom::EstimatedRunTimeInMinutes
//------------------------------------------------------------------------------
/// \brief Splits UGrid with 2D cells into quads by adding midpoints for each
/// edge and creating a quad by attaching adjacent midsides to the centroid.
/// \param[in] a_ugrid The input UGrid that contains only 2D cells.
/// \return The UGrid with cells split into quads.
//------------------------------------------------------------------------------
BSHP<XmUGrid> MeQuadBlossom::SplitToQuads(BSHP<XmUGrid> a_ugrid)
{
  VecPt3d points = a_ugrid->GetLocations();
  VecInt cells;

  int numPoints = a_ugrid->GetPointCount();
  VecVecAdjPointMidpoint edgeMidsides(numPoints);
  for (int pointIdx = 0; pointIdx < numPoints; ++pointIdx)
  {
    const Pt3d p0 = points[pointIdx];
    VecInt adjPoints;
    GetPointIdxsAttachedByEdge(a_ugrid, pointIdx, adjPoints);
    VecAdjPointMidpoint adjPoints2;
    for (auto otherIdx : adjPoints)
    {
      if (otherIdx > pointIdx)
      {
        const Pt3d& p1 = points[otherIdx];
        Pt3d midPoint = (p0 + p1) * 0.5;
        int midIdx = (int)points.size();
        points.push_back(midPoint);
        adjPoints2.push_back({ otherIdx, midIdx });
      }
    }
    // sort adjPoint2 by element 0
    std::sort(adjPoints2.begin(), adjPoints2.end());
    edgeMidsides[pointIdx] = adjPoints2;
  }

  // also get centroids and add to new points
  // track start of new centroid starts
  int numCells = a_ugrid->GetCellCount();
  for (int cellIdx = 0; cellIdx < numCells; ++cellIdx)
  {
    // TODO: Use new UGrid centroid code in xmsgrid
    VecInt cellPointIdxs = a_ugrid->GetCellPoints(cellIdx);
    VecPt3d cellPoints = a_ugrid->GetPointsLocations(cellPointIdxs);
    Pt3d centroid = gmComputeCentroid(cellPoints);
    int centroidIdx = (int)points.size();
    points.push_back(centroid);
    int numCellPoints = (int)cellPointIdxs.size();
    int lastIdx = cellPointIdxs.back();
    int midsideLast = FindMidpoint(edgeMidsides, cellPointIdxs.front(), lastIdx);
    for (size_t i = 0; i < numCellPoints; ++i)
    {
      int currentIdx = cellPointIdxs[i];
      int nextIdx = cellPointIdxs[(i + 1) % numCellPoints];
      int midside = FindMidpoint(edgeMidsides, currentIdx, nextIdx);

      ///    centroidIdx ---------- midside
      ///       |                       |
      ///    midsideLast ------------ currentIdx

      cells.push_back(XMU_QUAD);
      cells.push_back(4);
      cells.push_back(currentIdx);
      cells.push_back(midside);
      cells.push_back(centroidIdx);
      cells.push_back(midsideLast);

      midsideLast = midside;
    }
  }

  return XmUGrid::New(points, cells);
} // MeQuadBlossom::SplitToQuads

} // namespace xms

#if CXX_TEST
////////////////////////////////////////////////////////////////////////////////
// UNIT TESTS
////////////////////////////////////////////////////////////////////////////////

#include <xmsmesh/meshing/detail/MeQuadBlossom.t.h>

#include <xmscore/testing/TestTools.h>
#include <xmsinterp/triangulate/TrTriangulatorPoints.h>

//----- Namespace declaration --------------------------------------------------

using namespace xms;

////////////////////////////////////////////////////////////////////////////////
/// \class MeQuadBlossomUnitTests
/// \brief Contains unit tests for MeQuadBlossom and MeQuadBlossomImpl classes.
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief Test GetInteriorEdges function.
//------------------------------------------------------------------------------
void MeQuadBlossomUnitTests::testGetInteriorEdges()
{
  VecInt2d triangles;
  VecPt3d points;

  MeQuadBlossomImpl blossom(points, triangles);

  int p = 4;
  VecCellData chain = {{3, 0, 0}, {0, 1, 1}, {1, 2, 2}, {2, 5, 3}, {5, 8, 6}, {8, 7, 5}, {7, 3, 7}};
  VecEdge edges = blossom.GetInteriorEdges(p, chain);

  VecEdge expect = {{4, 5, 3, 6, 2, 8}, {4, 8, 6, 5, 5, 7}, {4, 7, 5, 7, 8, 3}};

  TS_ASSERT_EQUALS(expect, edges);

  chain = {{3, 0, 0}, {0, 1, 1}, {1, 2, 2}, {2, 5, 3}, {5, 8, 6}, {8, 7, 5}, {7, 3, 7}, {3, 0, 0}};
  edges = blossom.GetInteriorEdges(p, chain);

  // shouldn't chain result because last interior edge does not have increasing point indices;
  TS_ASSERT_EQUALS(expect, edges);

  chain = {{7, 3, 7}, {3, 0, 0}, {0, 1, 1}, {1, 2, 2}, {2, 5, 3}, {5, 8, 6}, {8, 7, 5}, {7, 3, 7}};
  edges = blossom.GetInteriorEdges(p, chain);

  expect = {{4, 5, 3, 6, 2, 8}, {4, 8, 6, 5, 5, 7}, {4, 7, 5, 7, 8, 3}};
  TS_ASSERT_EQUALS(expect, edges);

  chain = {{7, 3, 7}, {3, 0, 0}, {0, 1, 1}, {1, 2, 2}, {2, 5, 3}, {5, 8, 6}, {8, 7, 5}};
  edges = blossom.GetInteriorEdges(p, chain);

  expect = {{4, 5, 3, 6, 2, 8}, {4, 8, 6, 5, 5, 7}};
  TS_ASSERT_EQUALS(expect, edges);

  p = 5;
  chain = {{8, 4, 6}, {4, 2, 3}};
  edges = blossom.GetInteriorEdges(p, chain);
  expect = {};
  TS_ASSERT_EQUALS(expect, edges);

  p = 2;
  chain = {{5, 4, 3}, {4, 1, 2}};
  expect = {{2, 4, 3, 2, 5, 1}};
  edges = blossom.GetInteriorEdges(p, chain);
  TS_ASSERT_EQUALS(expect, edges);
} // MeQuadBlossomUnitTests::testGetInteriorEdges
//------------------------------------------------------------------------------
/// \brief Test GetBoundaryEdges function.
//------------------------------------------------------------------------------
void MeQuadBlossomUnitTests::testGetBoundaryEdges()
{
  VecInt2d triangles;
  VecPt3d points(8, Pt3d());
  MeQuadBlossomImpl blossom(points, triangles);

  int p = 0;
  VecCellData chain = {{1, 4, 1}, {4, 3, 0}};
  int cost = -1000;
  VecMeEdge extraEdges;
  VecInt2d extraPoints;
  Edge result = blossom.GetBoundaryEdges(p, chain, cost, extraEdges, extraPoints);
  Edge expect = {0, 3, 0, XM_NONE, 4, XM_NONE};

  TS_ASSERT_EQUALS(expect, result);
  TS_ASSERT(extraEdges.empty());
  TS_ASSERT(extraPoints.empty());

  p = 3;
  chain = {{0, 4, 0}, {4, 7, 7}, {7, 6, 4}};
  cost = -91;
  points.resize(8);
  result = blossom.GetBoundaryEdges(p, chain, cost, extraEdges, extraPoints);
  expect = {3, 6, 4, XM_NONE, 7, XM_NONE};
  VecMeEdge expectExtra = {{0, 4, -91}};
  VecInt2d expectExtraVertices = {{3, 4, 0, 6, 7}};

  TS_ASSERT_EQUALS(expect, result);
  TS_ASSERT_EQUALS(expectExtra, extraEdges);
  TS_ASSERT_EQUALS(expectExtraVertices, extraPoints);

  p = 6;
  chain = {{3, 7, 4}};
  extraEdges.clear();
  extraPoints.clear();
  points.resize(8);
  result = blossom.GetBoundaryEdges(p, chain, cost, extraEdges, extraPoints);
  expect = {6, 7, 4, XM_NONE, 3, XM_NONE};
  expectExtraVertices = {};

  TS_ASSERT_EQUALS(expect, result);
  TS_ASSERT(extraEdges.empty());
  TS_ASSERT(extraPoints.empty());

  p = 7;
  chain = {{4, 8, 5}};
  points.resize(8);
  result = blossom.GetBoundaryEdges(p, chain, cost, extraEdges, extraPoints);
  expect = {7, 8, 5, XM_NONE, 4, XM_NONE};

  TS_ASSERT_EQUALS(expect, result);
  TS_ASSERT(extraEdges.empty());
  TS_ASSERT(extraPoints.empty());

  chain = {{6, 3, 4}};
  points.resize(8);
  result = blossom.GetBoundaryEdges(p, chain, cost, extraEdges, extraPoints);
  expect = {7, 3, 4, XM_NONE, 6, XM_NONE};

  TS_ASSERT_EQUALS(expect, result);
  TS_ASSERT(extraEdges.empty());
  TS_ASSERT(extraPoints.empty());
} // MeQuadBlossomUnitTests::testGetBoundaryEdges
//------------------------------------------------------------------------------
/// \brief Test ProcessVertexChains function.
//------------------------------------------------------------------------------
void MeQuadBlossomUnitTests::testProcessVertexChains()
{
  VecInt2d triangles;
  VecPt3d points;
  MeQuadBlossomImpl blossom(points, triangles);

  int p = 4;
  VecVecCellData chains = {
    {{7, 3, 7}, {3, 0, 0}, {0, 1, 1}, {1, 2, 2}, {2, 5, 3}, {5, 8, 6}, {8, 7, 5}}};
  int extraCost = -1000;
  VecEdge interiorEdges;
  VecEdge boundaryEdges;
  VecMeEdge extraEdges;
  VecInt2d extraPoints;
  blossom.ProcessPointChains(p, chains, extraCost, interiorEdges, boundaryEdges, extraEdges,
                             extraPoints);
  VecEdge expectInterior = {{4, 5, 3, 6, 2, 8}, {4, 8, 6, 5, 5, 7}, {4, 7, 5, 7, 8, 3}};
  VecEdge expectBoundary;
  VecMeEdge expectExtra;
  VecInt2d expectExtraVertices;

  TS_ASSERT_EQUALS(expectInterior, interiorEdges);
  TS_ASSERT_EQUALS(expectBoundary, boundaryEdges);
  TS_ASSERT_EQUALS(expectExtra, extraEdges);
  TS_ASSERT_EQUALS(expectExtraVertices, extraPoints);

  p = 6;
  chains = {{{3, 7, 4}}};
  interiorEdges.clear();
  boundaryEdges.clear();
  extraEdges.clear();
  extraPoints.clear();
  blossom.ProcessPointChains(p, chains, extraCost, interiorEdges, boundaryEdges, extraEdges,
                             extraPoints);
  expectInterior = {};
  expectBoundary = {{6, 7, 4, XM_NONE, 3, XM_NONE}};
  expectExtra = {};
  expectExtraVertices = {};

  TS_ASSERT_EQUALS(expectInterior, interiorEdges);
  TS_ASSERT_EQUALS(expectBoundary, boundaryEdges);
  TS_ASSERT_EQUALS(expectExtra, extraEdges);
  TS_ASSERT_EQUALS(expectExtraVertices, extraPoints);

  p = 5;
  chains = {{{8, 4, 6}, {4, 2, 3}}};
  interiorEdges.clear();
  boundaryEdges.clear();
  extraEdges.clear();
  extraPoints.clear();
  blossom.ProcessPointChains(p, chains, extraCost, interiorEdges, boundaryEdges, extraEdges,
                             extraPoints);
  expectInterior = {};
  expectBoundary = {{5, 2, 3, XM_NONE, 4, XM_NONE}};
  expectExtra = {};
  expectExtraVertices = {};

  TS_ASSERT_EQUALS(expectInterior, interiorEdges);
  TS_ASSERT_EQUALS(expectBoundary, boundaryEdges);
  TS_ASSERT_EQUALS(expectExtra, extraEdges);
  TS_ASSERT_EQUALS(expectExtraVertices, extraPoints);

  p = 3;
  chains = {{{0, 4, 0}}, {{7, 13, 8}, {13, 6, 4}}};
  extraCost = 13;
  interiorEdges.clear();
  boundaryEdges.clear();
  extraEdges.clear();
  extraPoints.clear();
  blossom.ProcessPointChains(p, chains, extraCost, interiorEdges, boundaryEdges, extraEdges,
                             extraPoints);
  expectInterior = {{3, 13, 8, 4, 7, 6}};
  expectBoundary = {{3, 4, 0, XM_NONE, 0, XM_NONE}, {3, 6, 4, XM_NONE, 13, XM_NONE}};
  expectExtra = {};
  expectExtraVertices = {};

  TS_ASSERT_EQUALS(expectInterior, interiorEdges);
  TS_ASSERT_EQUALS(expectBoundary, boundaryEdges);
  TS_ASSERT_EQUALS(expectExtra, extraEdges);
  TS_ASSERT_EQUALS(expectExtraVertices, extraPoints);
} // MeQuadBlossomUnitTests::testProcessVertexChains
//------------------------------------------------------------------------------
/// \brief Test GetEtaAngle and GetEtaDistance functions.
//------------------------------------------------------------------------------
void MeQuadBlossomUnitTests::testEta()
{
  double h = 10.0 * sqrt(2.0);
  Pt3d p1 = {0, -h};
  Pt3d p2 = {0, h};
  Pt3d p0 = {-h, 0};
  Pt3d p3 = {h, 0};

  int angle_cost = GetEtaAngle(p0, p1, p2, p3);
  TS_ASSERT_EQUALS(1000, angle_cost);
  int distance_cost = GetEtaDistance(p0, p1, p2, p3);
  TS_ASSERT_EQUALS(1000, distance_cost);

  double a = XM_PI / 3.0;
  double hy = 10.0 * sin(a);
  double hx = 10 * cos(a);
  double a_d = a * 180.0 / XM_PI;
  p1 = {0, -hy};
  p2 = {0, hy};
  p0 = {-hx, 0};
  p3 = {hx, 0};

  angle_cost = GetEtaAngle(p0, p1, p2, p3);
  TS_ASSERT_EQUALS(667, angle_cost);
  distance_cost = GetEtaDistance(p0, p1, p2, p3);
  TS_ASSERT_EQUALS(500, distance_cost);

  angle_cost = GetEtaAngle(p1, p0, p3, p2);
  TS_ASSERT_EQUALS(667, angle_cost);
  distance_cost = GetEtaDistance(p1, p0, p3, p2);
  TS_ASSERT_EQUALS(500, distance_cost);

  a = XM_PI / 4.0;
  hy = 10.0 * sin(a);
  hx = 10.0 * cos(a);
  a_d = a * 180.0 / XM_PI;
  p1 = {0, -hy};
  p2 = {0, hy};
  p0 = {-hx, 0};
  p3 = {hx, 0};

  angle_cost = GetEtaAngle(p0, p1, p2, p3);
  TS_ASSERT_EQUALS(1000, angle_cost);
  distance_cost = GetEtaDistance(p0, p1, p2, p3);
  TS_ASSERT_EQUALS(1000, distance_cost);

  a = XM_PI / 6.0;
  hy = 10.0 * sin(a);
  hx = 10.0 * cos(a);
  a_d = a * 180.0 / XM_PI;
  p1 = {0, -hy};
  p2 = {0, hy};
  p0 = {-hx, 0};
  p3 = {hx, 0};

  angle_cost = GetEtaAngle(p0, p1, p2, p3);
  TS_ASSERT_EQUALS(667, angle_cost);
  distance_cost = GetEtaDistance(p0, p1, p2, p3);
  TS_ASSERT_EQUALS(500, distance_cost);

  a = XM_PI / 8.0;
  hy = 10.0 * sin(a);
  hx = 10.0 * cos(a);
  a_d = a * 180.0 / XM_PI;
  p1 = {0, -hy};
  p2 = {0, hy};
  p0 = {-hx, 0};
  p3 = {hx, 0};

  angle_cost = GetEtaAngle(p0, p1, p2, p3);
  TS_ASSERT_EQUALS(500, angle_cost);
  distance_cost = GetEtaDistance(p0, p1, p2, p3);
  TS_ASSERT_EQUALS(293, distance_cost);

  angle_cost = GetEtaAngle(p1, p0, p3, p2);
  TS_ASSERT_EQUALS(500, angle_cost);
  distance_cost = GetEtaDistance(p1, p0, p3, p2);
  TS_ASSERT_EQUALS(293, distance_cost);

  a = XM_PI / 12.0;
  hy = sin(a);
  hx = cos(a);
  a_d = a * 180.0 / XM_PI;
  p1 = {0, -hy};
  p2 = {0, hy};
  p0 = {-hx, 0};
  p3 = {hx, 0};

  angle_cost = GetEtaAngle(p0, p1, p2, p3);
  TS_ASSERT_EQUALS(333, angle_cost);
  distance_cost = GetEtaDistance(p0, p1, p2, p3);
  TS_ASSERT_EQUALS(134, distance_cost);

  angle_cost = GetEtaAngle(p1, p0, p3, p2);
  TS_ASSERT_EQUALS(333, angle_cost);
  distance_cost = GetEtaDistance(p1, p0, p3, p2);
  TS_ASSERT_EQUALS(134, distance_cost);

  // Larger scale does not matter. Relative scale gives same result.
  p1 = {0, -hy * 100.0};
  p2 = {0, hy * 100.0};
  p0 = {-hx * 100.0, 0};
  p3 = {hx * 100.0, 0};

  angle_cost = GetEtaAngle(p0, p1, p2, p3);
  TS_ASSERT_EQUALS(333, angle_cost);
  distance_cost = GetEtaDistance(p0, p1, p2, p3);
  TS_ASSERT_EQUALS(134, distance_cost);

  angle_cost = GetEtaAngle(p1, p0, p3, p2);
  TS_ASSERT_EQUALS(333, angle_cost);
  distance_cost = GetEtaDistance(p1, p0, p3, p2);
  TS_ASSERT_EQUALS(134, distance_cost);

  // Smaller scale does not matter. Relative scale gives same result.
  p1 = {0, -hy / 100.0};
  p2 = {0, hy / 100.0};
  p0 = {-hx / 100.0, 0};
  p3 = {hx / 100.0, 0};

  angle_cost = GetEtaAngle(p0, p1, p2, p3);
  TS_ASSERT_EQUALS(333, angle_cost);
  distance_cost = GetEtaDistance(p0, p1, p2, p3);
  TS_ASSERT_EQUALS(134, distance_cost);

  angle_cost = GetEtaAngle(p1, p0, p3, p2);
  TS_ASSERT_EQUALS(333, angle_cost);
  distance_cost = GetEtaDistance(p1, p0, p3, p2);
  TS_ASSERT_EQUALS(134, distance_cost);
} // MeQuadBlossomUnitTests::testEta
//------------------------------------------------------------------------------
/// \brief Test GetAngle function.
//------------------------------------------------------------------------------
void MeQuadBlossomUnitTests::testGetAngle()
{
  Pt3d pc = {17, -31};
  Pt3d p0 = {17 + 8, -31};
  Pt3d p1 = {17 + 8, -31 + 0};
  double angle = GetAngle(p0, pc, p1);
  TS_ASSERT_DELTA(0.0, angle, 1.0e-5);
  angle = GetAngle(p1, pc, p0);
  TS_ASSERT_DELTA(0.0, angle, 1.0e-5);
  p1 = {17 + 8, -31 + 8};
  angle = GetAngle(p0, pc, p1);
  TS_ASSERT_DELTA(XM_PI * 1.75, angle, 1.0e-5);
  angle = GetAngle(p1, pc, p0);
  TS_ASSERT_DELTA(XM_PI * 0.25, angle, 1.0e-5);
  p1 = {17 + 0, -31 + 8};
  angle = GetAngle(p0, pc, p1);
  TS_ASSERT_DELTA(XM_PI * 1.5, angle, 1.0e-5);
  angle = GetAngle(p1, pc, p0);
  TS_ASSERT_DELTA(XM_PI * 0.5, angle, 1.0e-5);
  p1 = {17 - 8, -31 + 8};
  angle = GetAngle(p0, pc, p1);
  TS_ASSERT_DELTA(XM_PI * 1.25, angle, 1.0e-5);
  angle = GetAngle(p1, pc, p0);
  TS_ASSERT_DELTA(XM_PI * 0.75, angle, 1.0e-5);
  p1 = {17 - 8, -31 + 0};
  angle = GetAngle(p0, pc, p1);
  TS_ASSERT_DELTA(XM_PI, angle, 1.0e-5);
  angle = GetAngle(p1, pc, p0);
  TS_ASSERT_DELTA(XM_PI, angle, 1.0e-5);
  p1 = {17 - 8, -31 - 8};
  angle = GetAngle(p0, pc, p1);
  TS_ASSERT_DELTA(XM_PI * 0.75, angle, 1.0e-5);
  angle = GetAngle(p1, pc, p0);
  TS_ASSERT_DELTA(XM_PI * 1.25, angle, 1.0e-5);
  p1 = {17 - 0, -31 - 8};
  angle = GetAngle(p0, pc, p1);
  TS_ASSERT_DELTA(XM_PI * 0.5, angle, 1.0e-5);
  angle = GetAngle(p1, pc, p0);
  TS_ASSERT_DELTA(XM_PI * 1.5, angle, 1.0e-5);
  p1 = {17 + 8, -31 - 8};
  angle = GetAngle(p0, pc, p1);
  TS_ASSERT_DELTA(XM_PI * 0.25, angle, 1.0e-5);
  angle = GetAngle(p1, pc, p0);
  TS_ASSERT_DELTA(XM_PI * 1.75, angle, 1.0e-5);
} // MeQuadBlossomUnitTests::testGetAngle
//------------------------------------------------------------------------------
/// \brief Test ExtractCellsAdjacentToPoint function.
//------------------------------------------------------------------------------
void MeQuadBlossomUnitTests::testExtractCellsAdjacentToPoint()
{
  VecInt2d triangles = {{0, 3, 4}, {0, 4, 1}, {1, 4, 2}, {2, 4, 5},
                        {3, 6, 7}, {4, 7, 8}, {4, 8, 5}, {3, 7, 4}};
  VecPt3d points(9, Pt3d());

  MeQuadBlossomImpl blossom(points, triangles);
  VecVecCellData vertexFaces;
  VecVecVecCellData sortedFacesChains;
  blossom.ExtractCellsAdjacentToPoint(vertexFaces, sortedFacesChains);

  VecVecCellData expected = {
    {{4, 3, 0}, {1, 4, 1}},
    {{4, 0, 1}, {2, 4, 2}},
    {{4, 1, 2}, {5, 4, 3}},
    {{0, 4, 0}, {7, 6, 4}, {4, 7, 7}},
    {{3, 0, 0}, {0, 1, 1}, {1, 2, 2}, {2, 5, 3}, {8, 7, 5}, {5, 8, 6}, {7, 3, 7}},
    {{4, 2, 3}, {8, 4, 6}},
    {{3, 7, 4}},
    {{6, 3, 4}, {4, 8, 5}, {3, 4, 7}},
    {{7, 4, 5}, {4, 5, 6}}};

  VecVecVecCellData expectedSorted = {
    {{{1, 4, 1}, {4, 3, 0}}},
    {{{2, 4, 2}, {4, 0, 1}}},
    {{{5, 4, 3}, {4, 1, 2}}},
    {{{0, 4, 0}, {4, 7, 7}, {7, 6, 4}}},
    {{{3, 0, 0}, {0, 1, 1}, {1, 2, 2}, {2, 5, 3}, {5, 8, 6}, {8, 7, 5}, {7, 3, 7}}},
    {{{8, 4, 6}, {4, 2, 3}}},
    {{{3, 7, 4}}},
    {{{6, 3, 4}, {3, 4, 7}, {4, 8, 5}}},
    {{{7, 4, 5}, {4, 5, 6}}}};

  TS_ASSERT_EQUALS(expected, vertexFaces);
  TS_ASSERT_EQUALS(expectedSorted, sortedFacesChains);

  triangles = {{0, 3, 4}, {0, 4, 1}, {1, 4, 2}, {2, 4, 5}, {3, 6, 7}, {4, 7, 8}, {4, 8, 5}};

  MeQuadBlossomImpl blossom2(points, triangles);
  blossom2.ExtractCellsAdjacentToPoint(vertexFaces, sortedFacesChains);

  expected = {{{4, 3, 0}, {1, 4, 1}},
              {{4, 0, 1}, {2, 4, 2}},
              {{4, 1, 2}, {5, 4, 3}},
              {{0, 4, 0}, {7, 6, 4}},
              {{3, 0, 0}, {0, 1, 1}, {1, 2, 2}, {2, 5, 3}, {8, 7, 5}, {5, 8, 6}},
              {{4, 2, 3}, {8, 4, 6}},
              {{3, 7, 4}},
              {{6, 3, 4}, {4, 8, 5}},
              {{7, 4, 5}, {4, 5, 6}}};

  // p: {p, p01=p10, f1, f0, f
  // 0: {0, 4, 0, 1, 3, 1}  boundary: 1, 0, 3: {0, 3, 0, 4}, {1, 0, 1, 4}
  // 1: {1, 4, 1, 2, 0, 2}  boundary:
  expectedSorted = {{{{1, 4, 1}, {4, 3, 0}}},
                    {{{2, 4, 2}, {4, 0, 1}}},
                    {{{5, 4, 3}, {4, 1, 2}}},
                    {{{0, 4, 0}}, {{7, 6, 4}}},
                    {{{3, 0, 0}, {0, 1, 1}, {1, 2, 2}, {2, 5, 3}, {5, 8, 6}, {8, 7, 5}}},
                    {{{8, 4, 6}, {4, 2, 3}}},
                    {{{3, 7, 4}}},
                    {{{6, 3, 4}}, {{4, 8, 5}}},
                    {{{7, 4, 5}, {4, 5, 6}}}};

  TS_ASSERT_EQUALS(expected, vertexFaces);
  TS_ASSERT_EQUALS(expectedSorted, sortedFacesChains);
} // MeQuadBlossomUnitTests::testExtractCellsAdjacentToPoint
//------------------------------------------------------------------------------
/// \brief Test GetAngle function.
//------------------------------------------------------------------------------
void MeQuadBlossomUnitTests::testGetEdges()
{
  {
    // single triangle
    VecInt2d triangles = {{2, 9, 6}};
    VecPt3d points(10, Pt3d());

    MeQuadBlossomImpl blossom(points, triangles);
    blossom.GetEdges(17);

    VecEdge expect_interior = {};
    VecEdge expect_boundary = {
      {9, 6, 0, XM_NONE, 2, XM_NONE},
      {6, 2, 0, XM_NONE, 9, XM_NONE},
      {2, 9, 0, XM_NONE, 6, XM_NONE},
    };
    VecMeEdge expect_extra = {};
    VecInt2d expectExtraVertices = {};
    TS_ASSERT_EQUALS(expect_interior, blossom.m_interiorEdges);
    TS_ASSERT_EQUALS(expect_boundary, blossom.m_boundaryEdges);
    TS_ASSERT_EQUALS(expect_extra, blossom.m_extraEdges);
    TS_ASSERT_EQUALS(expectExtraVertices, blossom.m_extraPoints);
  }

  {
    // two adjacent
    VecInt2d triangles = {{2, 9, 6}, {9, 4, 6}};
    VecPt3d points(10, Pt3d());

    MeQuadBlossomImpl blossom(points, triangles);
    blossom.GetEdges(17);

    VecEdge expect_interior = {{6, 9, 1, 0, 4, 2}};
    VecEdge expect_boundary = {{9, 4, 1, XM_NONE, 6, XM_NONE},
                               {6, 2, 0, XM_NONE, 9, XM_NONE},
                               {4, 6, 1, XM_NONE, 9, XM_NONE},
                               {2, 9, 0, XM_NONE, 6, XM_NONE}};
    VecMeEdge expect_extra = {};
    VecInt2d expectExtraVertices = {};
    TS_ASSERT_EQUALS(expect_interior, blossom.m_interiorEdges);
    TS_ASSERT_EQUALS(expect_boundary, blossom.m_boundaryEdges);
    TS_ASSERT_EQUALS(expect_extra, blossom.m_extraEdges);
    TS_ASSERT_EQUALS(expectExtraVertices, blossom.m_extraPoints);
  }

  {
    // three adjacent
    VecInt2d triangles = {{2, 9, 6}, {9, 4, 6}, {9, 3, 4}};
    VecPt3d points(10, Pt3d());

    MeQuadBlossomImpl blossom(points, triangles);
    blossom.GetEdges(17);

    VecEdge expect_interior = {
      {6, 9, 1, 0, 4, 2},
      {4, 9, 2, 1, 3, 6},
    };
    VecEdge expect_boundary = {
      {9, 3, 2, XM_NONE, 4, XM_NONE}, {6, 2, 0, XM_NONE, 9, XM_NONE},
      {4, 6, 1, XM_NONE, 9, XM_NONE}, {3, 4, 2, XM_NONE, 9, XM_NONE},
      {2, 9, 0, XM_NONE, 6, XM_NONE},
    };
    VecMeEdge expect_extra = {{0, 2, 17}};
    VecInt2d expectExtraVertices = {{9, 6, 2, 3, 4}};
    TS_ASSERT_EQUALS(expect_interior, blossom.m_interiorEdges);
    TS_ASSERT_EQUALS(expect_boundary, blossom.m_boundaryEdges);
    TS_ASSERT_EQUALS(expect_extra, blossom.m_extraEdges);
    TS_ASSERT_EQUALS(expectExtraVertices, blossom.m_extraPoints);
  }

  {
    // 3x3 grid
    VecInt2d triangles = {{4, 3, 0}, {1, 4, 0}, {2, 4, 1}, {5, 4, 2},
                          {7, 6, 3}, {8, 7, 4}, {5, 8, 4}, {4, 7, 3}};
    VecPt3d points(9, Pt3d());

    MeQuadBlossomImpl blossom(points, triangles);
    blossom.GetEdges(17);
    VecEdge expect_interior = {
      {4, 7, 7, 5, 3, 8}, // 0
      {4, 8, 5, 6, 7, 5}, // 1
      {4, 5, 6, 3, 8, 2}, // 2
      {3, 7, 4, 7, 6, 4}, // 3
      {3, 4, 7, 0, 7, 0}, // 4
      {2, 4, 2, 3, 1, 5}, // 5
      {1, 4, 1, 2, 0, 2}, // 6
      {0, 4, 0, 1, 3, 1}, // 7
    };
    VecEdge expect_boundary = {
      {8, 7, 5, XM_NONE, 4, XM_NONE}, // 0
      {7, 6, 4, XM_NONE, 3, XM_NONE}, // 1
      {6, 3, 4, XM_NONE, 7, XM_NONE}, // 2
      {5, 8, 6, XM_NONE, 4, XM_NONE}, // 3
      {3, 0, 0, XM_NONE, 4, XM_NONE}, // 4
      {2, 5, 3, XM_NONE, 4, XM_NONE}, // 5
      {1, 2, 2, XM_NONE, 4, XM_NONE}, // 6
      {0, 1, 1, XM_NONE, 4, XM_NONE}, // 7
    };
    VecMeEdge expect_extra = {{4, 0, 17}, {5, 4, 17}};
    VecInt2d expectExtraVertices = {{3, 7, 6, 0, 4}, {7, 4, 8, 6, 3}};
    TS_ASSERT_EQUALS(expect_interior, blossom.m_interiorEdges);
    TS_ASSERT_EQUALS(expect_boundary, blossom.m_boundaryEdges);
    TS_ASSERT_EQUALS(expect_extra, blossom.m_extraEdges);
    TS_ASSERT_EQUALS(expectExtraVertices, blossom.m_extraPoints);
  }

  {
    // 3x3 grid with a hole where face { 3, 4, 7} was in the earlier test
    // and with { 3, 7, 6 } split into { 3, 7, 9 } and { 3, 9, 6}
    VecInt2d triangles = {{4, 3, 0}, {1, 4, 0}, {2, 4, 1}, {5, 4, 2},
                          {3, 7, 9}, {8, 7, 4}, {5, 8, 4}, {3, 9, 6}};
    VecPt3d points(10, Pt3d());
    
    MeQuadBlossomImpl blossom(points, triangles);
    blossom.GetEdges(17);
    VecEdge expect_interior = {
      {4, 8, 5, 6, 7, 5}, // 0
      {4, 5, 6, 3, 8, 2}, // 1
      {3, 9, 7, 4, 6, 7}, // 2
      {2, 4, 2, 3, 1, 5}, // 3
      {1, 4, 1, 2, 0, 2}, // 4
      {0, 4, 0, 1, 3, 1}, // 5
    };
    VecEdge expect_boundary = {
      {9, 6, 7, XM_NONE, 3, XM_NONE}, // 0
      {8, 7, 5, XM_NONE, 4, XM_NONE}, // 1
      {7, 9, 4, XM_NONE, 3, XM_NONE}, // 2
      {7, 4, 5, XM_NONE, 8, XM_NONE}, // 3
      {6, 3, 7, XM_NONE, 9, XM_NONE}, // 4
      {5, 8, 6, XM_NONE, 4, XM_NONE}, // 5
      {4, 3, 0, XM_NONE, 0, XM_NONE}, // 6
      {3, 0, 0, XM_NONE, 4, XM_NONE}, // 7
      {3, 7, 4, XM_NONE, 9, XM_NONE}, // 8
      {2, 5, 3, XM_NONE, 4, XM_NONE}, // 9
      {1, 2, 2, XM_NONE, 4, XM_NONE}, // 10
      {0, 1, 1, XM_NONE, 4, XM_NONE}, // 11
    };
    VecMeEdge expect_extra = {{5, 0, 17}};
    VecInt2d expectExtraVertices = {{4, 8, 7, 3, 0}};
    TS_ASSERT_EQUALS(expect_interior, blossom.m_interiorEdges);
    TS_ASSERT_EQUALS(expect_boundary, blossom.m_boundaryEdges);
    TS_ASSERT_EQUALS(expect_extra, blossom.m_extraEdges);
    TS_ASSERT_EQUALS(expectExtraVertices, blossom.m_extraPoints);
  }
} // MeQuadBlossomUnitTests::testGetEdges
//------------------------------------------------------------------------------
/// \brief Test GetAngle function.
//------------------------------------------------------------------------------
void MeQuadBlossomUnitTests::testEliminateEdges()
{
  {
    VecInt2d triangles = {{4, 3, 0}, {1, 4, 0}, {2, 4, 1}, {5, 4, 2},
                          {7, 6, 3}, {8, 7, 4}, {5, 8, 4}, {4, 7, 3}};
    VecPt3d points(9, Pt3d());

    MeQuadBlossomImpl blossom(points, triangles);
    blossom.GetEdges(17);
    VecInt2d quads = blossom.EliminateEdges({1, 3, 5, 7});
    VecInt2d expect = {{0, 1, 4, 3}, {2, 5, 4, 1}, {3, 4, 7, 6}, {4, 5, 8, 7}};
    TS_ASSERT_EQUALS(expect, quads);
  }

  {
    VecInt2d triangles = {{4, 3, 0}, {1, 4, 0}, {2, 4, 1}, {5, 4, 2},
                          {3, 7, 9}, {8, 7, 4}, {5, 8, 4}, {3, 9, 6}};
    VecPt3d points(10, Pt3d());
    MeQuadBlossomImpl blossom(points, triangles);
    blossom.GetEdges(17);
    VecInt2d quads = blossom.EliminateEdges({2, 5, 3, 0});
    
    VecInt2d expect = {{0, 1, 4, 3}, {2, 5, 4, 1}, {4, 5, 8, 7}, {3, 7, 9, 6}};
    TS_ASSERT_EQUALS(expect, quads);
  }
} // MeQuadBlossomUnitTests::testEliminateEdges
//------------------------------------------------------------------------------
/// \brief Test simple triangle with three divisions on each side.
//------------------------------------------------------------------------------
void MeQuadBlossomUnitTests::testSimpleTriangle()
{
  // clang-format off
  VecPt3d points = {{0, 0, 0},  {10, 0, 0},  {20, 0, 0}, {30, 0, 0},
                    {0, 10, 0}, {10, 10, 0}, {20, 10, 0},
                    {0, 20, 0}, {10, 20, 0},
                    {0, 30, 0}};
  // clang-format on
  VecInt2d triangles = {{0, 1, 4}, {1, 5, 4}, {1, 2, 5}, {2, 6, 5}, {2, 3, 6},
                        {4, 5, 7}, {5, 8, 7}, {5, 6, 8}, {7, 8, 9}};
  // VecInt2d edges = {{0, 1, 15}, {1, 2, 10}, {2, 3, 15}, {3, 4, 10},
  //                   {1, 5, 10}, {3, 7, 10}, {5, 6, 15}, {6, 7, 10},
  //                   {6, 8, 10},
  //                   {0, 2, -10}, {2, 4, -10}, {4, 7, -10}, {7, 8, -10},
  //                   {8, 5, -10}, {5, 0, -10}};
  {
    MeQuadBlossomImpl blossom(points, triangles);
    int cost = -10;
    bool splitBoundaryPoints = false;
    bool useAngleCost = true;
    BSHP<XmUGrid> ugrid = blossom._MakeQuads(splitBoundaryPoints, useAngleCost, cost);
    const VecInt& cells = ugrid->GetCellstream();
    VecInt expectedCells = {
      XMU_QUAD, 4, 1, 5, 4, 0,  // 0
      XMU_QUAD, 4, 2, 6, 5, 1,  // 1
      XMU_TRIANGLE, 3, 2, 3, 6, // 2
      XMU_QUAD, 4, 5, 8, 7, 4,  // 3
      XMU_TRIANGLE, 3, 5, 6, 8, // 4
      XMU_TRIANGLE, 3, 7, 8, 9  // 5
    };
    TS_ASSERT_EQUALS(expectedCells, cells);
    const VecPt3d& actualPoints = ugrid->GetLocations();
    TS_ASSERT_DELTA_VECPT3D(points, actualPoints, 1.0e-5);
  }

  {
    MeQuadBlossomImpl blossom(points, triangles);
    int cost = -10;
    bool splitBoundaryPoints = true;
    bool useAngleCost = true;
    BSHP<XmUGrid> ugrid = blossom._MakeQuads(splitBoundaryPoints, useAngleCost, cost);
    const VecInt& cells = ugrid->GetCellstream();
    VecInt expectedCells = {
      XMU_QUAD, 4, 1, 5, 4, 0,  // 0
      XMU_QUAD, 4, 2, 10, 5, 1, // 1
      XMU_QUAD, 4, 2, 3, 6, 10, // 2
      XMU_QUAD, 4, 5, 8, 7, 4,  // 3
      XMU_QUAD, 4, 6, 8, 5, 10, // 4
      XMU_TRIANGLE, 3, 7, 8, 9  // 5
    };
    TS_ASSERT_EQUALS(expectedCells, cells);
    const VecPt3d& actualPoints = ugrid->GetLocations();
    VecPt3d expectedPoints = points;
    expectedPoints.push_back(Pt3d(50.0/3.0, 20.0/3.0, 0.0));
    TS_ASSERT_DELTA_VECPT3D(expectedPoints, actualPoints, 1.0e-5);
  }

  {
    // clang-format off
    VecPt3d points = { {-10, 0, 0}, {0, 0, 0}, {10, 0, 0}, {-5, 7.07, 0}, {5, 7.07, 0}, {0, 14.14, 0} };
    // clang-format on
    VecInt2d triangles = { {0, 1, 3}, {1, 2, 4}, {1, 4, 3}, {3, 4, 5} };

    MeQuadBlossomImpl blossom(points, triangles);
    int cost = -10;
    bool splitBoundaryPoints = true;
    bool useAngleCost = false;
    BSHP<XmUGrid> ugrid = blossom._MakeQuads(splitBoundaryPoints, useAngleCost, cost);
    const VecInt& cells = ugrid->GetCellstream();
    VecInt expectedCells = {
      XMU_QUAD, 4, 3, 0, 1, 6,  // 0
      XMU_QUAD, 4, 1, 2, 4, 6,  // 1
      XMU_QUAD, 4, 3, 6, 4, 5   // 2
    };
    TS_ASSERT_EQUALS(expectedCells, cells);
    const VecPt3d& actualPoints = ugrid->GetLocations();
    VecPt3d expectedPoints = points;
    expectedPoints.push_back(Pt3d(0, 14.14 / 3, 0));
    TS_ASSERT_DELTA_VECPT3D(expectedPoints, actualPoints, 1.0e-5);
  }
} // MeQuadBlossomUnitTests::testSimpleTriangle
//------------------------------------------------------------------------------
/// \brief Test simple quad with three divisions on each side.
//------------------------------------------------------------------------------
void MeQuadBlossomUnitTests::testSimpleQuad()
{
  // clang-format off
  VecPt3d points = {
    {0, 0, 0}, {10, 0, 0}, {20, 0, 0}, {30, 0, 0},
    {0, 10, 0}, {10, 10, 0}, {20, 10, 0}, {30, 10, 0},
    {0, 20, 0}, {10, 20, 0}, {20, 20, 0}, {30, 20, 0},
    {0, 30, 0}, {10, 30, 0}, {20, 30, 0}, {30, 30, 0}};
  VecInt2d triangles = {
    {0, 1, 4}, {1, 5, 4}, {1, 2, 5}, {2, 6, 5}, {2, 3, 6}, {3, 7, 6},
    {4, 9, 8}, {4, 5, 9}, {5, 10, 9}, {5, 6, 10}, {6, 7, 11}, {6, 11, 10},
    {8, 9, 12}, {9, 13, 12}, {9, 10, 13}, {10, 14, 13}, {10, 11, 14},
    {11, 15, 14}};

  {
    MeQuadBlossomImpl blossom(points, triangles);
    int cost = -10;
    bool splitBoundaryPoints = false;
    bool useAngleCost = true;
    BSHP<XmUGrid> ugrid = blossom._MakeQuads(splitBoundaryPoints, useAngleCost, cost);
    const VecInt& cells = ugrid->GetCellstream();
    VecInt expectedCells = {
      XMU_QUAD, 4, 1, 5, 4, 0,     // 0
      XMU_QUAD, 4, 2, 6, 5, 1,     // 1
      XMU_QUAD, 4, 3, 7, 6, 2,     // 2
      XMU_QUAD, 4, 4, 5, 9, 8,     // 3
      XMU_QUAD, 4, 5, 6, 10, 9,    // 4
      XMU_QUAD, 4, 6, 7, 11, 10,   // 5
      XMU_QUAD, 4, 9, 13, 12, 8,   // 6
      XMU_QUAD, 4, 10, 14, 13, 9,  // 7
      XMU_QUAD, 4, 11, 15, 14, 10, // 8
    };
    TS_ASSERT_EQUALS(expectedCells, cells);
    const VecPt3d& actualPoints = ugrid->GetLocations();
    VecPt3d expectedPoints = points;
    TS_ASSERT_DELTA_VECPT3D(expectedPoints, actualPoints, 1.0e-5);
  }

  {
    MeQuadBlossomImpl blossom(points, triangles);
    int cost = -10;
    bool splitBoundaryPoints = true;
    bool useAngleCost = true;
    BSHP<XmUGrid> ugrid = blossom._MakeQuads(splitBoundaryPoints, useAngleCost, cost);
    const VecInt& cells = ugrid->GetCellstream();
    VecInt expectedCells = {
      XMU_QUAD, 4, 1, 5, 4, 0,     // 0
      XMU_QUAD, 4, 2, 6, 5, 1,     // 1
      XMU_QUAD, 4, 3, 7, 6, 2,     // 2
      XMU_QUAD, 4, 4, 5, 9, 8,     // 3
      XMU_QUAD, 4, 5, 6, 10, 9,    // 4
      XMU_QUAD, 4, 6, 7, 11, 10,   // 5
      XMU_QUAD, 4, 9, 13, 12, 8,   // 6
      XMU_QUAD, 4, 10, 14, 13, 9,  // 7
      XMU_QUAD, 4, 11, 15, 14, 10, // 8
    };
    TS_ASSERT_EQUALS(expectedCells, cells);
    const VecPt3d& actualPoints = ugrid->GetLocations();
    VecPt3d expectedPoints = points;
    TS_ASSERT_DELTA_VECPT3D(expectedPoints, actualPoints, 1.0e-5);
  }
} // MeQuadBlossomUnitTests::testSimpleQuad
//------------------------------------------------------------------------------
/// \brief Test complex quad with three divisions on two opposite sides and
//  two on the horizontal bottom row increasing to six divisions on top row.
//------------------------------------------------------------------------------
void MeQuadBlossomUnitTests::testComplexQuad()
{
  // clang-format off
  VecPt3d points = {
    {-10, 0, 0}, {0, 0, 0}, {10, 0,0},
    {-15, 10, 0}, {-5, 10, 0}, {5, 10, 0}, {15, 10, 0},
    {-20, 20, 0}, {-10, 20, 0}, {0, 20, 0}, {10, 20, 0}, {20, 20, 0},
    {-25, 30, 0}, {-15, 30, 0}, {-5, 30, 0}, {5, 30, 0}, {15, 30, 0}, {25, 30, 0},
    {-30, 40, 0}, {-20, 40, 0}, {-10, 40, 0}, {0, 40, 0}, {10, 40, 0}, {20, 40, 0}, {30, 40, 0}};
  VecInt2d triangles = {
    {0, 4, 3}, {0, 1, 4}, {1, 5, 4}, {1, 2, 5}, {2, 6, 5},
    {3, 8, 7}, {3, 4, 8}, {4, 9, 8}, {4, 5, 9}, {5, 10, 9}, {5, 11, 10}, {5, 6, 11},
    {7, 13, 12}, {7, 8, 13}, {8, 14, 13}, {8, 9, 14}, {9, 15, 14}, {9, 10, 15}, {10, 16, 15}, {10, 11, 16}, {11, 17, 16},
    {12, 19, 18}, {12, 13, 19}, {13, 20, 19}, {13, 14, 20}, {14, 21, 20}, {14, 15, 21}, {15, 22, 21}, {15, 23, 22}, {15, 16, 23}, {16, 17, 23}, {17, 24, 23}};

  {
    MeQuadBlossomImpl blossom(points, triangles);
    int cost = -10;
    bool splitBoundaryPoints = false;
    bool useAngleCost = true;
    int numBoundaryEdges = blossom.PreMakeQuads();
    TS_ASSERT_EQUALS(16, numBoundaryEdges);
    BSHP<XmUGrid> ugrid = blossom._MakeQuads(splitBoundaryPoints, useAngleCost, cost);
    const VecInt& cells = ugrid->GetCellstream();
    VecInt expectedCells = {
      XMU_TRIANGLE, 3, 0, 4, 3,    //  0
      XMU_QUAD, 4, 1, 5, 4, 0,     //  1
      XMU_QUAD, 4, 2, 6, 5, 1,     //  2
      XMU_TRIANGLE, 3, 3, 8, 7,    //  3
      XMU_QUAD, 4, 4, 9, 8, 3,     //  4
      XMU_QUAD, 4, 5, 10, 9, 4,    //  5
      XMU_QUAD, 4, 5, 6, 11, 10,   //  6
      XMU_TRIANGLE, 3, 7, 13, 12,  //  7
      XMU_QUAD, 4, 8, 14, 13, 7,   //  8
      XMU_QUAD, 4, 9, 15, 14, 8,   //  9
      XMU_QUAD, 4, 10, 16, 15, 9,  // 10
      XMU_QUAD, 4, 11, 17, 16, 10, // 11
      XMU_TRIANGLE, 3, 12, 19, 18, // 12
      XMU_QUAD, 4, 13, 20, 19, 12, // 13
      XMU_QUAD, 4, 14, 21, 20, 13, // 14
      XMU_QUAD, 4, 15, 22, 21, 14, // 15
      XMU_QUAD, 4, 15, 16, 23, 22, // 16
      XMU_QUAD, 4, 17, 24, 23, 16  // 17
    };
    TS_ASSERT_EQUALS(expectedCells, cells);
    const VecPt3d& actualPoints = ugrid->GetLocations();
    VecPt3d expectedPoints = points;
    TS_ASSERT_DELTA_VECPT3D(expectedPoints, actualPoints, 1.0e-5);
  }

  {
    MeQuadBlossomImpl blossom(points, triangles);
    int cost = -10;
    bool splitBoundaryPoints = true;
    bool useAngleCost = true;
    BSHP<XmUGrid> ugrid = blossom._MakeQuads(splitBoundaryPoints, useAngleCost, cost);
    const VecInt& cells = ugrid->GetCellstream();
    VecInt expectedCells = {
      XMU_QUAD, 4, 3, 0, 4, 25,    //  0
      XMU_QUAD, 4, 1, 5, 4, 0,     //  1
      XMU_QUAD, 4, 2, 6, 5, 1,     //  2
      XMU_QUAD, 4, 8, 7, 3, 25,    //  3
      XMU_QUAD, 4, 4, 9, 8, 25,    //  4
      XMU_QUAD, 4, 5, 10, 9, 4,    //  5
      XMU_QUAD, 4, 5, 6, 11, 10,   //  6
      XMU_QUAD, 4, 12, 7, 13, 26,  //  7
      XMU_QUAD, 4, 8, 14, 13, 7,   //  8
      XMU_QUAD, 4, 9, 15, 14, 8,   //  9
      XMU_QUAD, 4, 10, 16, 15, 9,  // 10
      XMU_QUAD, 4, 11, 17, 16, 10, // 11
      XMU_QUAD, 4, 19, 18, 12, 26, // 12
      XMU_QUAD, 4, 13, 20, 19, 26, // 13
      XMU_QUAD, 4, 14, 21, 20, 13, // 14
      XMU_QUAD, 4, 15, 22, 21, 14, // 15
      XMU_QUAD, 4, 15, 16, 23, 22, // 16
      XMU_QUAD, 4, 17, 24, 23, 16, // 17
    };
    TS_ASSERT_EQUALS(expectedCells, cells);
    const VecPt3d& actualPoints = ugrid->GetLocations();
    VecPt3d expectedPoints = points;
    expectedPoints.push_back({-10, 13.3333, 0.0});
    expectedPoints.push_back({-20.0, 33.3333, 0.0});
    TS_ASSERT_DELTA_VECPT3D(expectedPoints, actualPoints, 1.0e-3);
  }
} // MeQuadBlossomUnitTests::testComplexQuad
//------------------------------------------------------------------------------
/// \brief Test MeQuadBlossom::SplitToQuads.
//------------------------------------------------------------------------------
void MeQuadBlossomUnitTests::testSplitToQuads()
{
  // clang-format off
  //! [snip_test_2DShapes]
  VecPt3d points = { { 0, 0, 0 }, { 10, 0, 0 }, { 20, 0, 0 }, { 30, 0, 0 },
                     {0, 10, 0},  {10, 10, 0}, {20, 10, 0}, {30, 10, 0},
                     {0, 20, 0}, {10, 20, 0}, {20, 20, 0}, {30, 20, 0}
  };

  // Cell type (5), number of points (3), point numbers, counterclockwise
  VecInt cells = {
    // row 1
    XMU_TRIANGLE, 3, 0, 5, 4,
    XMU_TRIANGLE, 3, 0, 1, 5,
    XMU_TRIANGLE, 3, 1, 2, 5,
    XMU_TRIANGLE, 3, 2, 6, 5,
    XMU_QUAD, 4, 2, 3, 7, 6,
    // row 2
    XMU_TRIANGLE, 3, 4, 5, 8,
    XMU_POLYGON, 3, 5, 9, 8,
    XMU_QUAD, 4, 5, 6, 10, 9,
    XMU_POLYGON, 4, 6, 7, 11, 10
  };
  //! [snip_test_2DShapes]
  // clang-format on

  BSHP<XmUGrid> ugrid = XmUGrid::New(points, cells);

  BSHP<XmUGrid> quadUGrid = MeQuadBlossom::SplitToQuads(ugrid);

  VecPt3d quadPoints = quadUGrid->GetLocations();
  // clang-format off
  VecPt3d expectedQuadPoints = {
    { 0, 0, 0 }, { 10, 0, 0 }, { 20, 0, 0 }, { 30, 0, 0 },
    {0, 10, 0},  {10, 10, 0}, {20, 10, 0}, {30, 10, 0},
    {0, 20, 0}, {10, 20, 0}, {20, 20, 0}, {30, 20, 0},
    // The split edges
    {5, 0, 0}, {0, 5, 0}, {5, 5, 0}, {15, 0, 0}, {10, 5, 0}, // 12-16
    {25, 0, 0}, {15, 5, 0}, {20, 5, 0}, {30, 5, 0}, {5, 10, 0}, // 17-21
    {0, 15, 0}, {15, 10, 0}, {5, 15, 0}, {10, 15, 0}, {25, 10, 0}, // 22-26
    {20, 15, 0}, {30, 15, 0}, {5, 20, 0}, {15, 20, 0}, {25, 20, 0}, // 27-31
    // The centroids
    {3.33333, 6.66667, 0}, {6.66667, 3.33333, 0}, // 32-33
    {13.33333, 3.33333, 0}, {16.66667, 6.66667, 0}, // 34-35
    {25, 5, 0}, {3.33333, 13.33333, 0}, {6.66667, 16.66667, 0}, // 36-38
    {15, 15, 0}, {25, 15, 0}  //39-40
  };
  TS_ASSERT_DELTA_VECPT3D(expectedQuadPoints, quadPoints, 1.0e-5);

  VecInt expectedQuadCells = {
    XMU_QUAD, 4, 0, 14, 32, 13,
    XMU_QUAD, 4, 5, 21, 32, 14,
    XMU_QUAD, 4, 4, 13, 32, 21,
    XMU_QUAD, 4, 0, 12, 33, 14,
    XMU_QUAD, 4, 1, 16, 33, 12,
    XMU_QUAD, 4, 5, 14, 33, 16,
    XMU_QUAD, 4, 1, 15, 34, 16,
    XMU_QUAD, 4, 2, 18, 34, 15,
    XMU_QUAD, 4, 5, 16, 34, 18,
    XMU_QUAD, 4, 2, 19, 35, 18,
    XMU_QUAD, 4, 6, 23, 35, 19,
    XMU_QUAD, 4, 5, 18, 35, 23,
    XMU_QUAD, 4, 2, 17, 36, 19,
    XMU_QUAD, 4, 3, 20, 36, 17,
    XMU_QUAD, 4, 7, 26, 36, 20,
    XMU_QUAD, 4, 6, 19, 36, 26,
    XMU_QUAD, 4, 4, 21, 37, 22,
    XMU_QUAD, 4, 5, 24, 37, 21,
    XMU_QUAD, 4, 8, 22, 37, 24,
    XMU_QUAD, 4, 5, 25, 38, 24,
    XMU_QUAD, 4, 9, 29, 38, 25,
    XMU_QUAD, 4, 8, 24, 38, 29,
    XMU_QUAD, 4, 5, 23, 39, 25,
    XMU_QUAD, 4, 6, 27, 39, 23,
    XMU_QUAD, 4, 10, 30, 39, 27,
    XMU_QUAD, 4, 9, 25, 39, 30,
    XMU_QUAD, 4, 6, 26, 40, 27,
    XMU_QUAD, 4, 7, 28, 40, 26,
    XMU_QUAD, 4, 11, 31, 40, 28,
    XMU_QUAD, 4, 10, 27, 40, 31
  };
  // clang-format on
  VecInt quadCells = quadUGrid->GetCellstream();
  TS_ASSERT_EQUALS_VEC(expectedQuadCells, quadCells);
} // MeQuadBlossomUnitTests::testSplitToQuads
//------------------------------------------------------------------------------
/// \brief Test function to get estimated time to run Quad Blossom algorithm.
//------------------------------------------------------------------------------
void MeQuadBlossomUnitTests::testEstimatedRunTime()
{
  TS_ASSERT_DELTA(0.0085, MeQuadBlossom::EstimatedRunTimeInMinutes(500), 0.001);
  TS_ASSERT_DELTA(0.0680, MeQuadBlossom::EstimatedRunTimeInMinutes(1000), 0.001);
  TS_ASSERT_DELTA(8.5000, MeQuadBlossom::EstimatedRunTimeInMinutes(5000), 0.001);
  TS_ASSERT_DELTA(68.0, MeQuadBlossom::EstimatedRunTimeInMinutes(10000), 0.001);
  TS_ASSERT_DELTA(68000.0, MeQuadBlossom::EstimatedRunTimeInMinutes(100000), 0.001);
} // MeQuadBlossomUnitTests::testEstimatedRunTime
//------------------------------------------------------------------------------
/// \brief Test PreMakeQuads function.
//------------------------------------------------------------------------------
void MeQuadBlossomUnitTests::testPreMakeQuads()
{
  {
    // single triangle
    //  2-----1
    //  |    /
    //  |   /
    //  |  /
    //  | /
    //  0
    VecInt2d triangles = {{0, 1, 2}};
    VecPt3d points = {{0, 0, 0}, {10, 10, 0}, {0, 10, 0}};

    MeQuadBlossomImpl blossom(points, triangles);
    int numBoundaryEdges = blossom.PreMakeQuads();
    TS_ASSERT_EQUALS(3, numBoundaryEdges);

    VecEdge expect_interior = {};
    VecEdge expect_boundary = {
      {2, 0, 0, XM_NONE, 1, XM_NONE},
      {1, 2, 0, XM_NONE, 0, XM_NONE},
      {0, 1, 0, XM_NONE, 2, XM_NONE},
    };
    VecMeEdge expect_extra = {};
    VecInt2d expectExtraVertices = {};
    TS_ASSERT_EQUALS(expect_interior, blossom.m_interiorEdges);
    TS_ASSERT_EQUALS(expect_boundary, blossom.m_boundaryEdges);
    TS_ASSERT_EQUALS(expect_extra, blossom.m_extraEdges);
    TS_ASSERT_EQUALS(expectExtraVertices, blossom.m_extraPoints);
    
    bool splitBoundaryPoints = true;
    bool useAngle = false;
    blossom.MakeQuads(splitBoundaryPoints, useAngle);
    
    VecInt2d faces = blossom.m_faces;
    TS_ASSERT_EQUALS(triangles, faces);
  }

  {
    // two adjacent
    //  2-----3
    //  |    /|
    //  |   / |
    //  |  /  |
    //  | /   |
    //  0-----1
    VecInt2d triangles = {{2, 0, 3}, {3, 0, 1}};
    VecPt3d points = {{0, 0, 0}, {10, 0, 0}, {0, 10, 0}, {10, 10, 0}};

    MeQuadBlossomImpl blossom(points, triangles);
    int numBoundaryEdges = blossom.PreMakeQuads();
    TS_ASSERT_EQUALS(4, numBoundaryEdges);

    VecEdge expect_interior = {{0, 3, 0, 1, 2, 1}};
    VecEdge expect_boundary = {{3, 2, 0, XM_NONE, 0, XM_NONE},
                               {2, 0, 0, XM_NONE, 3, XM_NONE},
                               {1, 3, 1, XM_NONE, 0, XM_NONE},
                               {0, 1, 1, XM_NONE, 3, XM_NONE}
    };
    VecMeEdge expect_extra = {};
    VecInt2d expectExtraVertices = {};
    TS_ASSERT_EQUALS(expect_interior, blossom.m_interiorEdges);
    TS_ASSERT_EQUALS(expect_boundary, blossom.m_boundaryEdges);
    TS_ASSERT_EQUALS(expect_extra, blossom.m_extraEdges);
    TS_ASSERT_EQUALS(expectExtraVertices, blossom.m_extraPoints);
    
    bool splitBoundaryPoints = true;
    bool useAngle = false;
    blossom.MakeQuads(splitBoundaryPoints, useAngle);
    
    VecInt2d faces = blossom.m_faces;
    VecInt2d expectedFaces = {{0, 1, 3, 2}};
    TS_ASSERT_EQUALS(expectedFaces, faces);
  }

  {
    // three adjacent
    VecInt2d triangles = {{2, 9, 6}, {9, 4, 6}, {9, 3, 4}};
    VecPt3d points(10, Pt3d());
    points[2] = {0, 0, 0};
    points[3] = {10, 20, 0};
    points[4] = {0, 20, 0};
    points[6] = {0, 10, 0};
    points[9] = {10, 10, 0};
    

    MeQuadBlossomImpl blossom(points, triangles);
    int numBoundaryEdges = blossom.PreMakeQuads();
    TS_ASSERT_EQUALS(5, numBoundaryEdges);

    VecEdge expect_interior = {
      {6, 9, 1, 0, 4, 2},
      {4, 9, 2, 1, 3, 6},
    };
    VecEdge expect_boundary = {
      {9, 3, 2, XM_NONE, 4, XM_NONE}, {6, 2, 0, XM_NONE, 9, XM_NONE},
      {4, 6, 1, XM_NONE, 9, XM_NONE}, {3, 4, 2, XM_NONE, 9, XM_NONE},
      {2, 9, 0, XM_NONE, 6, XM_NONE},
    };
    VecMeEdge expect_extra = {{0, 2, -1000}};
    VecInt2d expectExtraVertices = {{9, 6, 2, 3, 4}};
    TS_ASSERT_EQUALS(expect_interior, blossom.m_interiorEdges);
    TS_ASSERT_EQUALS(expect_boundary, blossom.m_boundaryEdges);
    TS_ASSERT_EQUALS(expect_extra, blossom.m_extraEdges);
    TS_ASSERT_EQUALS(expectExtraVertices, blossom.m_extraPoints);
    
    bool splitBoundaryPoints = true;
    bool useAngle = false;
    blossom.MakeQuads(splitBoundaryPoints, useAngle);
    
    VecInt2d faces = blossom.m_faces;
    VecInt2d expectedFaces = {{2, 9, 6}, {4, 6, 9, 3}};
    TS_ASSERT_EQUALS(expectedFaces, faces);
  }
} // MeQuadBlossomUnitTests::testPreMakeQuads

#endif // CXX_TEST
