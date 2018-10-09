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
#include <xmsmesh/meshing/MeBadQuadRemover.h>

// 3. Standard library headers
#include <cmath>
#include <numeric>

// 4. External library headers
#include <boost/utility.hpp>

// 5. Shared code headers
#include <xmscore/misc/DynBitset.h>
#include <xmscore/misc/XmError.h>
#include <xmscore/misc/xmstype.h>
#include <xmsgrid/ugrid/XmEdge.h>
#include <xmsgrid/ugrid/XmUGrid.h>
#include <xmsinterp/geometry/geoms.h>

// 6. Non-shared code headers

//----- Forward declarations ---------------------------------------------------

//----- External globals -------------------------------------------------------

//----- Namespace declaration --------------------------------------------------

namespace xms
{
//----- Constants / Enumerations -----------------------------------------------

//----- Classes / Structs ------------------------------------------------------

namespace
{
/// Contains the number of 3 edge points in the cell and index of the last three
/// edge point if any.
struct CellData
{
  /// Constructor.
  /// \param[in] a_num3EdgePoints The number of cell points that have 3 adjacent
  /// edges.
  /// \param[in] a_pointIdx The index of the last three edge point if any.
  explicit CellData(int a_num3EdgePoints = -999, int a_pointIdx = -1)
  : m_num3EdgePoints(a_num3EdgePoints)
  , m_pointIdx(a_pointIdx)
  {
  }

  int m_num3EdgePoints; ///< Count of non-boundary vertices in the face with 3 attached edges. -999
                        ///< means it hasn't been set. Less than 0 if face is not a quad. -4 is used
                        ///< to flag the cell as not collapsible.
  int m_pointIdx; ///< The index of the point in a quad that will collapse to its diagonal point if
                  ///< the quad is collapsible (or -1 if the collapse criteria is not met or not yet
                  ///< determined).
};

typedef std::vector<CellData> VecCellData; ///< Vector of CellData.

class MeBadQuadRemoverImpl : public MeBadQuadRemover
{
public:
  MeBadQuadRemoverImpl(BSHP<XmUGrid> a_ugrid);

  virtual BSHP<XmUGrid> RemoveBadQuads(double a_maxAspect = 0.7) override;

  // implementation helpers
  bool ReplacePoint(int a_ptIdx, int a_newPtIdx);
  void MovePoint(int a_ptIdx, const Pt3d& a_newPoint);
  void DeleteCell(int a_cellIdx);
  BSHP<XmUGrid> BuildUGridFromReplacedPoints();
  void CollapseFromPoint(int a_cellIdx, int a_pointIdx_w3, const VecInt& a_adjCells);
  void ComputeCellData(int a_cellIdx, double max_aspect);
  bool CanCollapse(int a_cellIdx, int pointIdx_w3, VecInt& a_adjCells);

private:
  typedef std::pair<int, Pt3d> MovedPoint;       ///< Index of point and new location.
  typedef std::vector<MovedPoint> MovedPointVec; ///< Vector of moved points.
  BSHP<XmUGrid> m_ugrid;                         ///< The input UGrid containing the bad quads.
  VecInt m_pointIdxMap;        ///< A mapping of original index to new point index.
  DynBitset m_cellsToDelete;   ///< True if a cell is to be deleted.
  MovedPointVec m_movedPoints; ///< List of points moved with new location.
  VecCellData m_cellsData;     ///< Vector of potentially collapsable quads.

  /// Number of adjacent points to m_ugrid->GetPoint(i). Negative if a boundary
  /// point.
  const VecInt m_adjPointCnts;
};

//----- Internal functions -----------------------------------------------------

//----- Class / Function definitions -------------------------------------------

//------------------------------------------------------------------------------
/// \brief Build a UGrid from a vector of points and a vector of faces.
/// \param[in] a_points The UGrid points.
/// \param[in] a_faces A vector of vectors of counter clock-wise oriented point
/// indices. Must be 3 or 4 indices.
/// \return The UGrid created from the points and faces.
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
/// \brief Get point indices attached to a UGrid point across edges.
/// \param[in] a_ugrid The UGrid.
/// \param[in] a_pointIdx The point index to get the adjacent points from.
/// \param[out] a_edgePoints The point indices that are adjacent to the point.
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
} // GetPointIdxsAttachedByEdge
//------------------------------------------------------------------------------
/// \brief Get The count of adjacent points to a given point. The result is
/// the negative count if the point is on a boundary.
/// \param[in] a_ugrid The UGrid containing the point.
/// \param[in] a_pointIdx The index of the point.
/// \return The number of edges eminating from each vertex (but negative if the
/// vertex is on the boundary.
//------------------------------------------------------------------------------
int GetAdjacentPointCount(BSHP<XmUGrid> a_ugrid, int a_pointIdx)
{
  VecInt adjacentPoints;
  GetPointIdxsAttachedByEdge(a_ugrid, a_pointIdx, adjacentPoints);
  bool isBoundary = false;
  VecInt adjacentCells;
  for (auto adjacentPtIdx : adjacentPoints)
  {
    a_ugrid->GetPointsAdjacentCells(a_pointIdx, adjacentPtIdx, adjacentCells);
    if (adjacentCells.size() == 1)
    {
      isBoundary = true;
      break;
    }
  }

  int count = (int)adjacentPoints.size();
  return isBoundary ? -count : count;
} // GetAdjacentPointCount
//------------------------------------------------------------------------------
/// \brief Get the adjacent point counts for the points in a UGrid.
/// \param[in] a_ugrid The UGrid to get the count for.
/// \return The number of edges eminating from each point (but negative if the
/// point is on the boundary.
//------------------------------------------------------------------------------
VecInt GetAdjacentPointCounts(BSHP<XmUGrid> a_ugrid)
{
  VecInt counts(a_ugrid->PointCount());
  int numPoints = a_ugrid->PointCount();
  for (int pointIdx = 0; pointIdx < numPoints; ++pointIdx)
  {
    counts[pointIdx] = GetAdjacentPointCount(a_ugrid, pointIdx);
  }
  return counts;
} // GetAdjacentPointCounts
////////////////////////////////////////////////////////////////////////////////
/// \class MeBadQuadRemoverImpl
/// \brief Identifies and removes badly formed quads. Badly formed quads
/// include those
/// - that have a non-boundary point that is adjacent to only 2 cells
/// - that have only one non-boundary point that is adjacent to 3 cells
/// - that have exactly two non-boundary points opposite each other that
///   both have exactly 3 adjacent cells
/// and that have a narrow aspect ratio of the diagonal including the 2 or 3
/// adjacent cell points and the other diagonal.
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief Constructor.
/// \param[in] a_ugrid The UGrid to remove badly formed quads from.
//------------------------------------------------------------------------------
MeBadQuadRemoverImpl::MeBadQuadRemoverImpl(BSHP<XmUGrid> a_ugrid)
: m_ugrid(a_ugrid)
, m_pointIdxMap(a_ugrid->PointCount(), -1)
, m_adjPointCnts(GetAdjacentPointCounts(a_ugrid))
{
  int cellCount = a_ugrid->GetCellCount();
  m_cellsToDelete.resize(cellCount);
  m_cellsData.resize(cellCount);
} // MeBadQuadRemoverImpl::MeBadQuadRemoverImpl
//------------------------------------------------------------------------------
/// \brief Remove bad quads and return a reconstructed UGrid with them removed.
/// \param[in] a_maxAspect The maximum aspect ratio for the diagonals.
/// \return The reconstructed UGrid with the bad quads removed.
//------------------------------------------------------------------------------
BSHP<XmUGrid> MeBadQuadRemoverImpl::RemoveBadQuads(double a_maxAspect)
{
  if (a_maxAspect != 0.0)
  {
    a_maxAspect *= a_maxAspect;
  }

  int cellCnt = m_ugrid->GetCellCount();

  for (int a_cellIdx = 0; a_cellIdx < cellCnt; ++a_cellIdx)
  {
    CellData& cellData = m_cellsData[a_cellIdx];
    if (cellData.m_num3EdgePoints == -999)
    {
      ComputeCellData(a_cellIdx, a_maxAspect);
    }
  }

  while (true)
  {
    int collapseCnt = 0;
    for (int a_cellIdx = 0; a_cellIdx < cellCnt; ++a_cellIdx)
    {
      if (!m_cellsToDelete[a_cellIdx])
      {
        CellData& data = m_cellsData[a_cellIdx];
        if (data.m_num3EdgePoints == 1)
        {
          int collapsablePtIdx = data.m_pointIdx;
          VecInt adjCells;
          if (collapsablePtIdx != -1 && CanCollapse(a_cellIdx, collapsablePtIdx, adjCells))
          {
            ++collapseCnt;
            CollapseFromPoint(a_cellIdx, collapsablePtIdx, adjCells);
          }
          // m_cellsData[a_cellIdx].m_num3EdgeVertices = -4;
        }
      }
    }
    if (collapseCnt == 0)
    {
      break;
    }
  }

  BSHP<XmUGrid> newUgrid = BuildUGridFromReplacedPoints();
  return newUgrid;
} // MeBadQuadRemoverImpl::RemoveBadQuads
//------------------------------------------------------------------------------
/// \brief Create a mapping to redirect all references from one point index to
/// another in the new UGrid.
/// \param[in] a_ptIdx The point to replace.
/// \param[in] a_newPtIdx The point to replace a_ptIdx with.
/// \return true if the neither point has already been replaced.
//------------------------------------------------------------------------------
bool MeBadQuadRemoverImpl::ReplacePoint(int a_ptIdx, int a_newPtIdx)
{
  int& oldIndex = m_pointIdxMap[a_ptIdx];
  int& newIndex = m_pointIdxMap[a_newPtIdx];
  if (oldIndex == -1 && newIndex == -1)
  {
    m_pointIdxMap[a_ptIdx] = a_newPtIdx;
    return true;
  }
  return false;
} // ReplacePoint
//------------------------------------------------------------------------------
/// \brief Record an action to move a point to a new location in the new UGrid.
/// \param[in] a_ptIdx The point to move.
/// \param[in] a_newPoint The new point location.
//------------------------------------------------------------------------------
void MeBadQuadRemoverImpl::MovePoint(int a_ptIdx, const Pt3d& a_newPoint)
{
  m_movedPoints.push_back({a_ptIdx, a_newPoint});
} // MeBadQuadRemoverImpl::MovePoint
//------------------------------------------------------------------------------
/// \brief Record an action to delete a cell from the old UGrid when creating
/// the new UGrid.
/// \param[in] a_cellIdx The index of the cell to delete.
//------------------------------------------------------------------------------
void MeBadQuadRemoverImpl::DeleteCell(int a_cellIdx)
{
  m_cellsToDelete[a_cellIdx] = true;
} // MeBadQuadRemoverImpl::DeleteCell
//------------------------------------------------------------------------------
/// \brief Build a ugrid given the recorded actions to take including deleting
/// cells, moving points, and replacing points.
/// \return The new UGrid with the actions performed.
//------------------------------------------------------------------------------
BSHP<XmUGrid> MeBadQuadRemoverImpl::BuildUGridFromReplacedPoints()
{
  VecPt3d oldPoints = m_ugrid->GetLocations();
  for (auto& movedPoint : m_movedPoints)
  {
    oldPoints[movedPoint.first] = movedPoint.second;
  }

  VecPt3d newPoints;
  int numPoints = (int)m_ugrid->PointCount();
  newPoints.reserve(numPoints);
  int currPtIdx = 0;
  VecInt newPointIdxLookupTable(numPoints, -1);
  // - - - - - - - before mapping any
  // 1 - 5 - - - 3 lookup after mapping 0 to 1; 2 to 5; 6 to 3
  for (int pointIdx = 0; pointIdx < numPoints; ++pointIdx)
  {
    if (m_pointIdxMap[pointIdx] == -1)
    {
      newPoints.push_back(oldPoints[pointIdx]);
      newPointIdxLookupTable[pointIdx] = currPtIdx;
      ++currPtIdx;
    }
  }
  // - 0 - 1 2 3 - newPointIdxLookupTable after the first pass

  for (int pointIdx = 0; pointIdx < numPoints; ++pointIdx)
  {
    int idx = m_pointIdxMap[pointIdx];
    if (idx != -1)
    {
      newPointIdxLookupTable[pointIdx] = newPointIdxLookupTable[idx];
    }
  }
  // 0 * 3 * * * 1 newPointIdxLookupTable after second pass (* means unchanged)
  // 0 0 3 1 2 3 1 newPointIdxLookupTable after second pass
  VecInt cells;
  int numCells = m_ugrid->GetCellCount();
  cells.reserve(6 * numCells);
  for (int cellIdx = 0; cellIdx < numCells; ++cellIdx)
  {
    if (!m_cellsToDelete[cellIdx])
    {
      int cellType = m_ugrid->GetCellType(cellIdx);
      XM_ASSERT(m_ugrid->GetCellDimension(cellIdx) == 2);
      cells.push_back(cellType);
      cells.push_back(m_ugrid->GetCellEdgeCount(cellIdx));
      VecInt cellPoints = m_ugrid->GetCellPoints(cellIdx);
      for (auto pointIdx : cellPoints)
      {
        cells.push_back(newPointIdxLookupTable[pointIdx]);
      }
    }
  }
  BSHP<XmUGrid> newUGrid = XmUGrid::New(newPoints, cells);
  return newUGrid;
} // MeBadQuadRemoverImpl::BuildUGridFromReplacedPoints
//------------------------------------------------------------------------------
/// \brief Collapses one cell on itself by moving one point to its opposite
/// diagonal point.
/// \param[in] a_cellIdx The index of the cell to collapse.
/// \param[in] a_pointIdx_w3 The UGrid index of the point to collapse.
/// \param[in] a_adjCells The other cells adjacent to a_pointIdx_w3.
//------------------------------------------------------------------------------
void MeBadQuadRemoverImpl::CollapseFromPoint(int a_cellIdx,
                                             int a_pointIdx_w3,
                                             const VecInt& a_adjCells)
{
  VecInt pointIdxs = m_ugrid->GetCellPoints(a_cellIdx);
  auto it = std::find(pointIdxs.begin(), pointIdxs.end(), a_pointIdx_w3);
  XM_ASSERT(it != pointIdxs.end());
  int position = int(it - pointIdxs.begin());
  int diagonalPtIdx = pointIdxs[(position + 2) % (int)pointIdxs.size()];
  if (ReplacePoint(a_pointIdx_w3, diagonalPtIdx))
  {
    int diagonalCnt = m_adjPointCnts[diagonalPtIdx];
    int vCnt = m_adjPointCnts[a_pointIdx_w3];
    XM_ASSERT(vCnt == 3);
    if (diagonalCnt > 0)
    {
      const VecPt3d& points = m_ugrid->GetLocations();
      double sumWt = (double)(diagonalCnt + vCnt);
      double ptWt = double(vCnt) / sumWt;
      double diagonalPtWt = double(diagonalCnt) / sumWt;
      Pt3d movedPos = diagonalPtWt * points[diagonalPtIdx] + ptWt * points[a_pointIdx_w3];
      MovePoint(diagonalPtIdx, movedPos);
    }
    DeleteCell(a_cellIdx);
    for (auto adjCell : a_adjCells)
    {
      m_cellsData[adjCell].m_num3EdgePoints = -4;
    }
  }
} // MeBadQuadRemoverImpl::CollapseFromPoint
//------------------------------------------------------------------------------
/// \brief Determines if a cell is badly shaped and can be collapsed.
///
/// If the cell has points only 2 adjacent cells, record the actions to collapse
/// it. It if has exactly 2 points that are opposite each other with 3 adjacent
/// cells, record the actions to collapse it.  If it has exactly 1 three-cell
/// point, record that for possible future collapse. In all other cases, flag
/// the cell as not collapsable.  If a cell is collapsed, flag the cells that
/// are adjacent to its collapsed point as not collapsable.
/// \param[in] a_cellIdx The index of the cell to check.
/// \param[in] a_maxAspect The maximum aspect ratio for the diagonals.
//------------------------------------------------------------------------------
void MeBadQuadRemoverImpl::ComputeCellData(int a_cellIdx, double max_aspect)
{
  VecInt pointIdxs = m_ugrid->GetCellPoints(a_cellIdx);
  VecPt3d points = m_ugrid->GetPointsLocations(pointIdxs);
  int adjPointCnt = (int)points.size();

  if (adjPointCnt != 4)
  {
    CellData data(-adjPointCnt, -1);
    m_cellsData[a_cellIdx] = data;
    return;
  }

  // Compute the aspect ratio of diagonals
  VecBool aspectOk = {true, true, true, true};
  VecDbl diagonals(2, 0.0);
  if (max_aspect > 0.0)
  {
    double d0 = gmXyDistanceSquared(points[0], points[2]);
    double d1 = gmXyDistanceSquared(points[1], points[3]);
    diagonals[0] = d0;
    diagonals[1] = d1;
    bool even = d0 / d1 <= max_aspect;
    bool odd = d1 / d0 <= max_aspect;
    aspectOk = {even, odd, even, odd};
  }

  int bits = 0;
  int threes = 0;
  int pointIdx_w3 = XM_NONE;
  for (int i = 0; i < 4; ++i)
  {
    int pointIdx = pointIdxs[i];
    int adjPointCnt = m_adjPointCnts[pointIdx];
    bool on_boundary = adjPointCnt < 0;
    adjPointCnt = abs(adjPointCnt);
    if (adjPointCnt == 2 && !on_boundary)
    {
      int opposingIdx = pointIdxs[(i + 2) % 4];
      Pt3d opposingPt = m_ugrid->GetPointLocation(opposingIdx);
      int adjCellIdx = m_ugrid->GetCell2dEdgeAdjacentCell(a_cellIdx, i);
      VecInt adjPointIdxs = m_ugrid->GetCellPoints(adjCellIdx);
      if (adjPointIdxs.size() == 3)
      {
        // adjacent cell is a triangle
        // +
        // |\ \
        // | \   \
        // |  +    +
        // | /   /
        // |/ /
        // +
        if (ReplacePoint(opposingIdx, pointIdxs[i]))
        {
          MovePoint(pointIdxs[i], opposingPt);
          DeleteCell(a_cellIdx);
        }
        m_cellsData[a_cellIdx].m_num3EdgePoints = -4;
        m_cellsData[adjCellIdx].m_num3EdgePoints = -4;
        return;
      }
      else if (adjPointIdxs.size() == 4)
      {
        // Compute the other diagonal if the edges were deleted
        // adjOpposingPt ---> +------+
        //                    |     /|
        //                    |    / |
        //                    |   + <--- at this point and cell to right
        //                    |  /   |
        //                    | /    |
        //                    +------+ <--- opposingPt

        auto it = std::find(adjPointIdxs.begin(), adjPointIdxs.end(), pointIdx);
        int position = int(it - adjPointIdxs.begin());
        int adjOpposingIdx = adjPointIdxs[(position + 2) % (int)adjPointIdxs.size()];
        // int adjOpposingIdx = pointIdxs[(position + 2) % (int)adjPointIdxs.size()];
        Pt3d adjOpposingPt = m_ugrid->GetPointLocation(adjOpposingIdx);
        // Pt3d adjOpposingPt = m_ugrid->GetPoint(adjPointIdxs.at(adjOpposingIdx));
        double d0 = gmXyDistanceSquared(adjOpposingPt, opposingPt);
        double d1 = diagonals[(i + 1) & 0x1];
        if (d0 / d1 < 1.0)
        {
          if (ReplacePoint(pointIdxs[i], opposingIdx))
          {
            DeleteCell(a_cellIdx);
            m_cellsData[adjCellIdx].m_num3EdgePoints = -4;
            m_cellsData[a_cellIdx].m_num3EdgePoints = -4;
          }
          return;
        }
      }
    }

    if (adjPointCnt == 2 && on_boundary)
    {
      threes += 1;
      if (aspectOk[i])
      {
        bits |= 0x1;
      }
    }

    if (adjPointCnt == 3 && !on_boundary)
    {
      threes += 1;
      if (aspectOk[i])
      {
        bits |= 0x1;
        pointIdx_w3 = pointIdx;
      }
    }

    bits <<= 1;
  } // end for

  bits >>= 1;

  if (pointIdx_w3 != XM_NONE && threes == 2 &&
      (bits == BOOST_BINARY(0101) || bits == BOOST_BINARY(1010)))
  {
    VecInt adjCells = m_ugrid->GetPointAdjacentCells(pointIdx_w3);
    CollapseFromPoint(a_cellIdx, pointIdx_w3, adjCells);
    return;
  }

  m_cellsData[a_cellIdx] = CellData(threes, pointIdx_w3);
} // MeBadQuadRemoverImpl::Compute3VertexCellData
//------------------------------------------------------------------------------
/// \brief Determines a cell can collapse from a given point.
/// \param[in] a_cellIdx The index of the cell to check.
/// \param[in] a_pointIdx_w3 The index of the point with 3 adjacent cells.
/// \param[out] a_adjCells The other two cells adjacent to point.
//------------------------------------------------------------------------------
bool MeBadQuadRemoverImpl::CanCollapse(int a_cellIdx, int pointIdx_w3, VecInt& a_adjCells)
{
  int adjPointCnt = m_adjPointCnts[pointIdx_w3];
  bool on_boundary = adjPointCnt < 0;
  if (on_boundary)
  {
    return false;
  }

  adjPointCnt = std::abs(adjPointCnt);
  XM_ASSERT(adjPointCnt == 3 && !on_boundary);

  a_adjCells = m_ugrid->GetPointAdjacentCells(pointIdx_w3);
  for (auto adjCell : a_adjCells)
  {
    if (adjCell != a_cellIdx)
    {
      const CellData& fd = m_cellsData[adjCell];
      int fd0 = fd.m_num3EdgePoints;
      if (fd0 < 0 || fd0 > 1)
      {
        return false;
      }
    }
  }

  return true;
} // MeBadQuadRemoverImpl::CanCollapse

} // namespace

////////////////////////////////////////////////////////////////////////////////
/// \class MeBadQuadRemover
/// \see MeBadQuadRemoverImpl
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief Create new MeBadQuadRemover.
/// \param[in] a_ugrid The UGrid to remove badly formed quads from.
/// \return The new MeBadQuadRemover.
//------------------------------------------------------------------------------
BSHP<MeBadQuadRemover> MeBadQuadRemover::New(BSHP<XmUGrid> a_ugrid)
{
  BSHP<MeBadQuadRemover> badQuadRemover(new MeBadQuadRemoverImpl(a_ugrid));
  return badQuadRemover;
} // MeBadQuadRemover::New
//------------------------------------------------------------------------------
/// \brief Constructor.
//------------------------------------------------------------------------------
MeBadQuadRemover::MeBadQuadRemover()
{
} // MeBadQuadRemover::MeBadQuadRemover
//------------------------------------------------------------------------------
/// \brief Destructor.
//------------------------------------------------------------------------------
MeBadQuadRemover::~MeBadQuadRemover()
{
} // MeBadQuadRemover::~MeBadQuadRemover

} // namespace xms

#if CXX_TEST
////////////////////////////////////////////////////////////////////////////////
// UNIT TESTS
////////////////////////////////////////////////////////////////////////////////

#include <xmsmesh/meshing/MeBadQuadRemover.t.h>

#include <xmscore/testing/TestTools.h>

//----- Namespace declaration --------------------------------------------------

using namespace xms;

////////////////////////////////////////////////////////////////////////////////
/// \class MeBadQuadRemoverUnitTests
/// \brief Tests for MeBadQuadRemover class.
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief Test GetAdjacentPointCounts.
//------------------------------------------------------------------------------
void MeBadQuadRemoverUnitTests::testGetAdjacentPointCounts()
{
  //     0-----1-----2
  //     |  0  |  1  |
  //     3-----4-----5
  //     |  2  |  3  |
  //     6-----7-----8

  BSHP<xms::XmUGrid> grid = TEST_XmUGridSimpleQuad();
  VecInt counts = GetAdjacentPointCounts(grid);
  VecInt expectedCounts = {-2, -3, -2, -3, 4, -3, -2, -3, -2};
  TS_ASSERT_EQUALS_VEC(expectedCounts, counts);
} // MeBadQuadRemoverUnitTests::testGetAdjacentPointCounts
//------------------------------------------------------------------------------
/// \brief Test ReplacePoint and BuildUGridFromReplacedPoints.
//------------------------------------------------------------------------------
void MeBadQuadRemoverUnitTests::testReplacePoints()
{
  {
    //   10------11----12----4
    //    |       |    /\    |
    //    |       |   /  \   |
    //    5-------6--7    8--9
    //    |      /|   \  /   |
    //    |  /13  |    \/    |
    //    0 ------1-----2----3
    VecPt3d points = {
      {0, 0, 0},   {10, 0, 0},  {20, 0, 0},  {30, 0, 0}, {30, 20, 0}, {0, 10, 0},  {10, 10, 0},
      {15, 10, 0}, {25, 10, 0}, {30, 10, 0}, {0, 20, 0}, {10, 20, 0}, {20, 20, 0}, {5, 3, 0},
    };
    VecInt cells = {XMU_QUAD, 4, 0, 1, 6,  13, XMU_QUAD, 4, 0, 13, 6,  5,
                    XMU_QUAD, 4, 1, 2, 7,  6,  XMU_QUAD, 4, 2, 3,  9,  8,
                    XMU_QUAD, 4, 5, 6, 11, 10, XMU_QUAD, 4, 6, 7,  12, 11,
                    XMU_QUAD, 4, 2, 8, 12, 7,  XMU_QUAD, 4, 8, 9,  4,  12};

    BSHP<XmUGrid> ugrid = XmUGrid::New(points, cells);
    MeBadQuadRemoverImpl replacer(ugrid);

    replacer.MovePoint(8, {20, 10, 0});

    TS_ASSERT(replacer.ReplacePoint(7, 8));
    replacer.DeleteCell(6);
    TS_ASSERT(replacer.ReplacePoint(13, 1));
    replacer.DeleteCell(0);

    VecPt3d expectPoints = {{0, 0, 0},   {10, 0, 0}, {20, 0, 0},  {30, 0, 0},
                            {30, 20, 0}, {0, 10, 0}, {10, 10, 0}, {20, 10, 0},
                            {30, 10, 0}, {0, 20, 0}, {10, 20, 0}, {20, 20, 0}};

    BSHP<XmUGrid> newUGrid = replacer.BuildUGridFromReplacedPoints();

    TS_ASSERT_DELTA_VECPT3D(expectPoints, newUGrid->GetLocations(), 1.0e-5);

    // 0 1 2 3 4 5 6 7  quad indices
    // x           x    flagged to delete
    VecInt expectCells = {// XMU_QUAD, 4, 0, 1, 6, 13,
                          XMU_QUAD, 4, 0, 1, 6, 5, XMU_QUAD, 4, 1, 2, 7, 6, XMU_QUAD, 4, 2, 3, 8, 7,
                          XMU_QUAD, 4, 5, 6, 10, 9, XMU_QUAD, 4, 6, 7, 11, 10,
                          // XMU_QUAD, 4, 2, 8, 12, 7,
                          XMU_QUAD, 4, 7, 8, 4, 11};
    //   10------11----12----4
    //    |       |    /\    |
    //    |       |   /  \   |
    //    5-------6--7*   8--9
    //    |      /|   \  /   |
    //    |  /13* |    \/    |
    //    0 ------1-----2----3

    //    9------10----11----4
    //    |       |     |    |
    //    |       |     |    |
    //    5-------6-----7----8
    //    |       |     |    |
    //    |       |     |    |
    //    0 ------1-----2----3

    TS_ASSERT_EQUALS(expectCells, newUGrid->GetCellStream());
  }
} // MeBadQuadRemoverUnitTests::testDeletePoints
//------------------------------------------------------------------------------
/// \brief Testing quad removal on a more complex quad UGrid. The original mesh
/// has two rows of quads with some 3 diamond shaped quads.  One if on the left
/// boundary.  One is roughly in the middle.  One is near the right side.  There
/// are also two pairs of quads where each pair is separated by 2 sequential
/// shared edges.  One has an aspect ratio to collapse.  The other does not.
/// Then toward the right is a quad that has 3 points with 3 edges attached. It
/// should not collapse.  One the bottom right is a quad with 4 3 edge points,
/// and it also should not collapse.
//------------------------------------------------------------------------------
void MeBadQuadRemoverUnitTests::testCollapse()
{
  /* The points below were generated in python in the following manner
    points = []
    for y in range(20, -10, -10):
        for x in range(0, 100, 10):
            points.append([x,y])

    points[10][0] -= 3
    points[11][0] -= 3
    points[12][0] -= 3
    points[13][0] -= 3

    points[14][0] += 3
    points[15][0] += 1.5

    points[17][0] -= 1.5
    points[18][0] -= 3
    points[19][0] += 5

    points.append([3, 10])    #30
    points.append([33, 10])   #31
    points.append([82, 10])    #32
    points.append([90, 12])  #33
    points.append([90, 8])  #34
    points.append([15, 15])  #35
    points.append([15, 5])  #36

    points.append([100, 10])  #37
    points.append([70, -10])  #38
    points.append([80, -10])  #39
    points.append([73, -3])  #40
    points.append([76, -4])  #41
    points.append([77, -7])  #42
    points.append([74, -6])  #43
  */

  VecInt2d faces = {{10, 20, 30, 0},  {30, 11, 1, 0},   {35, 12, 2, 1},   {11, 12, 35, 1},
                    {12, 13, 3, 2},   {31, 14, 4, 3},   {14, 15, 5, 4},   {15, 16, 6, 5},
                    {16, 17, 7, 6},   {17, 18, 8, 7},   {32, 33, 9, 8},   {20, 21, 11, 30},
                    {36, 22, 12, 11}, {21, 22, 36, 11}, {22, 23, 13, 12}, {13, 23, 31, 3},
                    {23, 24, 14, 31}, {24, 25, 15, 14}, {25, 26, 16, 15}, {26, 27, 17, 16},
                    {27, 28, 18, 17}, {18, 28, 32, 8},  {28, 29, 34, 32}, {34, 19, 33, 32},
                    {33, 19, 37, 9},  {34, 29, 37, 19}, {40, 41, 28, 27}, {41, 42, 39, 28},
                    {42, 43, 38, 39}, {38, 43, 40, 27}, {43, 42, 41, 40}};

  VecPt3d points = {// row 1 (0-9)
                    {0, 20, 0},
                    {10, 20, 0},
                    {20, 20, 0},
                    {30, 20, 0},
                    {40, 20, 0},
                    {50, 20, 0},
                    {60, 20, 0},
                    {70, 20, 0},
                    {80, 20, 0},
                    {90, 20, 0},
                    // row 2 (10-19)
                    {-3, 10, 0},
                    {7, 10, 0},
                    {17, 10, 0},
                    {27, 10, 0},
                    {43, 10, 0},
                    {51.5, 10, 0},
                    {60, 10, 0},
                    {68.5, 10, 0},
                    {77, 10, 0},
                    {95, 10, 0},
                    // row 3 (20-29)
                    {0, 0, 0},
                    {10, 0, 0},
                    {20, 0, 0},
                    {30, 0, 0},
                    {40, 0, 0},
                    {50, 0, 0},
                    {60, 0, 0},
                    {70, 0, 0},
                    {80, 0, 0},
                    {90, 0, 0},
                    // additional (30-43)
                    {3, 10, 0},
                    {33, 10, 0},
                    {82, 10, 0},
                    {90, 12, 0},
                    {90, 8, 0},
                    {15, 15, 0},
                    {15, 5, 0},
                    {100, 10, 0},
                    {70, -10, 0},
                    {80, -10, 0},
                    {73, -3, 0},
                    {76, -4, 0},
                    {77, -7, 0},
                    {74, -6, 0}};

  BSHP<XmUGrid> ugridIn = BuildUGrid(points, faces);
  BSHP<MeBadQuadRemover> remover = MeBadQuadRemover::New(ugridIn);
  BSHP<XmUGrid> collapsedUGrid = remover->RemoveBadQuads(0.7);
  VecPt3d expectedPoints = {// top row (0-9)
                            {0, 20, 0},
                            {10, 20, 0},
                            {20, 20, 0},
                            {30, 20, 0},
                            {40, 20, 0},
                            {50, 20, 0},
                            {60, 20, 0},
                            {70, 20, 0},
                            {80, 20, 0},
                            {90, 20, 0},
                            // row 2 (10-18)
                            {-3, 10, 0},
                            {7, 10, 0},
                            {17, 10, 0},
                            {30, 10, 0},
                            {43, 10, 0},
                            {51.5, 10, 0},
                            {60, 10, 0},
                            {68.5, 10, 0},
                            {95, 10, 0},
                            // bottom row (19-28)
                            {0, 0, 0},
                            {10, 0, 0},
                            {20, 0, 0},
                            {30, 0, 0},
                            {40, 0, 0},
                            {50, 0, 0},
                            {60, 0, 0},
                            {70, 0, 0},
                            {80, 0, 0},
                            {90, 0, 0},
                            // extra points in the middle row
                            {79.8571, 10, 0},
                            {90, 12, 0},
                            {90, 8, 0},
                            {15, 15, 0},
                            {100, 10, 0},
                            // extra points at the bottom
                            {70, -10, 0},
                            {80, -10, 0},
                            {73, -3, 0},
                            {76, -4, 0},
                            {77, -7, 0},
                            {74, -6, 0}};
  VecPt3d actualPoints = collapsedUGrid->GetLocations();
  TS_ASSERT_DELTA_VECPT3D(expectedPoints, actualPoints, 1.0e-4);

  VecInt expectedCells = {
    9, 4, 10, 11, 1,  0,  9, 4, 32, 12, 2,  1,  9, 4, 11, 12, 32, 1,  9, 4, 12, 13, 3,  2,
    9, 4, 13, 14, 4,  3,  9, 4, 14, 15, 5,  4,  9, 4, 15, 16, 6,  5,  9, 4, 16, 17, 7,  6,
    9, 4, 17, 29, 8,  7,  9, 4, 29, 30, 9,  8,  9, 4, 19, 20, 11, 10, 9, 4, 20, 21, 12, 11,
    9, 4, 21, 22, 13, 12, 9, 4, 22, 23, 14, 13, 9, 4, 23, 24, 15, 14, 9, 4, 24, 25, 16, 15,
    9, 4, 25, 26, 17, 16, 9, 4, 26, 27, 29, 17, 9, 4, 27, 28, 31, 29, 9, 4, 31, 18, 30, 29,
    9, 4, 30, 18, 33, 9,  9, 4, 31, 28, 33, 18, 9, 4, 36, 37, 27, 26, 9, 4, 37, 38, 35, 27,
    9, 4, 38, 39, 34, 35, 9, 4, 34, 39, 36, 26, 9, 4, 39, 38, 37, 36};
  VecInt actualCells = collapsedUGrid->GetCellStream();
  TS_ASSERT_EQUALS(expectedCells, actualCells);
} // MeBadQuadRemoverUnitTests::testCollapse
//------------------------------------------------------------------------------
/// \brief Test simple mesh with one quad and one triangle that share two
/// adjacent edges. Expected to result in one triangle.
//------------------------------------------------------------------------------
void MeBadQuadRemoverUnitTests::testCollapseQuadTri()
{
  VecInt2d faces = {{0, 2, 1, 3}, {0, 1, 2}};

  VecPt3d points = {
    {-10, 0, 0},
    {10, 0, 0}, // row 1
    {0, 10, 0}, // row 2
    {0, 20, 0}  // row 3
  };

  BSHP<XmUGrid> ugridIn = BuildUGrid(points, faces);
  BSHP<MeBadQuadRemover> remover = MeBadQuadRemover::New(ugridIn);
  BSHP<XmUGrid> collapsedUGrid = remover->RemoveBadQuads(0.7);
  VecPt3d expectedPoints = {{-10, 0, 0}, {10, 0, 0}, {0, 20, 0}};

  VecPt3d actualPoints = collapsedUGrid->GetLocations();
  TS_ASSERT_DELTA_VECPT3D(expectedPoints, actualPoints, 1.0e-4);

  VecInt expectedCells = {XMU_TRIANGLE, 3, 0, 1, 2};

  VecInt actualCells = collapsedUGrid->GetCellStream();
  TS_ASSERT_EQUALS(expectedCells, actualCells);
} // MeBadQuadRemoverUnitTests::testCollapseQuadTri

#endif // CXX_TEST
