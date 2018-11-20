//------------------------------------------------------------------------------
/// \file
/// \ingroup meshing_detail
/// \copyright (C) Copyright Aquaveo 2018. Distributed under FreeBSD License
/// (See accompanying file LICENSE or https://aqaveo.com/bsd/license.txt)
//------------------------------------------------------------------------------
#pragma once

//----- Included files ---------------------------------------------------------
#include <list>
#include <set>
#include <vector>
#include <xmscore/points/pt.h>
#include <xmscore/misc/boost_defines.h>

//----- Forward declarations ---------------------------------------------------

//----- Namespace declaration --------------------------------------------------
namespace xms
{
// typedef boost::unordered_set<size_t> setIdx;
typedef std::set<size_t> SetIdx; ///< typedef for shorter declaration

//----- Constants / Enumerations -----------------------------------------------

//----- Structs / Classes ------------------------------------------------------
/// \brief Utility class to work with polygon paving
class MePolyPts
{
public:
  MePolyPts();

  double& XyTol();
  // Pt3d &Offset();
  std::vector<Pt3d>& Pts();
  BSHP<std::vector<Pt3d>> PtsSharedPointer();

  std::vector<size_t> HashPts();
  std::vector<size_t> SegmentsForCleanPolyOffset();
  void IntersectSegs(const std::vector<size_t>& a_segs);
  std::list<size_t> SequenceWithIntersects(std::vector<size_t>& a_segs);
  void CheckIntersectTwoSegs(size_t a_i,
                             size_t a_j,
                             const std::vector<size_t>& a_iSeg,
                             const std::vector<size_t>& a_jSeg);
  void CalcLoopsForCleanPolyOffset(std::list<size_t>& a_sequence,
                                   std::list<std::vector<size_t>>& a_loops);
  void RemoveBackwardLoopsForCleanPolyOffset(std::list<std::vector<size_t>>& a_loops, int a_pType);
  void ClassifyLoopsFromInPolyAndRemoveInvalid(std::list<std::vector<size_t>>& a_loops,
                                               std::vector<int>& a_loopType);

  bool PolyInsideOfPoly(const std::vector<size_t>& a_poly, const std::vector<size_t>& a_polyToTest);
  bool PolyInsideOfPoly(const std::vector<Pt3d>& a_poly, const std::vector<size_t>& a_polyToTest);
  size_t IdxFromPt3d(const Pt3d& a_pt);

private:
  class impl;
  BSHP<impl> m_p; ///< implementation class
};

//----- Function prototypes ----------------------------------------------------

} // namespace xms
