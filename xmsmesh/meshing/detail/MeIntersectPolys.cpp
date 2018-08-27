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
#include <xmsmesh/meshing/detail/MeIntersectPolys.h>

// 3. Standard library headers

// 4. External library headers
#pragma warning(push)
#pragma warning(disable : 4512) // boost code: no assignment operator
#pragma warning(disable : 4244) // boost code: possible loss of data
#pragma warning(disable : 4127) // boost code: conditional expression is constant
#pragma warning(disable : 4267) // boost code: size_t to const int
#include <boost/geometry/geometry.hpp>
#pragma warning(pop)
#include <boost/unordered_set.hpp>

// 5. Shared code headers
#include <xmsinterp/geometry/GmBoostTypes.h> // GmBstPoly3d, XmBstRing
#include <xmsmesh/meshing/detail/MePolyOffsetter.h>
#include <xmscore/misc/XmError.h>

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
namespace bg = boost::geometry;

#define T_TOL 1e-13   ///< tolerance used in multipoly intersector
#define TOLERANCE 1e9 ///< tolerance used in PolyOffsetter
} // unnamed namespace

class MeIntersectPolys::impl
{
public:
  impl() {}

  void SetupInIn(const std::vector<MePolyOffsetterOutput>& a_offsets, double a_xyTol);
  void SetupInOut(const MePolyOffsetterOutput& a_offsets, double a_xyTol);
  void InInTrivialPolyCases();
  void InInDoIntersection();
  void InOutDoIntersection();
  void FillOutput(MePolyOffsetterOutput& a_out);
  void ClassifyPolys(const MePolyOffsetterOutput& a_input,
                     std::vector<std::vector<size_t>>& a_output);

  void OnSetup();
  void ClassifyDisjointPolys(SetIdx& a_delPolys);
  void ClassifyPolysInsideOfPolys(SetIdx& a_delPolys);
  void UpdateInPolyLoopsFromHash(const std::vector<size_t>& a_idx);
  void InOutCalcStartingStack();
  std::vector<size_t> InOutFindIntersectingPolys(size_t a_poly);
  void FillPolysInsideOfPolys(std::set<size_t>& a_oPoly,
                              std::set<size_t>& a_psetUsed,
                              std::vector<std::vector<size_t>>& a_output);

  GmBstPoly3d BoostPoly(size_t a_loopIdx);
  bool BoostPolyUnion(size_t a_i, size_t a_j);
  bool BoostPolySubtract(size_t a_i, size_t a_j);

  MePolyPts m_polyPts;                      ///< polygon points class
  std::vector<std::vector<size_t>> m_loops; ///< indexes to points that make polygons
  std::vector<int> m_loopTypes;        ///< type of polygon (OUTSIDE_POLY, INSIDE_POLY, NEWOUT_POLY)
  std::vector<GmBstBox3d> m_envelopes; ///< extents of polygons
  std::list<size_t> m_stack;           ///< stack for processing the polygons
  MePolyOffsetterOutput m_out;         ///< newpolygons created by theis class

  std::vector<GmBstPoly3d> m_bPolys;    ///< boost::geometry::polygon
  std::vector<GmBstPoly3d> m_bOutPolys; ///< boost::geometry::polygon
  std::vector<int> m_bOutPolyType; ///< type of polygon (OUTSIDE_POLY, INSIDE_POLY, NEWOUT_POLY)

  SetIdx m_inOutOpolys; ///< outside polygons

  // used by intersect IN wth IN
  boost::unordered_set<std::pair<size_t, size_t>>
    m_polyInsideOfPoly; ///< index of polygon that is inside of another polygon
};
//----- Internal functions -----------------------------------------------------

//----- Class / Function definitions -------------------------------------------
//------------------------------------------------------------------------------
/// \brief Returns the envelop of a polygon.
/// \param a_loop Vector of point indexes defining a polygon.
/// \param a_pts Vector of locations. The indexes in a_loop refer to this
/// vector.
/// \return Returns the envelop of a polygon.
//------------------------------------------------------------------------------
static GmBstBox3d iCalcPolyEnvelope(const std::vector<size_t>& a_loop,
                                    const std::vector<Pt3d>& a_pts)
{
  GmBstBox3d b;
  if (a_loop.empty() || a_pts.empty())
    return b;

  b.min_corner() = b.max_corner() = a_pts[a_loop[0]];
  b.min_corner().z = b.max_corner().z = 0;
  for (size_t j = 1; j < a_loop.size(); ++j)
  {
    const Pt3d& p(a_pts[a_loop[j]]);
    if (b.min_corner().x > p.x)
      b.min_corner().x = p.x;
    if (b.max_corner().x < p.x)
      b.max_corner().x = p.x;
    if (b.min_corner().y > p.y)
      b.min_corner().y = p.y;
    if (b.max_corner().y < p.y)
      b.max_corner().y = p.y;
  }
  return b;
} // iCalcPolyEnvelope
//------------------------------------------------------------------------------
/// \brief Returns true only if the envelopes overlap (not touch)
/// \param a_i Index to an envelope from a polygon.
/// \param a_j Index to an envelope from a polygon.
/// \param a_iEnv Vector of envelopes used by the a_i variable.
/// \param a_jEnv Vector of envelopes used by the a_j variable.
/// \return Returns true only if the envelopes overlap (not touch)
//------------------------------------------------------------------------------
static bool iEnvelopesOverlap(size_t a_i,
                              size_t a_j,
                              const std::vector<GmBstBox3d>& a_iEnv,
                              const std::vector<GmBstBox3d>& a_jEnv)
{
  const Pt3d &iMin(a_iEnv[a_i].min_corner()), &iMax(a_iEnv[a_i].max_corner()),
    &jMin(a_jEnv[a_j].min_corner()), &jMax(a_jEnv[a_j].max_corner());
  if (iMax.x > jMin.x && iMin.x < jMax.x && iMax.y > jMin.y && iMin.y < jMax.y)
    return true;
  return false;
} // iEnvelopesOverlap
//------------------------------------------------------------------------------
/// \brief Returns true if envelope a_i is inside of envelope a_j
/// \param a_i Index to an envelope from a polygon.
/// \param a_j Index to an envelope from a polygon.
/// \param a_envelopes Vector of envelopes used by the a_i & a_j variables.
/// \return Returns true if envelope a_i is inside of envelope a_j
//------------------------------------------------------------------------------
static bool iEnvelopeInsideOfEnvelope(size_t a_i,
                                      size_t a_j,
                                      const std::vector<GmBstBox3d>& a_envelopes)
{
  const Pt3d &iMin(a_envelopes[a_i].min_corner()), &iMax(a_envelopes[a_i].max_corner()),
    &jMin(a_envelopes[a_j].min_corner()), &jMax(a_envelopes[a_j].max_corner());
  if (jMin.x >= iMin.x && jMin.y >= iMin.y && jMax.x <= iMax.x && jMax.y <= iMax.y)
    return true;
  return false;
} // iEnvelopeInsideOfEnvelope

//------------------------------------------------------------------------------
/// \brief Constructor
//------------------------------------------------------------------------------
MeIntersectPolys::MeIntersectPolys()
: m_p(new MeIntersectPolys::impl())
{
} // MpInsidePolys::MpInsidePolys
//------------------------------------------------------------------------------
/// \brief Destructor
//------------------------------------------------------------------------------
MeIntersectPolys::~MeIntersectPolys()
{
  if (m_p)
    delete (m_p);
} // MpInsidePolys::MpInsidePolys
//------------------------------------------------------------------------------
/// \brief calls implementation
/// \param a_offsets: Vector of outputs from MePolyOffsetter::Offset. These
/// would be the offsets from paving from OUTSIDE_POLY and INSIDE_POLY polygons.
/// \param a_xyTol: tolerance for floating point comparisons
//------------------------------------------------------------------------------
void MeIntersectPolys::SetupInIn(const std::vector<MePolyOffsetterOutput>& a_offsets,
                                 double a_xyTol)
{
  m_p->SetupInIn(a_offsets, a_xyTol);
} // InsidePolys::SetupInIn
//------------------------------------------------------------------------------
/// \brief calls implementation
/// \param a_offsets: The output from MePolyCleaner::IntersectCleanInPolys
/// \param a_xyTol: tolerance for floating point comparisons
//------------------------------------------------------------------------------
void MeIntersectPolys::SetupInOut(const MePolyOffsetterOutput& a_offsets, double a_xyTol)
{
  m_p->SetupInOut(a_offsets, a_xyTol);
} // InsidePolys::SetupInOut
//------------------------------------------------------------------------------
/// \brief Calculates the envelope of all of the polygons
//------------------------------------------------------------------------------
void MeIntersectPolys::CalcEnvelopes()
{
  m_p->m_envelopes.resize(0);
  for (size_t i = 0; i < m_p->m_loops.size(); ++i)
  {
    m_p->m_envelopes.push_back(iCalcPolyEnvelope(m_p->m_loops[i], m_p->m_polyPts.Pts()));
  }
} // InsidePolys::CalcEnvelopes
//------------------------------------------------------------------------------
/// \brief calls implementation
//------------------------------------------------------------------------------
void MeIntersectPolys::InInTrivialPolyCases()
{
  m_p->InInTrivialPolyCases();
} // InsidePolys::InInTrivialPolyCases
//------------------------------------------------------------------------------
/// \brief calls implementation
//------------------------------------------------------------------------------
void MeIntersectPolys::InInDoIntersection()
{
  m_p->InInDoIntersection();
} // InsidePolys::InInDoIntersection
//------------------------------------------------------------------------------
/// \brief calls implementation
//------------------------------------------------------------------------------
void MeIntersectPolys::InOutDoIntersection()
{
  m_p->InOutDoIntersection();
} // InsidePolys::InOutDoIntersection
//------------------------------------------------------------------------------
/// \brief calls implementation
/// \param a_out: Holds the output polygons from the operation performed (either
/// InInDoIntersection or InOutDoIntersection
//------------------------------------------------------------------------------
void MeIntersectPolys::FillOutput(MePolyOffsetterOutput& a_out)
{
  m_p->FillOutput(a_out);
} // MpInsidePolys::FillOutput
//------------------------------------------------------------------------------
/// \brief calls implementation
/// \param a_input: Input polygons
/// \param a_output: 2D Vector where the first dimension is the number of
/// "outside" polygons.
//------------------------------------------------------------------------------
void MeIntersectPolys::ClassifyPolys(const MePolyOffsetterOutput& a_input,
                                     std::vector<std::vector<size_t>>& a_output)
{
  m_p->ClassifyPolys(a_input, a_output);
} // MpInsidePolys::ClassifyPolys

////////////////////////////////////////////////////////////////////////////////
/// \class MeIntersectPolys::impl
/// \brief Does polygon intersection for MePolyCleaner
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief Setup for the InInDoIntersection
/// \param a_offsets Vector of outputs from MePolyOffsetter::Offset. These
/// would be the offsets from paving from OUTSIDE_POLY and INSIDE_POLY polygons.
/// \param a_xyTol tolerance for floating point comparisons
//------------------------------------------------------------------------------
void MeIntersectPolys::impl::SetupInIn(const std::vector<MePolyOffsetterOutput>& a_offsets,
                                       double a_xyTol)
{
  m_polyPts.XyTol() = a_xyTol;
  // get all points of INSIDE_POLY polys and put into a single
  // array and update indexes for segment definitions.
  std::vector<Pt3d>& pts(m_polyPts.Pts());
  for (size_t i = 0; i < a_offsets.size(); ++i)
  {
    const MePolyOffsetterOutput& o(a_offsets[i]);
    // add the points and INSIDE_POLY loops
    size_t start = pts.size();
    pts.insert(pts.end(), o.m_pts.begin(), o.m_pts.end());
    for (size_t j = 0; j < o.m_loops.size(); ++j)
    {
      if (MePolyOffsetter::OUTSIDE_POLY == o.m_loopTypes[j] ||
          MePolyOffsetter::NEWOUT_POLY == o.m_loopTypes[j])
      {
        m_out.m_loopTypes.push_back(o.m_loopTypes[j]);
        m_out.m_loops.push_back(o.m_loops[j]);
        for (size_t k = 0; start > 0 && k < m_out.m_loops.back().size(); ++k)
        {
          m_out.m_loops.back()[k] += start;
        }
        continue;
      }

      m_loops.push_back(o.m_loops[j]);
      for (size_t k = 0; start > 0 && k < m_loops.back().size(); ++k)
      {
        m_loops.back()[k] += start;
      }
    }
  }
  OnSetup();
} // InsidePolys::impl::SetupInIn
//------------------------------------------------------------------------------
/// \brief Setup member variables. Called by SetupInIn and SetupInOut
//------------------------------------------------------------------------------
void MeIntersectPolys::impl::OnSetup()
{
  // hash the points and handle duplicates
  std::vector<size_t> idx(m_polyPts.HashPts());
  // update inside poly loops with indexes from hashing
  UpdateInPolyLoopsFromHash(idx);
  // create boost polygons
  for (size_t i = 0; i < m_loops.size(); ++i)
  {
    m_bPolys.push_back(BoostPoly(i));
  }
} // MeIntersectPolys::impl::OnSetup
//------------------------------------------------------------------------------
/// \brief Setup for InOutDoIntersection
/// \param a_offsets The output from MePolyCleaner::IntersectCleanInPolys
/// \param a_xyTol tolerance for floating point comparisons
//------------------------------------------------------------------------------
void MeIntersectPolys::impl::SetupInOut(const MePolyOffsetterOutput& a_offsets, double a_xyTol)
{
  m_polyPts.XyTol() = a_xyTol;
  // get all points of OUTSIDE_POLY & INSIDE_POLY polys and put into a single
  // array and update indexes for segment definitions.
  std::vector<Pt3d>& pts(m_polyPts.Pts());
  const MePolyOffsetterOutput& o(a_offsets);
  pts.insert(pts.end(), o.m_pts.begin(), o.m_pts.end());
  for (size_t j = 0; j < o.m_loops.size(); ++j)
  {
    if (MePolyOffsetter::NEWOUT_POLY == o.m_loopTypes[j])
    {
      m_out.m_loopTypes.push_back(o.m_loopTypes[j]);
      m_out.m_loops.push_back(o.m_loops[j]);
      continue;
    }
    m_loops.push_back(o.m_loops[j]);
    m_loopTypes.push_back(o.m_loopTypes[j]);
  }
  OnSetup();
} // MeIntersectPolys::impl::SetupInOut
//------------------------------------------------------------------------------
/// \brief Looks for polygons that are completely disjoint based on their
/// envelopes. Also finds polygons that are inside of other polygons. If an
/// INSIDE_POLY is inside of another INSIDE_POLY then it will be deleted.
/// However, it is possible for an INSIDE_POLY to be inside of a NEWOUT_POLY.
/// If this is the case then we don't delete that polygon.
//------------------------------------------------------------------------------
void MeIntersectPolys::impl::InInTrivialPolyCases()
{
  SetIdx delPolys;
  ClassifyDisjointPolys(delPolys);
  ClassifyPolysInsideOfPolys(delPolys);
  ClassifyDisjointPolys(delPolys);

  // sort the remaining polys based on area of envelope and build a stack
  // to process
  std::multimap<double, size_t> areaMap;
  std::vector<GmBstBox3d>& env(m_envelopes);
  for (size_t i = 0; i < env.size(); ++i)
  {
    if (delPolys.find(i) == delPolys.end())
    {
      GmBstBox3d& b(env[i]);
      double area = (b.max_corner().x - b.min_corner().x) * (b.max_corner().y - b.min_corner().y);
      areaMap.insert(std::make_pair(area, i));
    }
  }
  std::multimap<double, size_t>::reverse_iterator rit = areaMap.rbegin();
  for (; rit != areaMap.rend(); ++rit)
  {
    m_stack.push_back(rit->second);
  }
} // InsidePolys::impl::InInTrivialPolyCases
//------------------------------------------------------------------------------
/// \brief Uses a stack of INSIDE_POLY polygons to intersect the polygons with
/// one another. The stack is sorted by the size of the polygon envelope. The
/// largest envelopes are processed first. If 2 polygons do intersect then a
/// new polygon is formed that is the union of the 2 polygons. The 2 polygons
/// are removed from the stack and the new polygon is put on the front of the
/// stack and intersected with the remaining polygons. If a polygon does not
/// intersect any other polygon then it is removed from the stack and moved to
/// the output.
//------------------------------------------------------------------------------
void MeIntersectPolys::impl::InInDoIntersection()
{
  while (!m_stack.empty())
  {
    std::list<size_t>::iterator it = m_stack.begin();
    size_t idx = *it, idx2;

    bool processed = false;
    std::list<size_t>::iterator it2 = it;
    ++it2;
    while (!processed && it2 != m_stack.end())
    {
      idx2 = *it2;
      // find potentially intersecting poly
      if (iEnvelopesOverlap(idx, idx2, m_envelopes, m_envelopes))
      { // process the 2 polygons. If they intersect a new poly is created and
        // the original 2 are removed from the stack. A new OUTSIDE_POLY can
        // also be created in this process; it is moved to the output
        processed = BoostPolyUnion(idx, idx2);
      }
      if (processed)
      {
        m_stack.erase(it);
        m_stack.erase(it2);
      }
      else
        ++it2;
    }

    if (!processed)
    { // move this loop to the output
      m_stack.pop_front();
      m_bOutPolys.push_back(m_bPolys[idx]);
      m_bOutPolyType.push_back(MePolyOffsetter::INSIDE_POLY);
    }
  }
} // InsidePolys::impl::InInDoIntersection
//------------------------------------------------------------------------------
/// \brief A stack is created of INSIDE_POLY polygons. Each of these polygons
/// is checked to see if it intersects with OUTSIDE_POLY polygons. If they
/// intersect then the INSIDE_POLY is subtracted from the OUTSIDE_POLY and the
/// result from that operation can be multiple new OUTSIDE_POLY polygons. The
/// old OUTSIDE_POLY is removed from consideration for intersections and the
/// results from the polygon subtraction are now considered for intersection
/// with other INSIDE_POLY polygons. Any INSIDE_POLY polygon that intersects an
/// OUTSIDE_POLY is deleted. Otherwise the INSIDE_POLY moves to the output.
//------------------------------------------------------------------------------
void MeIntersectPolys::impl::InOutDoIntersection()
{
  InOutCalcStartingStack();
  while (!m_stack.empty())
  {
    std::list<size_t>::iterator it = m_stack.begin();
    size_t idx = *it;
    std::vector<size_t> polys(InOutFindIntersectingPolys(idx));

    bool processed(false);
    for (size_t i = 0; i < polys.size(); ++i)
    {
      if (BoostPolySubtract(polys[i], idx))
      {
        processed = true;
      }
    }
    m_stack.pop_front();
    if (!processed)
    { // move this loop to the output
      m_bOutPolys.push_back(m_bPolys[idx]);
      m_bOutPolyType.push_back(MePolyOffsetter::INSIDE_POLY);
    }
  }
  // move the left over OUTSIDE_POLY to the output
  SetIdx::iterator it(m_inOutOpolys.begin()), iend(m_inOutOpolys.end());
  for (; it != iend; ++it)
  {
    m_bOutPolys.push_back(m_bPolys[*it]);
    m_bOutPolyType.push_back(MePolyOffsetter::OUTSIDE_POLY);
  }
} // MeIntersectPolys::impl::InOutDoIntersection
//------------------------------------------------------------------------------
/// \brief Fills the output class
/// \param a_out Holds the output polygons from the operation performed (either
/// InInDoIntersection or InOutDoIntersection
//------------------------------------------------------------------------------
void MeIntersectPolys::impl::FillOutput(MePolyOffsetterOutput& a_out)
{
  a_out = m_out;
  // put the boost polys into the output
  std::vector<size_t> loop;
  for (size_t i = 0; i < m_bOutPolys.size(); ++i)
  {
    GmBstPoly3d& p(m_bOutPolys[i]);
    loop.resize(0);
    if (m_bOutPolyType[i] == MePolyOffsetter::INSIDE_POLY)
    {
      for (size_t j = 1; j < p.outer().size(); ++j)
      {
        Pt3d& pt(p.outer()[j]);
        size_t idx = m_polyPts.IdxFromPt3d(pt);
        loop.push_back(idx);
      }
      std::reverse(loop.begin(), loop.end());
    }
    else
    {
      for (size_t j = 0; j < p.outer().size() - 1; ++j)
      {
        Pt3d& pt(p.outer()[j]);
        size_t idx = m_polyPts.IdxFromPt3d(pt);
        loop.push_back(idx);
      }
    }
    a_out.m_loops.push_back(loop);
    a_out.m_loopTypes.push_back(m_bOutPolyType[i]);
  }
  a_out.m_pts = m_polyPts.Pts();
} // MpInsidePolys::impl::FillOutput
//------------------------------------------------------------------------------
/// \brief Finds INSIDE_POLY that are inside of NEWOUT_POLY or OUTSIDE_POLY.
/// Used by MeIntersectPolys::impl::ClassifyPolys.
/// \param a_oPoly A set of indices to "outside" polygons
/// (either NEWOUT_POLY or OUTSIDE_POLY)
/// \param a_psetUsed A sect of indices to polygons that have already been
/// classified
/// \param a_output 2D Vector where the first dimension is the number of
/// "outside" polygons.
//------------------------------------------------------------------------------
void MeIntersectPolys::impl::FillPolysInsideOfPolys(std::set<size_t>& a_oPoly,
                                                    std::set<size_t>& a_psetUsed,
                                                    std::vector<std::vector<size_t>>& a_output)
{
  std::vector<Pt3d>& pts(m_polyPts.Pts());
  std::set<size_t>::iterator it(a_oPoly.begin()), itend(a_oPoly.end()), usedEnd(a_psetUsed.end());
  for (; it != itend; ++it)
  {
    GmBstPoly3d poly = BoostPoly(*it);
    a_output.push_back(std::vector<size_t>());
    a_output.back().push_back(*it);
    // for this method, no polygons should intersect so we just need to see
    // if the first point of a polygon is inside of this polygon
    for (size_t i = 0; i < m_loops.size(); ++i)
    {
      if (a_psetUsed.find(i) != usedEnd)
        continue;

      Pt3d& p(pts[m_loops[i].front()]);
      if (bg::within(p, poly))
      {
        a_output.back().push_back(i);
        a_psetUsed.insert(i);
      }
    }
  }
} // MeIntersectPolys::impl::FillPolysInsideOfPolys
//------------------------------------------------------------------------------
/// \brief Classifies the polygons. see PolyClassifierImpl::Classify
/// \param a_input Input polygons
/// \param a_output 2D Vector where the first dimension is the number of
/// "outside" polygons.
//------------------------------------------------------------------------------
void MeIntersectPolys::impl::ClassifyPolys(const MePolyOffsetterOutput& a_input,
                                           std::vector<std::vector<size_t>>& a_output)
{
  a_output.resize(0);
  m_polyPts.Pts() = a_input.m_pts;
  m_loops = a_input.m_loops;
  m_loopTypes = a_input.m_loopTypes;

  std::set<size_t> psetNewOut, psetOut, psetUsed;
  // first find all NEWOUT_POLY polygons
  for (size_t i = 0; i < m_loopTypes.size(); ++i)
  {
    if (MePolyOffsetter::NEWOUT_POLY == m_loopTypes[i])
    {
      psetNewOut.insert(i);
      psetUsed.insert(i);
    }
    else if (MePolyOffsetter::OUTSIDE_POLY == m_loopTypes[i])
    {
      psetOut.insert(i);
      psetUsed.insert(i);
    }
  }
  // find INSIDE_POLY that are inside of NEWOUT_POLY
  FillPolysInsideOfPolys(psetNewOut, psetUsed, a_output);
  // find INSIDE_POLY that are inside of OUTSIDE_POLY
  FillPolysInsideOfPolys(psetOut, psetUsed, a_output);
} // MeIntersectPolys::impl::ClassifyPolys
//------------------------------------------------------------------------------
/// \brief Classifies polygons that are completely disjoint from other polygons
/// based on the envelopes of the polygons.
/// \param a_delPolys A set of indexes to polygons that have been found to not
/// intersect any other polygons.
//------------------------------------------------------------------------------
void MeIntersectPolys::impl::ClassifyDisjointPolys(SetIdx& a_delPolys)
{
  std::vector<GmBstBox3d>& env(m_envelopes);
  // find polys that are completely out of other poly envelopes
  std::vector<int> overlaps(env.size(), 0);
  for (size_t i = 0; i < env.size(); ++i)
  {
    if (a_delPolys.find(i) != a_delPolys.end())
      continue;

    for (size_t j = 0; j < env.size(); j++)
    {
      if (i == j)
        continue;
      if (overlaps[j])
        continue;
      if (a_delPolys.find(j) != a_delPolys.end())
        continue;

      if (iEnvelopesOverlap(i, j, env, env))
      {
        overlaps[i] = overlaps[j] = 1;
      }
    }
  }
  for (size_t i = 0; i < overlaps.size(); ++i)
  { // anything that doesn't overlap move to the output and mark the poly
    // as deleted
    if (!overlaps[i] && a_delPolys.find(i) == a_delPolys.end())
    {
      m_out.m_loops.push_back(m_loops[i]);
      m_out.m_loopTypes.push_back(MePolyOffsetter::INSIDE_POLY);
      a_delPolys.insert(i);
    }
  }
} // InsidePolys::impl::ClassifyDisjointPolys
//------------------------------------------------------------------------------
/// \brief If an INSIDE_POLY is inside of another INSIDE_POLY then it will be
/// deleted. But if it is inside of a NEWOUT_POLY then don't delete it.
/// \param a_delPolys A set of indexes to polygons that have been found to not
/// intersect any other polygons.
//------------------------------------------------------------------------------
void MeIntersectPolys::impl::ClassifyPolysInsideOfPolys(SetIdx& a_delPolys)
{
  std::vector<GmBstBox3d>& env(m_envelopes);
  for (size_t i = 0; i < env.size(); ++i)
  {
    if (a_delPolys.find(i) != a_delPolys.end())
      continue;

    for (size_t j = 0; j < env.size(); j++)
    {
      if (i == j)
        continue;
      if (a_delPolys.find(j) != a_delPolys.end())
        continue;

      if (iEnvelopeInsideOfEnvelope(i, j, env))
      {
        if (m_polyPts.PolyInsideOfPoly(m_loops[i], m_loops[j]))
        {
          m_polyInsideOfPoly.insert(std::make_pair(i, j));
          // make sure the poly is not inside of an NEWOUT_POLY
          bool inNewOut(false);
          for (size_t k = 0; !inNewOut && k < m_out.m_loops.size(); ++k)
          {
            if (m_out.m_loopTypes[k] == MePolyOffsetter::NEWOUT_POLY)
            {
              std::vector<size_t> idx(m_out.m_loops[k]);
              std::reverse(idx.begin(), idx.end());
              if (m_polyPts.PolyInsideOfPoly(idx, m_loops[j]))
              {
                inNewOut = true;
              }
            }
          }
          if (!inNewOut)
            a_delPolys.insert(j);
        }
      }
    }
  }
} // MeIntersectPolys::impl::ClassifyPolysInsideOfPolys
//------------------------------------------------------------------------------
/// \brief Updates the indexes that define the polygons after the point
/// locations have been hashed to find duplicate locations.
/// \param a_idx Vector with new indexes based on the hashing of point
/// locations. If there are no duplicates the values would be 0,1,2,3,...
/// If there are duplicates then the values might be 0,1,2,0,2,5,6,0...
//------------------------------------------------------------------------------
void MeIntersectPolys::impl::UpdateInPolyLoopsFromHash(const std::vector<size_t>& a_idx)
{
  size_t idx;
  std::vector<std::vector<size_t>>& loops(m_loops);
  for (size_t i = 0; i < loops.size(); ++i)
  {
    for (size_t j = 0; j < loops[i].size(); ++j)
    {
      idx = loops[i][j];
      loops[i][j] = a_idx[idx];
    }
  }
} // InsidePolys::impl::UpdateInPolyLoopsFromHash
//------------------------------------------------------------------------------
/// \brief Sets up the stack for the InOutDoIntersection method. A stack of
/// INSIDE_POLY polygons is created and sorted by polygon envelope from largest
/// to smallest. A set of OUTSIDE_POLY polygons is created. These are candidates
/// for intersection with the INSIDE_POLY polygons.
//------------------------------------------------------------------------------
void MeIntersectPolys::impl::InOutCalcStartingStack()
{
  std::vector<size_t> inPolys;
  for (size_t i = 0; i < m_loopTypes.size(); ++i)
  {
    if (MePolyOffsetter::OUTSIDE_POLY == m_loopTypes[i])
      m_inOutOpolys.insert(i);
    else if (MePolyOffsetter::INSIDE_POLY == m_loopTypes[i])
      inPolys.push_back(i);
  }
  // sort the polys by envelope size
  std::multimap<double, size_t> areaMap;
  for (size_t i = 0; i < inPolys.size(); ++i)
  {
    GmBstBox3d& b(m_envelopes[inPolys[i]]);
    double area = (b.max_corner().x - b.min_corner().x) * (b.max_corner().y - b.min_corner().y);
    areaMap.insert(std::make_pair(area, inPolys[i]));
  }
  std::multimap<double, size_t>::reverse_iterator rit = areaMap.rbegin();
  for (; rit != areaMap.rend(); ++rit)
  {
    m_stack.push_back(rit->second);
  }
} // MeIntersectPolys::impl::InOutCalcStartingStack
//------------------------------------------------------------------------------
/// \brief Finds OUTSIDE_POLY polygons that potentially intersect the
/// INSIDE_POLY "a_poly". Used by InOutDoIntersection.
/// \param a_poly Index to an INSIDE_POLY polygon.
/// \return Vector of polygon indexes for the polys that intersect.
//------------------------------------------------------------------------------
std::vector<size_t> MeIntersectPolys::impl::InOutFindIntersectingPolys(size_t a_poly)
{
  std::vector<size_t> outPolys;
  // find potentially intersecting OUTSIDE_POLY
  std::vector<GmBstBox3d>& env(m_envelopes);
  SetIdx::iterator it(m_inOutOpolys.begin()), iend(m_inOutOpolys.end());
  for (; it != iend; ++it)
  {
    if (iEnvelopesOverlap(a_poly, *it, env, env))
    {
      if (!m_polyPts.PolyInsideOfPoly(m_bPolys[*it].outer(), m_loops[a_poly]))
      {
        outPolys.push_back(*it);
      }
    }
  }
  return outPolys;
} // MeIntersectPolys::impl::InOutFindIntersectingPolys
//------------------------------------------------------------------------------
/// \brief Creates a boost polygon from a vector of point indexes
/// \param a_loopIdx Index to a polygon. The point indexes that make up the
/// polygon will be in the m_loops[m_loopIdx] variable. These indexes refer to
/// locations in the MePolyPts::Pts() vector.
/// \return A boost polygon class.
//------------------------------------------------------------------------------
GmBstPoly3d MeIntersectPolys::impl::BoostPoly(size_t a_loopIdx)
{
  GmBstPoly3d poly;
  std::vector<Pt3d>& pts(m_polyPts.Pts());
  std::vector<size_t>& l(m_loops[a_loopIdx]);
  if (m_loopTypes.empty() || m_loopTypes[a_loopIdx] == MePolyOffsetter::INSIDE_POLY)
  {
    std::vector<size_t>::reverse_iterator it(l.rbegin()), iend(l.rend());
    bg::exterior_ring(poly).push_back(pts[l[0]]);
    for (; it != iend; ++it)
    {
      bg::exterior_ring(poly).push_back(pts[*it]);
    }
  }
  else
  {
    std::vector<size_t>::iterator it(l.begin()), iend(l.end());
    for (; it != iend; ++it)
    {
      bg::exterior_ring(poly).push_back(pts[*it]);
    }
    bg::exterior_ring(poly).push_back(pts[l[0]]);
  }
  return poly;
} // MeInsidePolys::impl::BoostPoly
//------------------------------------------------------------------------------
/// \brief Performs a union operation on 2 polygons. Used by InInDoIntersection
/// \param a_i Index to a polygon (m_loops[a_i]).
/// \param a_j Index to a polygon (m_loops[a_j]).
/// \return true if the union resulted in a new polygon.
//------------------------------------------------------------------------------
bool MeIntersectPolys::impl::BoostPolyUnion(size_t a_i, size_t a_j)
{
  std::vector<GmBstPoly3d> out;
  GmBstPoly3d &iPoly(m_bPolys[a_i]), &jPoly(m_bPolys[a_j]);
  std::pair<size_t, size_t> p1(a_i, a_j), p2(a_j, a_i);
  if (m_polyInsideOfPoly.find(p1) != m_polyInsideOfPoly.end() ||
      m_polyInsideOfPoly.find(p2) != m_polyInsideOfPoly.end())
    return false;
  bg::union_(iPoly, jPoly, out);
  if (out.size() != 1)
    return false;
  m_stack.push_back(m_bPolys.size());
  m_bPolys.push_back(out[0]);
  m_bPolys.back().inners().clear();
  std::vector<Pt3d>& pts(m_bPolys.back().outer());
  std::vector<size_t> loop(pts.size());
  for (size_t i = 0; i < loop.size(); ++i)
    loop[i] = i;
  m_envelopes.push_back(iCalcPolyEnvelope(loop, pts));
  // put any inner polygons into the output
  for (size_t i = 0; i < out[0].inners().size(); ++i)
  {
    GmBstPoly3d innerPoly;
    std::vector<Pt3d>& vi(out[0].inners()[i]);
    std::vector<Pt3d>::reverse_iterator it(vi.rbegin()), iend(vi.rend());
    for (; it != iend; ++it)
    {
      bg::exterior_ring(innerPoly).push_back(*it);
    }
    m_bOutPolys.push_back(innerPoly);
    m_bOutPolyType.push_back(MePolyOffsetter::NEWOUT_POLY);
  }
  return true;
} // MeInsidePolys::impl::BoostPolyUnion
//------------------------------------------------------------------------------
/// \brief Performs a subtraction with 2 polygons. Used by InOutDoIntersection
/// \param a_i Index to a polygon (m_loops[a_i]).
/// \param a_j Index to a polygon (m_loops[a_j]). a_j polygon is subtracted
/// from a_i polygon
/// \return true if the intersection resulted in a new polygon
//------------------------------------------------------------------------------
bool MeIntersectPolys::impl::BoostPolySubtract(size_t a_i, size_t a_j)
{
  std::vector<GmBstPoly3d> out;
  GmBstPoly3d &iPoly(m_bPolys[a_i]), &jPoly(m_bPolys[a_j]);
  if (!bg::intersects(iPoly, jPoly))
    return false;
  bg::difference(iPoly, jPoly, out);
  m_inOutOpolys.erase(a_i);
  for (size_t i = 0; i < out.size(); ++i)
  {
    m_inOutOpolys.insert(m_bPolys.size());
    m_bPolys.push_back(out[i]);
    m_bPolys.back().inners().clear();
    std::vector<Pt3d>& pts(m_bPolys.back().outer());
    std::vector<size_t> loop(pts.size());
    for (size_t i = 0; i < loop.size(); ++i)
      loop[i] = i;
    m_envelopes.push_back(iCalcPolyEnvelope(loop, pts));
    m_loopTypes.push_back(MePolyOffsetter::OUTSIDE_POLY);
  }
  return true;
} // MeInsidePolys::impl::BoostPolyUnion

} // namespace xms

#ifdef CXX_TEST
////////////////////////////////////////////////////////////////////////////////

#include <xmsmesh/meshing/detail/MeIntersectPolys.t.h>

#include <xmscore/testing/TestTools.h>

// namespace xms {
using namespace xms;

////////////////////////////////////////////////////////////////////////////////
/// \class MeIntersectPolysUnitTests
/// \brief tester for the MeIntersectPolys class
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief tests classifying polygons
//------------------------------------------------------------------------------
void MeIntersectPolysUnitTests::testClassify()
{
  // x =     0                4                 8
  //
  // y=4     5----------------------------------6
  //        |                                  |
  //        |                    15------14    |
  //        |   1-----------2     |      |     |
  //        |   |           |    12------13    |
  //        |   |  19---18  |                  |
  // y=2     |   |   |   |   |     8------11    |
  //        |   |  16---17  |     |      |     |
  //        |   |           |     |      |     |
  //        |   0-----------3     9------10    |
  //        |                                  |
  // y=0     4----------------------------------7
  //
  // Poly 0,1,2,3 is a NEWOUT_POLY, 4,5,6,7 is OUTSIDE_POLY
  // all the others are INSIDE_POLY
  //
  MePolyOffsetterOutput o, o1;
  o.m_pts = {{1, 1, 0},   {1, 3, 0},     {4, 3, 0},     {4, 1, 0},     {0, 0, 0},
             {0, 4, 0},   {8, 4, 0},     {8, 0, 0},     {5, 1, 0},     {7, 1, 0},
             {7, 2, 0},   {5, 2, 0},     {5, 3, 0},     {7, 3, 0},     {7, 3.5, 0},
             {5, 3.5, 0}, {1.5, 1.5, 0}, {3.5, 1.5, 0}, {3.5, 2.5, 0}, {1.5, 2.5, 0}};
  std::vector<size_t> v0 = {0, 1, 2, 3};
  o.m_loops.push_back(v0);
  std::vector<size_t> v1 = {4, 5, 6, 7};
  o.m_loops.push_back(v1);
  std::vector<size_t> v2 = {8, 9, 10, 11};
  o.m_loops.push_back(v2);
  std::vector<size_t> v3 = {12, 13, 14, 15};
  o.m_loops.push_back(v3);
  std::vector<size_t> v4 = {16, 17, 18, 19};
  o.m_loops.push_back(v4);
  o.m_loopTypes.assign(5, MePolyOffsetter::INSIDE_POLY);
  o.m_loopTypes[0] = MePolyOffsetter::NEWOUT_POLY;
  o.m_loopTypes[1] = MePolyOffsetter::OUTSIDE_POLY;
  std::vector<std::vector<size_t>> output;
  MeIntersectPolys::impl pc;
  pc.ClassifyPolys(o, output);
  TS_ASSERT_EQUALS(2, output.size());
  if (2 != output.size())
    return;
  std::vector<size_t> baseIdx = {0, 4};
  TS_ASSERT_EQUALS_VEC(baseIdx, output[0]);
  baseIdx = {1, 2, 3};
  TS_ASSERT_EQUALS_VEC(baseIdx, output[1]);
} // MeIntersectPolysUnitTests::testClassify

//} // namespace xms
#endif
