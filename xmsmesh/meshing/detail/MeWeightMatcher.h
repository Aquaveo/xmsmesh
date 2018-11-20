#pragma once
//------------------------------------------------------------------------------
/// \file
/// \ingroup meshing_detail
/// \copyright (C) Copyright Aquaveo 2018. Distributed under FreeBSD License
/// (See accompanying file LICENSE or https://aqaveo.com/bsd/license.txt)
//------------------------------------------------------------------------------

//----- Included files ---------------------------------------------------------

// 3. Standard library headers
#include <vector>

// 4. External library headers

// 5. Shared code headers
#include <xmscore/misc/base_macros.h>
#include <xmscore/misc/boost_defines.h>
#include <xmscore/stl/vector.h>

//----- Forward declarations ---------------------------------------------------

//----- Namespace declaration --------------------------------------------------

namespace xms
{
//----- Constants / Enumerations -----------------------------------------------

//----- Structs / Classes ------------------------------------------------------

/// Representation of an edge or pseudo edge between two faces with a weight.
/// The two faces may be adjacent by sharing a common interior edge of the
/// mesh or they may be adjacent to a common boundary vertex of the mesh but
/// separated by at least one other face in between them that also shares that
/// boundary edge.
struct MeEdge
{
  /// Constructor.
  /// \param[in] a_f0 Index of one of the adjacent faces.
  /// \param[in] a_f1 Index of one of the adjacent faces.
  /// \param[in] a_weight The weight assigned to the edge.
  MeEdge(int a_f0, int a_f1, int a_weight)
  : m_f0(a_f0)
  , m_f1(a_f1)
  , m_weight(a_weight)
  {
  }

  /// Equals operator for test assertions.
  /// \param[in] a_rhs The ege to compare against.
  /// \return true if the edges are equal.
  bool operator==(const MeEdge& a_rhs) const
  {
    return m_f0 == a_rhs.m_f0 && m_f1 == a_rhs.m_f1 && m_weight == a_rhs.m_weight;
  }

  int m_f0;     ///< Index of one of the adjacent faces.
  int m_f1;     ///< Index of one of the adjacent faces.
  int m_weight; ///< The weight assigned to the edge.

}; // struct MeEdge

typedef std::vector<MeEdge> VecMeEdge; ///< Vector of MeEdge

class MeWeightMatcher
{
public:
  static BSHP<MeWeightMatcher> New();
  MeWeightMatcher();
  virtual ~MeWeightMatcher();

  /// \cond
  virtual VecInt MatchWeights(const VecMeEdge& a_edges, bool maxcardinality = false) = 0;

private:
  XM_DISALLOW_COPY_AND_ASSIGN(MeWeightMatcher);
  /// \endcond
}; // class MeWeightMatcher

//----- Function prototypes ----------------------------------------------------

} // namespace xms
