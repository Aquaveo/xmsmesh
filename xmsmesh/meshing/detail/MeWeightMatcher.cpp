//------------------------------------------------------------------------------
/// \file
/// \ingroup meshing_detail
/// \copyright (C) Copyright Aquaveo 2018. Distributed under FreeBSD License
/// (See accompanying file LICENSE or https://aqaveo.com/bsd/license.txt)
//------------------------------------------------------------------------------

// defines for debugging
//#define DEBUG
//#define _DEBUG 1
//#define CHECK_DELTA
//#define CHECK_OPTIMUM
#ifdef DEBUG
#include <xmsmesh/meshing/detail/MeWeightMatcher.t.h>
#endif

//----- Included files ---------------------------------------------------------

// 1. Precompiled header

// 2. My own header
#include <xmsmesh/meshing/detail/MeWeightMatcher.h>

// 3. Standard library headers
#include <numeric>

// 4. External library headers

// 5. Shared code headers
#include <xmscore/misc/XmError.h>
#include <xmscore/stl/vector.h>

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

namespace
{
////////////////////////////////////////////////////////////////////////////////
/// \class VecIntPy
/// \brief Vector of ints to allow negative indexing and mimic other python
/// functions.
////////////////////////////////////////////////////////////////////////////////
class VecIntPy
{
public:
  /// Constructor.
  /// \param[in] a_size The desired size.
  /// \param[in] a_value The default value of elements.
  explicit VecIntPy(int a_size = 0, int a_value = 0)
  : m_vec(a_size, a_value)
  {
  }

  /// Copy constructor from VecInt.
  /// \param[in] a_vec The vector to copy values from.
  VecIntPy(const VecInt& a_vec)
  : m_vec(a_vec)
  {
  }

  /// Copy constructor.
  /// \param[in] a_vec The vector to copy values from.
  VecIntPy(const VecIntPy& a_vec)
  : m_vec(a_vec.m_vec)
  {
  }

  /// Assignment operator.
  /// \param[in] a_rhs The vector to assign values from.
  const VecIntPy& operator=(const VecIntPy& a_rhs)
  {
    m_vec = a_rhs.m_vec;
    return *this;
  }

  /// Indexing operator returning a l-value to the element.
  /// \param[in] a_idx The index to extract the value from. Negative index
  /// retrieves values backwards from end with -1 being the last element.
  /// \return The indexed element.
  int& operator[](int a_idx)
  {
    if (a_idx < 0)
    {
      a_idx = (int)m_vec.size() + a_idx;
    }
    XM_ASSERT(a_idx >= 0 && a_idx < m_vec.size());
    return m_vec[a_idx];
  }

  /// Indexing operator returning a copy of the element.
  /// \param[in] a_idx The index to extract the value from. Negative index
  /// retrieves values backwards from end with -1 being the last element.
  /// \return The indexed element.
  int operator[](int a_idx) const
  {
    if (a_idx < 0)
    {
      a_idx = (int)m_vec.size() + a_idx;
    }
    XM_ASSERT(a_idx >= 0 && a_idx < m_vec.size());
    return m_vec[a_idx];
  }

  /// Convert array to a VecInt.
  /// \return A constant reference to a VecInt of values.
  const VecInt& ToVecInt() const { return m_vec; }

  /// Get the size of the vector.
  /// \return The size of the vector.
  int size() const { return (int)m_vec.size(); }

  /// Is the vector empty?
  /// \return true if the vector is empty.
  bool empty() const { return m_vec.empty(); }

  /// Find the index of a value in the vector. Value must be in the vector.
  /// \param[in] a_value The value to find in the vector.
  /// \return The index of the value.
  int index(int a_value)
  {
    auto it = std::find(m_vec.begin(), m_vec.end(), a_value);
    XM_ASSERT(it != m_vec.end());
    return (int)(it - m_vec.begin());
  }

  /// Get the minimum value of the vector. Assumes the vector is not empty.
  /// \return The minimum value of all elements.
  int min() const
  {
    auto it = std::min_element(m_vec.begin(), m_vec.end());
    XM_ASSERT(it != m_vec.end());
    return *it;
  }

  /// Get the minimum value of the vector between two given indices. Assumes the
  /// vector contains those indices.
  /// \param[in] a_start The first index of the range of elements.
  /// \param[in] a_end One past the last index of the range of elements.
  /// \return The minimum value of all elements in the range.
  int min(int a_start, int a_end) const
  {
    auto it = std::min_element(m_vec.begin() + a_start, m_vec.begin() + a_end);
    XM_ASSERT(a_start < a_end);
    XM_ASSERT(a_start >= 0 && a_end <= (int)m_vec.size());
    XM_ASSERT(it != m_vec.begin() + a_end);
    return *it;
  }

  /// Get the maximum value of the vector. Assumes the vector is not empty.
  /// \return The maximum value of all elements.
  int max() const
  {
    auto it = std::max_element(m_vec.begin(), m_vec.end());
    XM_ASSERT(it != m_vec.end());
    return *it;
  }

  /// Get the maximum value of the vector between two given indices. Assumes the
  /// vector contains those indices.
  /// \param[in] a_start The first index of the range of elements.
  /// \param[in] a_end One past the last index of the range of elements.
  /// \return The maximum value of all elements in the range.
  int max(int a_start, int a_end) const
  {
    auto it = std::max_element(m_vec.begin() + a_start, m_vec.begin() + a_end);
    XM_ASSERT(a_start < a_end);
    XM_ASSERT(a_start >= 0 && a_end <= (int)m_vec.size());
    XM_ASSERT(it != m_vec.begin() + a_end);
    return *it;
  }

  /// Append a value to the vector.
  /// \param[in] a_value The value to append.
  void push_back(int a_value) { m_vec.push_back(a_value); }

  /// Get the last item of the vector. Assumes the vector is not empty.
  /// \return The last item of the vector.
  int back()
  {
    XM_ASSERT(!m_vec.empty());
    return m_vec.back();
  }

  /// Remove the last item from the vector and get its value. Assume vector is
  /// not empty.
  /// \return The last item of the vector.
  int pop()
  {
    XM_ASSERT(!m_vec.empty());
    int back = m_vec.back();
    m_vec.pop_back();
    return back;
  }

  /// Appends the values of a vector to the end of the vector.
  /// \param[in] a_vec The vector containing values to append.
  void append(const VecIntPy& a_vec)
  {
    m_vec.insert(m_vec.end(), a_vec.m_vec.begin(), a_vec.m_vec.end());
  }

  /// Remove all elements from the vector.
  void clear() { m_vec.clear(); }

  /// Resize the vector to a given size and default value. Only new elements
  /// beyond the current size get the new value. If the size is smaller elements
  /// are removed.
  /// \param[in] a_size The new size of the vector.
  /// \param[in] a_default The default value for any new elements.
  void resize(int a_size, int a_default) { m_vec.resize(a_size, a_default); }

  /// Reserve space in the vector to minimize future allocations.
  /// \param[in] a_size The number of elements to reserve space for.
  void reserve(int a_size) { m_vec.reserve(a_size); }

  /// Resize a vector and fill with a elements with counted values
  /// \param[in] a_size The new size of the vector.
  /// \param[in] a_start_at The value to start counting from at the first
  /// element.
  void resize_to_count(int a_size, int a_start_at = 0)
  {
    m_vec.resize(a_size);
    std::iota(m_vec.begin(), m_vec.end(), a_start_at);
  }

  /// Reverse the elements of the vector.
  void reverse() { std::reverse(m_vec.begin(), m_vec.end()); }

  /// Rotate the elements of the vector so a given element is now the first.
  /// \param[in] a_frontIdx The element to be rotated to the first.
  void rotate(int a_frontIdx)
  {
    std::rotate(m_vec.begin(), m_vec.begin() + a_frontIdx, m_vec.end());
  }

private:
  VecInt m_vec; ///< The vector storage.
};

////////////////////////////////////////////////////////////////////////////////
/// \class VecInt2dPy
/// \brief Vector of vector of ints to allow negative indexing and mimic other
/// python functions.
////////////////////////////////////////////////////////////////////////////////
class VecInt2dPy
{
public:
  /// Constructor.
  /// \param[in] a_size The desired size.
  /// \param[in] a_value The default value of elements.
  explicit VecInt2dPy(int a_size = 0, VecIntPy a_value = VecIntPy())
  : m_vec(a_size, a_value)
  {
  }

  /// Copy constructor.
  /// \param[in] a_vec The vector to copy values from.
  VecInt2dPy(const VecInt2dPy& a_vec)
  : m_vec(a_vec.m_vec)
  {
  }

  /// Copy constructor from VecInt2d.
  /// \param[in] a_vec The vector to copy values from.
  VecInt2dPy(const VecInt2d& a_vec)
  : m_vec()
  {
    for (auto& vec : a_vec)
    {
      m_vec.push_back(VecIntPy(vec));
    }
  }

  /// Indexing operator returning a l-value to the element.
  /// \param[in] a_idx The index to extract the value from. Negative index
  /// retrieves values backwards from end with -1 being the last element.
  /// \return The indexed element.
  VecIntPy& operator[](int a_idx)
  {
    if (a_idx < 0)
    {
      a_idx = (int)m_vec.size() + a_idx;
    }
    XM_ASSERT(a_idx >= 0 && a_idx < m_vec.size());
    return m_vec[a_idx];
  }

  /// Indexing operator returning a copy of the element.
  /// \param[in] a_idx The index to extract the value from. Negative index
  /// retrieves values backwards from end with -1 being the last element.
  /// \return The indexed element.
  const VecIntPy& operator[](int a_idx) const
  {
    if (a_idx < 0)
    {
      a_idx = (int)m_vec.size() + a_idx;
    }
    XM_ASSERT(a_idx >= 0 && a_idx < m_vec.size());
    return m_vec[a_idx];
  }

  /// Append a value to the vector.
  /// \param[in] a_value The value to append.
  void push_back(const VecIntPy& a_value) { m_vec.push_back(a_value); }

  /// Get the size of the vector.
  /// \return The size of the vector.
  int size() const { return (int)m_vec.size(); }

  /// Is the vector empty?
  /// \return true if the vector is empty.
  bool empty() const { return m_vec.empty(); }

  /// Resize the vector to a given size and default value. Only new elements
  /// beyond the current size get the new value. If the size is smaller elements
  /// are removed.
  /// \param[in] a_size The new size of the vector.
  /// \param[in] a_default The default value for any new elements.
  void resize(int a_size, const VecIntPy& a_vec = VecIntPy()) { m_vec.resize(a_size, a_vec); }

private:
  std::vector<VecIntPy> m_vec; ///< The vector storage.
};

class MeWeightMatcherImpl : public MeWeightMatcher
{
public:
  MeWeightMatcherImpl();

  void Init(const VecMeEdge& a_edges);
  virtual VecInt MatchWeights(const VecMeEdge& a_edges, bool a_maxCardinality = false) override;

  int slack(int a_edge);
  VecIntPy blossomLeaves(int b);
  void assignLabel(int w, int t, int p);
  int scanBlossom(int v, int w);
  void addBlossom(int base, int k);
  void expandBlossom(int b, bool endstage);
  void augmentBlossom(int b, int v);
  void augmentMatching(int k);
  void verifyOptimum(bool maxcardinality);
  void checkDelta2();
  void checkDelta3();

private:
  const MeEdge& GetEdge(int a_idx);

  VecMeEdge m_edges; ///< Edges are numbered 0 .. (nedge-1).

  int m_nEdge;     ///< Number of edges in m_edges.
  int m_nVertex;   ///< Number of points in the mesh.
  int m_maxWeight; ///< Maximum of weights in m_edges.

  /// Edge endpoints are numbered 0 .. (2*m_nEdge-1), such that endpoints
  /// (2*k) and (2*k+1) both belong to edge k.
  /// If p is an edge endpoint,
  /// m_endPoint[p] is the vertex to which endpoint p is attached.
  /// Not modified by the algorithm.
  VecIntPy m_endPoint;

  /// If v is a vertex,
  /// m_neighbEnd[v] is the list of remote endpoints of the edges attached to v.
  /// Not modified by the algorithm.
  VecInt2dPy m_neighbEnd;

  /// If v is a vertex,
  /// m_mate[v] is the remote endpoint of its matched edge, or -1 if it is
  /// single (i.e. m_endPoint[m_mate[v]] is v's partner vertex).
  /// Initially all vertices are single; updated during augmentation.
  VecIntPy m_mate;

  /// If b is a top-level blossom,
  /// m_label[b] is 0 if b is unlabeled (free);
  ///             1 if b is an S-vertex/blossom;
  ///             2 if b is a T-vertex/blossom.
  /// The label of a vertex is found by looking at the label of its
  /// top-level containing blossom.
  /// If v is a vertex inside a T-blossom,
  /// m_label[v] is 2 iff v is reachable from an S-vertex outside the blossom.
  /// Labels are assigned during a stage and reset after each augmentation.
  VecIntPy m_label;

  /// If b is a labeled top-level blossom,
  /// m_labelEnd[b] is the remote endpoint of the edge through which b obtained
  /// its label, or -1 if b's base vertex is single.
  /// If v is a vertex inside a T-blossom and m_label[v] == 2,
  /// m_labelEnd[v] is the remote endpoint of the edge through which v is
  /// reachable from outside the blossom.
  VecIntPy m_labelEnd;

  /// If v is a vertex,
  /// m_inBlossom[v] is the top-level blossom to which v belongs.
  /// If v is a top-level vertex, v is itself a blossom (a trivial blossom)
  /// and m_inBlossom[v] == v.
  /// Initially all vertices are top-level trivial blossoms.
  VecIntPy m_inBlossom;

  /// If b is a sub-blossom,
  /// m_blossomParent[b] is its immediate parent (sub-)blossom.
  /// If b is a top-level blossom, m_blossomParent[b] is -1.
  VecIntPy m_blossomParent;

  /// If b is a non-trivial (sub-)blossom,
  /// m_blossomChildren[b] is an ordered list of its sub-blossoms, starting with
  /// the base and going round the blossom.
  /// m_blossomChildren = (2 * m_nVertex) * [ None ]
  VecInt2dPy m_blossomChildren;

  /// If b is a (sub-)blossom,
  /// m_blossomBase[b] is its base VERTEX (i.e. recursive sub-blossom).
  /// if m_nVertex is 4: [0, 1, 2, 3, -1, -1, -1, -1]
  VecIntPy m_blossomBase;

  /// If b is a non-trivial (sub-)blossom,
  /// m_blossomEndPts[b] is a list of endpoints on its connecting edges,
  /// such that m_blossomEndPts[b][i] is the local endpoint of
  /// m_blossomChildren[b][i] on the edge that connects it to
  /// m_blossomChildren[b][wrap(i+1)].
  VecInt2dPy m_blossomEndPts;

  /// If v is a free vertex (or an unreached vertex inside a T-blossom),
  /// m_bestEdge[v] is the edge to an S-vertex with least slack,
  /// or -1 if there is no such edge.
  /// If b is a (possibly trivial) top-level S-blossom,
  /// m_bestEdge[b] is the least-slack edge to a different S-blossom,
  /// or -1 if there is no such edge.
  /// This is used for efficient computation of delta2 and delta3.
  VecIntPy m_bestEdge;

  /// If b is a non-trivial top-level S-blossom,
  /// m_blossomBestEdges[b] is a list of least-slack edges to neighbouring
  /// S-blossoms, or None if no such list has been computed yet.
  /// This is used for efficient computation of delta3.
  VecInt2dPy m_blossomBestEdges;

  /// List of currently unused blossom numbers.
  /// if m_nVertex is 4: [4, 5, 6, 7]
  VecIntPy m_unusedBlossoms;

  /// If v is a vertex,
  /// m_dualVar[v] = 2 * u(v) where u(v) is the v's variable in the dual
  /// optimization problem (multiplication by two ensures integer values
  /// throughout the algorithm if all edge weights are integers).
  /// If b is a non-trivial blossom,
  /// m_dualVar[b] = z(b) where z(b) is b's variable in the dual optimization
  /// problem.
  /// if m_nVertex is 4 and m_maxWeight is w: [w, w, w, w, 0, 0, 0, 0]
  VecIntPy m_dualVar;

  /// If m_allowEdge[k] is true, edge k has zero slack in the optimization
  /// problem; if m_allowEdge[k] is false, the edge's slack may or may not
  /// be zero.
  VecBool m_allowEdge;
  VecIntPy m_queue; ///< Queue of newly discovered S-vertices.
};

////////////////////////////////////////////////////////////////////////////////
/// \class MeWeightMatcherImpl
/// \brief Class to determine which interior edges to remove and which boundary
/// points to split to convert a triangular mesh to quads.
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief Constructor
//------------------------------------------------------------------------------
MeWeightMatcherImpl::MeWeightMatcherImpl()
{
} // MeWeightMatcherImpl::MeWeightMatcherImpl
//------------------------------------------------------------------------------
/// \brief Initialize class member values.
/// \param[in] a_edges A vector of adjacent faces and weights.
//------------------------------------------------------------------------------
void MeWeightMatcherImpl::Init(const VecMeEdge& a_edges)
{
  m_edges = a_edges;
  m_nEdge = (int)a_edges.size();
  m_nVertex = 0;
  m_maxWeight = 0;

  // Count vertices.
  // Find the maximum edge weight.

  // Edges are numbered 0 .. (m_nEdge-1).
  // Edge endpoints are numbered 0 .. (2*m_nEdge-1), such that endpoints
  // (2*k) and (2*k+1) both belong to edge k.
  // If p is an edge endpoint,
  // m_endPoint[p] is the vertex to which endpoint p is attached.
  // Not modified by the algorithm.
  m_endPoint.reserve(m_nEdge * 2);

  for (auto& edge : m_edges)
  {
    const int& i = edge.m_f0;
    const int& j = edge.m_f1;
    XM_ASSERT(i >= 0 && j >= 0 && i != j);
    if (i >= m_nVertex)
      m_nVertex = i + 1;
    if (j >= m_nVertex)
      m_nVertex = j + 1;
    m_maxWeight = std::max(m_maxWeight, edge.m_weight);

    m_endPoint.push_back(i);
    m_endPoint.push_back(j);
  }

  // If v is a vertex,
  // m_neighbEnd[v] is the list of remote endpoints of the edges attached to v.
  // Not modified by the algorithm.
  m_neighbEnd.resize(m_nVertex, VecIntPy());
  for (int k = 0; k < (int)m_edges.size(); ++k)
  {
    const MeEdge& edge = m_edges[k];
    const int& i = edge.m_f0;
    const int& j = edge.m_f1;
    m_neighbEnd[i].push_back(2 * k + 1);
    m_neighbEnd[j].push_back(2 * k);
  }

  // If v is a vertex,
  // m_mate[v] is the remote endpoint of its matched edge, or -1 if it is single
  // (i.e. endpoint[m_mate[v]] is v's partner vertex).
  // Initially all vertices are single; updated during augmentation.
  m_mate.resize(m_nVertex, -1);

  // If b is a top-level blossom,
  // m_label[b] is 0 if b is unlabeled (free);
  //             1 if b is an S-vertex/blossom;
  //             2 if b is a T-vertex/blossom.
  // The label of a vertex is found by looking at the label of its
  // top-level containing blossom.
  // If v is a vertex inside a T-blossom,
  // m_label[v] is 2 iff v is reachable from an S-vertex outside the blossom.
  // Labels are assigned during a stage and reset after each augmentation.
  m_label.resize(2 * m_nVertex, 0);

  // If b is a labeled top-level blossom,
  // m_labelEnd[b] is the remote endpoint of the edge through which b obtained
  // its label, or -1 if b's base vertex is single.
  // If v is a vertex inside a T-blossom and m_label[v] == 2,
  // m_labelEnd[v] is the remote endpoint of the edge through which v is
  // reachable from outside the blossom.
  m_labelEnd.resize(2 * m_nVertex, -1);

  // If v is a vertex,
  // m_inBlossom[v] is the top-level blossom to which v belongs.
  // If v is a top-level vertex, v is itself a blossom (a trivial blossom)
  // and m_inBlossom[v] == v.
  // Initially all vertices are top-level trivial blossoms.
  // if m_nVertex is 4: [0, 1, 2, 3]
  m_inBlossom.resize_to_count(m_nVertex);

  // If b is a sub-blossom,
  // m_blossomParent[b] is its immediate parent (sub-)blossom.
  // If b is a top-level blossom, m_blossomParent[b] is -1.
  m_blossomParent.resize(2 * m_nVertex, -1);

  // If b is a non-trivial (sub-)blossom,
  // m_blossomChildren[b] is an ordered list of its sub-blossoms, starting with
  // the base and going round the blossom.
  // m_blossomChildren = (2 * m_nVertex) * [ None ]
  m_blossomChildren.resize(2 * m_nVertex, VecInt());

  // If b is a (sub-)blossom,
  // m_blossomBase[b] is its base VERTEX (i.e. recursive sub-blossom).
  // if m_nVertex is 4: [0, 1, 2, 3, -1, -1, -1, -1]
  m_blossomBase.resize_to_count(m_nVertex);
  m_blossomBase.resize(2 * m_nVertex, -1);

  // If b is a non-trivial (sub-)blossom,
  // m_blossomEndPts[b] is a list of endpoints on its connecting edges,
  // such that m_blossomEndPts[b][i] is the local endpoint of
  // m_blossomChildren[b][i] on the edge that connects it to
  // m_blossomChildren[b][wrap(i+1)].
  m_blossomEndPts.resize(2 * m_nVertex, VecIntPy());

  // If v is a free vertex (or an unreached vertex inside a T-blossom),
  // m_bestEdge[v] is the edge to an S-vertex with least slack,
  // or -1 if there is no such edge.
  // If b is a (possibly trivial) top-level S-blossom,
  // m_bestEdge[b] is the least-slack edge to a different S-blossom,
  // or -1 if there is no such edge.
  // This is used for efficient computation of delta2 and delta3.
  m_bestEdge.resize(2 * m_nVertex, -1);

  // If b is a non-trivial top-level S-blossom,
  // m_blossomBestEdges[b] is a list of least-slack edges to neighbouring
  // S-blossoms, or None if no such list has been computed yet.
  // This is used for efficient computation of delta3.
  // m_blossomBestEdges = (2 * m_nVertex) * [ None ]
  m_blossomBestEdges.resize(2 * m_nVertex, VecIntPy());

  // List of currently unused blossom numbers.
  // if m_nVertex is 4: [4, 5, 6, 7]
  m_unusedBlossoms.resize_to_count(m_nVertex, m_nVertex);

  // If v is a vertex,
  // m_dualVar[v] = 2 * u(v) where u(v) is the v's variable in the dual
  // optimization problem (multiplication by two ensures integer values
  // throughout the algorithm if all edge weights are integers).
  // If b is a non-trivial blossom,
  // m_dualVar[b] = z(b) where z(b) is b's variable in the dual optimization
  // problem.
  // if m_nVertex is 4 and m_maxWeight is w: [w, w, w, w, 0, 0, 0, 0]
  m_dualVar.resize(m_nVertex, m_maxWeight);
  m_dualVar.resize(m_nVertex * 2, 0);

  // If m_allowEdge[k] is true, edge k has zero slack in the optimization
  // problem; if m_allowEdge[k] is false, the edge's slack may or may not
  // be zero.
  m_allowEdge.resize(m_nEdge, false);
} // MeWeightMatcherImpl::Init
//------------------------------------------------------------------------------
/// \brief Determine which interior edges to remove and which boundary points to
/// split.
/// \param[in] a_edges The edges to operate on. Each edge identifies a pair of
/// adjacent faces and weight to prioritize against other edges.
/// \param[in] a_maxCardinality When true identify boundary points to split
/// (extend into the interior to create a new edge).
/// \return A vector v of integers where v[i] is -1 the edge wasn't removed.
/// Otherwise, the edge between face i and face v[i] should be removed.
//------------------------------------------------------------------------------
VecInt MeWeightMatcherImpl::MatchWeights(const VecMeEdge& a_edges,
                                         bool a_maxCardinality /* = false*/)
{
  Init(a_edges);

  // Deal swiftly with empty graphs.
  if (m_edges.empty())
  {
    return VecInt();
  }

  m_mate.clear();
  m_mate.resize(m_nVertex, -1);
  int BOGUS_SLACK = 0xdeadbeef;

  // Main loop: continue until no further improvement is possible.
  for (int t = 0; t < m_nVertex; ++t)
  {
    // Each iteration of this loop is a "stage".
    // A stage finds an augmenting path and uses that to improve
    // the matching.
#ifdef DEBUG
    printf("STAGE %d\n", t);
#endif

    // Remove labels from top-level blossoms/vertices.
    m_label.clear();
    m_label.resize(2 * m_nVertex, 0);

    // Forget all about least-slack edges.
    m_bestEdge.clear();
    m_bestEdge.resize(2 * m_nVertex, -1);

    for (int i = m_nVertex; i < m_nVertex * 2; ++i)
    {
      m_blossomBestEdges[i].clear();
    }

    // Loss of labeling means that we can not be sure that currently
    // allowable edges remain allowable througout this stage.
    m_allowEdge.clear();
    m_allowEdge.resize(m_nEdge, false);

    // Make m_queue empty.
    m_queue.clear();

    // Label single blossoms/vertices with S and put them in the m_queue.
    for (int v = 0; v < m_nVertex; ++v)
    {
      if (m_mate[v] == -1 && m_label[m_inBlossom[v]] == 0)
      {
        assignLabel(v, 1, -1);
      }
    }
    // Loop until we succeed in augmenting the matching.
    bool augmented = false;
    while (true)
    {
      // Each iteration of this loop is a "substage".
      // A substage tries to find an augmenting path;
      // if found, the path is used to improve the matching and
      // the stage ends. If there is no augmenting path, the
      // primal-dual method is used to pump some slack out of
      // the dual variables.
#ifdef DEBUG
      printf("SUBSTAGE\n");
#endif
      // Continue labeling until all vertices which are reachable
      // through an alternating path have got a label.
      while (!m_queue.empty() && !augmented)
      {
        // Take an S vertex from the m_queue.
        int v = m_queue.pop();
#ifdef DEBUG
        printf("POP v=%d\n", v);
#endif
        XM_ASSERT(m_label[m_inBlossom[v]] == 1);

        // Scan its neighbours:
        for (int idx = 0; idx < m_neighbEnd[v].size(); ++idx)
        {
          int& p = m_neighbEnd[v][idx];
          int k = p / 2;
          int w = m_endPoint[p];
          // w is a neighbour to v
          if (m_inBlossom[v] == m_inBlossom[w])
          {
            // this edge is internal to a blossom; ignore it
            continue;
          }
          int kslack = BOGUS_SLACK;
          bool valid_kslack = false;
          if (!m_allowEdge[k])
          {
            kslack = slack(k);
            valid_kslack = true;
            if (kslack <= 0)
            {
              // edge k has zero slack => it is allowable
              m_allowEdge[k] = true;
            }
          }

          if (m_allowEdge[k])
          {
            if (m_label[m_inBlossom[w]] == 0)
            {
              // (C1) w is a free vertex;
              // label w with T and label its mate with S (R12).
              assignLabel(w, 2, p ^ 1);
            }
            else if (m_label[m_inBlossom[w]] == 1)
            {
              // (C2) w is an S-vertex (not in the same blossom);
              // follow back-links to discover either an
              // augmenting path or a new blossom.
              int base = scanBlossom(v, w);
              if (base >= 0)
              {
                // Found a new blossom; add it to the blossom
                // bookkeeping and turn it into an S-blossom.
                addBlossom(base, k);
              }
              else
              {
                // Found an augmenting path; augment the
                // matching and end this stage.
                augmentMatching(k);
                augmented = true;
                break;
              }
            }
            else if (m_label[w] == 0)
            {
              // w is inside a T-blossom, but w itself has not
              // yet been reached from outside the blossom;
              // mark it as reached (we need this to relabel
              // during T-blossom expansion).
              XM_ASSERT(m_label[m_inBlossom[w]] == 2);
              m_label[w] = 2;
              m_labelEnd[w] = p ^ 1;
            }
          }
          else if (m_label[m_inBlossom[w]] == 1)
          {
            // keep track of the least-slack non-allowable edge to
            // a different S-blossom.
            int b = m_inBlossom[v];
            if (m_bestEdge[b] == -1 || (valid_kslack && kslack < slack(m_bestEdge[b])))
            {
              m_bestEdge[b] = k;
            }
          }
          else if (m_label[w] == 0)
          {
            // w is a free vertex (or an unreached vertex inside
            // a T-blossom) but we can not reach it yet;
            // keep track of the least-slack edge that reaches w.
            if (m_bestEdge[w] == -1 || (valid_kslack && kslack < slack(m_bestEdge[w])))
            {
              m_bestEdge[w] = k;
            }
          }
        }
      }
      if (augmented)
      {
        break;
      }

      // There is no augmenting path under these constraints;
      // compute delta and reduce slack in the optimization problem.
      // (Note that our vertex dual variables, edge slacks and delta's
      // are pre-multiplied by two.)
      int deltatype = -1;
      int NONE = -1;
      int deltaedge = NONE;
      int delta = NONE;
      int deltablossom = NONE;

      // Verify data structures for delta2/delta3 computation.
#ifdef CHECK_DELTA
      checkDelta2();
      checkDelta3();
#endif

      // Compute delta1: the minumum value of any vertex dual.
      if (!a_maxCardinality)
      {
        deltatype = 1;
        delta = m_dualVar.min(0, m_nVertex);
      }

      // Compute delta2: the minimum slack on any edge between
      // an S-vertex and a free vertex.
      int dslack = BOGUS_SLACK;
      for (int v = 0; v < m_nVertex; ++v)
      {
        if (m_label[m_inBlossom[v]] == 0 && m_bestEdge[v] != -1)
        {
          dslack = slack(m_bestEdge[v]);
          if (deltatype == -1 || dslack < delta)
          {
            delta = dslack;
            deltatype = 2;
            deltaedge = m_bestEdge[v];
          }
        }
      }

      // Compute delta3: half the minimum slack on any edge between
      // a pair of S-blossoms.
      for (int b = 0; b < 2 * m_nVertex; ++b)
      {
        if (m_blossomParent[b] == -1 && m_label[b] == 1 && m_bestEdge[b] != -1)
        {
          int kslack = slack(m_bestEdge[b]);
          XM_ASSERT((kslack % 2) == 0);
          dslack = kslack / 2;
          if (deltatype == -1 || (dslack < delta))
          {
            delta = dslack;
            deltatype = 3;
            deltaedge = m_bestEdge[b];
          }
        }
      }

      // Compute delta4: minimum z variable of any T-blossom.
      for (int b = m_nVertex; b < 2 * m_nVertex; ++b)
      {
        if (m_blossomBase[b] >= 0 && m_blossomParent[b] == -1 && m_label[b] == 2 &&
            (deltatype == -1 || m_dualVar[b] < delta))
        {
          delta = m_dualVar[b];
          deltatype = 4;
          deltablossom = b;
        }
      }

      if (deltatype == -1)
      {
        // No further improvement possible; max-cardinality optimum
        // reached. Do a final delta update to make the optimum
        // verifyable.
        XM_ASSERT(a_maxCardinality);
        deltatype = 1;
        delta = std::max(0, m_dualVar.min(0, m_nVertex));
      }

      // Update dual variables according to delta.
      for (int v = 0; v < m_nVertex; ++v)
      {
        if (m_label[m_inBlossom[v]] == 1)
        {
          // S-vertex: 2*u = 2*u - 2*delta
          m_dualVar[v] -= delta;
        }
        else if (m_label[m_inBlossom[v]] == 2)
        {
          // T-vertex: 2*u = 2*u + 2*delta
          m_dualVar[v] += delta;
        }
      }
      for (int b = m_nVertex; b < 2 * m_nVertex; ++b)
      {
        if (m_blossomBase[b] >= 0 && m_blossomParent[b] == -1)
        {
          if (m_label[b] == 1)
          {
            // top-level S-blossom: z = z + 2*delta
            m_dualVar[b] += delta;
          }
          else if (m_label[b] == 2)
          {
            // top-level T-blossom: z = z - 2*delta
            m_dualVar[b] -= delta;
          }
        }
      }

      // Take action at the point where minimum delta occurred.
#ifdef DEBUG
      printf("delta%d=%d\n", deltatype, delta);
#endif
      if (deltatype == 1)
      {
        // No further improvement possible; optimum reached.
        break;
      }
      else if (deltatype == 2)
      {
        // Use the least-slack edge to continue the search.
        m_allowEdge[deltaedge] = true;
        auto& edge = GetEdge(deltaedge);
        int i = edge.m_f0;
        int j = edge.m_f1;
        if (m_label[m_inBlossom[i]] == 0)
        {
          std::swap(i, j);
        }
        XM_ASSERT(m_label[m_inBlossom[i]] == 1);
        m_queue.push_back(i);
      }
      else if (deltatype == 3)
      {
        // Use the least-slack edge to continue the search.
        m_allowEdge[deltaedge] = true;
        auto& edge = GetEdge(deltaedge);
        int i = edge.m_f0;
        XM_ASSERT(m_label[m_inBlossom[i]] == 1);
        m_queue.push_back(i);
      }
      else if (deltatype == 4)
      {
        // Expand the least-z blossom.
        expandBlossom(deltablossom, false);
      }
      // End of a this substage.
    }

    // Stop when no more augmenting path can be found.
    if (!augmented)
    {
      break;
    }

    // End of a stage; expand all S-blossoms which have m_dualVar = 0.
    for (int b = m_nVertex; b < 2 * m_nVertex; ++b)
    {
      if (m_blossomParent[b] == -1 && m_blossomBase[b] >= 0 && m_label[b] == 1 && m_dualVar[b] == 0)
      {
        expandBlossom(b, true);
      }
    }
  }

  // Verify that we reached the optimum solution.
#ifdef CHECK_OPTIMUM
  verifyOptimum(a_maxCardinality);
#endif

  // Transform m_mate[] such that m_mate[v] is the vertex to which v is paired.
  for (int v = 0; v < m_nVertex; ++v)
  {
    if (m_mate[v] >= 0)
    {
      m_mate[v] = m_endPoint[m_mate[v]];
    }
  }
  for (int v = 0; v < m_nVertex; ++v)
  {
    XM_ASSERT(m_mate[v] == -1 || m_mate[m_mate[v]] == v);
  }
  return m_mate.ToVecInt();
} // MeWeightMatcherImpl::MatchWeights
//------------------------------------------------------------------------------
/// \brief Return 2 * slack of a_edge (does not work inside blossoms).
/// \param[in] a_k The index of the edge to calculate the slack from.
/// \return 2 * slack of the edge (does not work inside blossoms).
//------------------------------------------------------------------------------
int MeWeightMatcherImpl::slack(int a_k)
{
  if (a_k < 0)
  {
    a_k = (int)m_edges.size() + a_k;
  }
  XM_ASSERT(a_k < m_edges.size());
  const MeEdge& edge = GetEdge(a_k);
  const int& i = edge.m_f0;
  const int& j = edge.m_f1;
  const int& wt = edge.m_weight;
  return m_dualVar[i] + m_dualVar[j] - 2 * wt;
} // MeWeightMatcherImpl::slack
//------------------------------------------------------------------------------
/// \brief Generate the leaf vertices of a blossom.
/// \param[in] a_b A Blossom index.
/// \return Indices of vertices the child blossoms.
//------------------------------------------------------------------------------
VecIntPy MeWeightMatcherImpl::blossomLeaves(int a_b)
{
  VecIntPy results;
  if (a_b < m_nVertex)
  {
    results.push_back(a_b);
    return results;
  }

  results = m_blossomChildren[a_b];
#ifdef DEBUG
  printf("blossomLeaves(%d) - blossomchilds[b] = %s\n", a_b, TS_AS_STRING(results.ToVecInt()));
#endif

  while (results.max() > m_nVertex)
  {
    // replace big values
    VecIntPy newresults;
    for (int idx = 0; idx < results.size(); ++idx)
    {
      int elem = results[idx];
      if (elem < m_nVertex)
      {
        newresults.push_back(elem);
      }
      else
      {
        newresults.append(m_blossomChildren[elem]);
      }
    }
    results = newresults;
  }
  return results;
} // MeWeightMatcherImpl::blossomLeaves
//------------------------------------------------------------------------------
/// \brief Assign label t to the top-level blossom containing vertex w and
/// record the fact that w was reached through the edge with remote endpoint p.
/// \param[in] a_w A vertex index.
/// \param[in] a_t The label to assign.
/// \param[in] a_p The remote end point.
//------------------------------------------------------------------------------
void MeWeightMatcherImpl::assignLabel(int a_w, int a_t, int a_p)
{
#ifdef DEBUG
  printf("assignLabel(%d,%d,%d)\n", a_w, a_t, a_p);
#endif
  int b = m_inBlossom[a_w];
  XM_ASSERT(m_label[a_w] == 0 && m_label[b] == 0);
  m_label[a_w] = m_label[b] = a_t;
  m_labelEnd[a_w] = m_labelEnd[b] = a_p;
  m_bestEdge[a_w] = m_bestEdge[b] = -1;
  if (a_t == 1)
  {
    // b became an S-vertex/blossom; add it(s vertices) to the m_queue.
    VecIntPy bleaves = blossomLeaves(b);
    m_queue.append(bleaves);
#ifdef DEBUG
    printf("PUSH %s\n", TS_AS_STRING(bleaves.ToVecInt()));
#endif
  }
  else if (a_t == 2)
  {
    // b became a T-vertex/blossom; assign label S to its mate.
    // (If b is a non-trivial blossom, its base is the only vertex
    // with an external mate.)
    int base = m_blossomBase[b];
    int mateBase = m_mate[base];
    XM_ASSERT(mateBase >= 0);
    assignLabel(m_endPoint[mateBase], 1, mateBase ^ 1);
  }
} // MeWeightMatcherImpl::assignLabel
//------------------------------------------------------------------------------
/// \brief Trace back from vertices v and w to discover either a new blossom
/// or an augmenting path.
/// \param[in] a_v A vertex index
/// \param[in] a_w Another vertex index
/// \return The base vertex of the new blossom or -1.
//------------------------------------------------------------------------------
int MeWeightMatcherImpl::scanBlossom(int a_v, int a_w)
{
#ifdef DEBUG
  printf("scanBlossom(%d, %d)\n", a_v, a_w);
#endif
  // Trace back from v and w, placing beadcrumbs as we go.
  VecIntPy path;
  int base = -1;
  while (a_v != -1 || a_w != -1)
  {
    // Look for a breadcrumb in v's blossom or put a new breadcrumb.
    int b = m_inBlossom[a_v];
    if (m_label[b] & 4)
    {
      base = m_blossomBase[b];
      break;
    }
    XM_ASSERT(m_label[b] == 1);
    path.push_back(b);
    m_label[b] = 5;
    // Trace one step back.
    XM_ASSERT(m_labelEnd[b] == m_mate[m_blossomBase[b]]);
    if (m_labelEnd[b] == -1)
    {
      // The base of blossom b is single; stop stracing this path.
      a_v = -1;
    }
    else
    {
      a_v = m_endPoint[m_labelEnd[b]];
      b = m_inBlossom[a_v];
      XM_ASSERT(m_label[b] == 2);
      // b is a T-blossom; trace one more step back.
      XM_ASSERT(m_labelEnd[b] >= 0);
      a_v = m_endPoint[m_labelEnd[b]];
    }
    // Swap v and w so that we alternate between both paths.
    if (a_w != -1)
    {
      std::swap(a_v, a_w);
    }
  }
  // Remove breadcrumbs.
  for (int idx = 0; idx < path.size(); ++idx)
  {
    int b = path[idx];
    m_label[b] = 1;
  }
  // Return base vertex, if we found one.
  return base;
} // MeWeightMatcherImpl::scanBlossom
//------------------------------------------------------------------------------
/// \brief Construct a new blossom with given base, containing edge k which
/// connects a pair of S vertices. Label the new blossom as S; set its dual
/// variable to zero; relabel its T-vertices to S and add them to the m_queue.
/// \param[in] a_base A vertex intex
/// \param[in] a_k An edge index
/// \return The base vertex of the new blossom or -1.
//------------------------------------------------------------------------------
void MeWeightMatcherImpl::addBlossom(int a_base, int a_k)
{
  const MeEdge& edge = GetEdge(a_k);
  int v = edge.m_f0;
  int w = edge.m_f1;
  int bb = m_inBlossom[a_base];
  int bv = m_inBlossom[v];
  int bw = m_inBlossom[w];
  // Create blossom.
  int b = m_unusedBlossoms.pop();
#ifdef DEBUG
  printf("addBlossom(%d,%d) (v=%d w=%d) -> %d\n", a_base, a_k, v, w, b);
#endif
  m_blossomBase[b] = a_base;
  m_blossomParent[b] = -1;
  m_blossomParent[bb] = b;
  // Make list of sub-blossoms and their interconnecting edge endpoints.
  m_blossomChildren[b].clear();
  VecIntPy& path = m_blossomChildren[b];
#ifdef DEBUG
  printf("blossomchilds[%d] cleared\n", b);
#endif
  m_blossomEndPts[b].clear();
  VecIntPy& endps = m_blossomEndPts[b];
  // Trace back from v to base.
  while (bv != bb)
  {
    // Add bv to the new blossom.
    m_blossomParent[bv] = b;
    path.push_back(bv);
    endps.push_back(m_labelEnd[bv]);
    XM_ASSERT(m_label[bv] == 2 ||
              (m_label[bv] == 1 && m_labelEnd[bv] == m_mate[m_blossomBase[bv]]));
    // Trace one step back.
    XM_ASSERT(m_labelEnd[bv] >= 0);
    v = m_endPoint[m_labelEnd[bv]];
    bv = m_inBlossom[v];
  }
  // Reverse lists, add endpoint that connects the pair of S vertices.
  path.push_back(bb);
  path.reverse();
  endps.reverse();
  endps.push_back(2 * a_k);
  // Trace back from w to base.
  while (bw != bb)
  {
    // Add bw to the new blossom.
    m_blossomParent[bw] = b;
    path.push_back(bw);
    endps.push_back(m_labelEnd[bw] ^ 1);
    XM_ASSERT(m_label[bw] == 2 ||
              (m_label[bw] == 1 && m_labelEnd[bw] == m_mate[m_blossomBase[bw]]));
    // Trace one step back;
    XM_ASSERT(m_labelEnd[bw] >= 0);
    w = m_endPoint[m_labelEnd[bw]];
    bw = m_inBlossom[w];
  }
  // Set label to S.
  XM_ASSERT(m_label[bb] == 1);
  m_label[b] = 1;
  m_labelEnd[b] = m_labelEnd[bb];
  // Set dual variable to zero.
  m_dualVar[b] = 0;
  // Relabel vertices.

  VecIntPy bleaves = blossomLeaves(b);
  for (int idx = 0; idx < bleaves.size(); ++idx)
  {
    v = bleaves[idx];
    if (m_label[m_inBlossom[v]] == 2) // is T-vertex
    {
      // This T-vertex now turns into an S-vertex because it becomes
      // part of an S-blossom; add it to the m_queue.
      m_queue.push_back(v);
    }
    m_inBlossom[v] = b;
  }

  // Compute m_blossomBestEdges[b].
  VecIntPy bestedgeto;
  bestedgeto.resize(2 * m_nVertex, -1);
  for (int idx = 0; idx < path.size(); ++idx)
  {
    int bv = path[idx];
    VecInt2dPy nblists;
    if (m_blossomBestEdges[bv].empty())
    {
      // This subblossom does not have a list of least-slack edges;
      // get the information from the vertices.
      VecIntPy nblist;
      for (int idx2 = 0; idx2 < m_neighbEnd[v].size(); ++idx2)
      {
        int p = m_neighbEnd[v][idx2];
        nblist.push_back(p / 2);
      }
      nblists.resize(blossomLeaves(bv).size(), nblist);
    }
    else
    {
      // Walk this subblossom's least-slack edges.
      nblists.resize(1);
      nblists[0] = m_blossomBestEdges[bv];
    }
    for (int idx3 = 0; idx3 < nblists.size(); ++idx3)
    {
      auto& nblist = nblists[idx3];
      for (int idx4 = 0; idx4 < nblist.size(); ++idx4)
      {
        int k = nblist[idx4];
        const MeEdge& edge = GetEdge(k);
        int i = edge.m_f0;
        int j = edge.m_f1;
        if (m_inBlossom[j] == b)
        {
          std::swap(i, j);
        }
        int bj = m_inBlossom[j];
        if (bj != b && m_label[bj] == 1 &&
            (bestedgeto[bj] == -1 || slack(k) < slack(bestedgeto[bj])))
        {
          bestedgeto[bj] = k;
        }
      }
    }
    // Forget about least-slack edges of the subblossom.
    m_blossomBestEdges[bv].clear();
    m_bestEdge[bv] = -1;
  }
  VecIntPy bk;
  for (int idx = 0; idx < bestedgeto.size(); ++idx)
  {
    int k = bestedgeto[idx];
    if (k != -1)
    {
      bk.push_back(k);
    }
  }
  m_blossomBestEdges[b] = bk; // [ k for k in bestedgeto if k != -1 ]
  // Select m_bestEdge[b]
  m_bestEdge[b] = -1;
  for (int idx = 0; idx < bk.size(); ++idx)
  {
    int k = bk[idx];
    if (m_bestEdge[b] == -1 || slack(k) < slack(m_bestEdge[b]))
    {
      m_bestEdge[b] = k;
    }
  }
#ifdef DEBUG
  printf("blossomchilds[%d]=%s\n", b, TS_AS_STRING(m_blossomChildren[b].ToVecInt()));
#endif
} // MeWeightMatcherImpl::addBlossom
//------------------------------------------------------------------------------
/// \brief Expand the given top-level blossom.
/// \param[in] a_b The blossom.
/// \param[in] a_endStage Flag to determine if recursing.
//------------------------------------------------------------------------------
void MeWeightMatcherImpl::expandBlossom(int a_b, bool a_endStage)
{
#ifdef DEBUG
  printf("expandBlossom(%d,%d) %s\n", a_b, a_endStage,
         TS_AS_STRING(m_blossomChildren[a_b].ToVecInt()));
#endif
  // Convert sub-blossoms into top-level blossoms.
  for (int idx = 0; idx < m_blossomChildren[a_b].size(); ++idx)
  {
    int s = m_blossomChildren[a_b][idx];
    m_blossomParent[s] = -1;
    if (s < m_nVertex)
    {
      m_inBlossom[s] = s;
    }
    else if (a_endStage && m_dualVar[s] == 0)
    {
      // Recursively expand this sub-blossom.
      expandBlossom(s, a_endStage);
    }
    else
    {
      VecIntPy bleaves = blossomLeaves(s);
      for (int idx2 = 0; idx2 < bleaves.size(); ++idx2)
      {
        int v = bleaves[idx2];
        m_inBlossom[v] = s;
      }
    }
  }

  // If we expand a T-blossom during a stage, its sub-blossoms must be
  // relabeled.
  if (!a_endStage && m_label[a_b] == 2)
  {
    // Start at the sub-blossom through which the expanding
    // blossom obtained its label, and relabel sub-blossoms untili
    // we reach the base.
    // Figure out through which sub-blossom the expanding blossom
    // obtained its label initially.
    XM_ASSERT(m_labelEnd[a_b] >= 0);
    int entrychild = m_inBlossom[m_endPoint[m_labelEnd[a_b] ^ 1]];
    // Decide in which direction we will go round the blossom.
    int j = m_blossomChildren[a_b].index(entrychild);
    // Assume start index is even; go backward.
    int jstep = -1;
    int endptrick = 1;
    if (j & 1)
    {
      // Start index is odd; go forward and wrap.
      j -= (int)m_blossomChildren[a_b].size();
      jstep = 1;
      endptrick = 0;
    }

    // Move along the blossom until we get to the base.
    int p = m_labelEnd[a_b];
    while (j != 0)
    {
      // Relabel the T-sub-blossom.
      m_label[m_endPoint[p ^ 1]] = 0;
      m_label[m_endPoint[m_blossomEndPts[a_b][j - endptrick] ^ endptrick ^ 1]] = 0;
      assignLabel(m_endPoint[p ^ 1], 2, p);
      // Step to the next S-sub-blossom and note its forward endpoint.
      m_allowEdge[m_blossomEndPts[a_b][j - endptrick] / 2] = true;
      j += jstep;
      p = m_blossomEndPts[a_b][j - endptrick] ^ endptrick;
      // Step to the next T-sub-blossom.
      m_allowEdge[p / 2] = true;
      j += jstep;
    }
    // Relabel the base T-sub-blossom WITHOUT stepping through to
    // its mate (so don't call assignLabel).
    int bv = m_blossomChildren[a_b][j];
    m_label[m_endPoint[p ^ 1]] = m_label[bv] = 2;
    m_labelEnd[m_endPoint[p ^ 1]] = m_labelEnd[bv] = p;
    m_bestEdge[bv] = -1;
    // Continue along the blossom until we get back to entrychild.
    j += jstep;
    while (m_blossomChildren[a_b][j] != entrychild)
    {
      // Examine the vertices of the sub-blossom to see whether
      // it is reachable from a neighbouring S-vertex outside the
      // expanding blossom.
      bv = m_blossomChildren[a_b][j];
      if (m_label[bv] == 1)
      {
        // This sub-blossom just got label S through one of its
        // neighbours; leave it.
        j += jstep;
        continue;
      }
      VecIntPy bleaves = blossomLeaves(bv);
      XM_ASSERT(!bleaves.empty());
      int v; // get where zero or last item in bleaves
      for (int idx2 = 0; idx2 < bleaves.size(); ++idx2)
      {
        v = bleaves[idx2];
        if (m_label[v] != 0)
        {
          break;
        }
      }
      // If the sub-blossom contains a reachable vertex, assign
      // label T to the sub-blossom.
      if (m_label[v] != 0)
      {
        XM_ASSERT(m_label[v] == 2);
        XM_ASSERT(m_inBlossom[v] == bv);
        m_label[v] = 0;
        m_label[m_endPoint[m_mate[m_blossomBase[bv]]]] = 0;
        assignLabel(v, 2, m_labelEnd[v]);
      }
      j += jstep;
    }
  }
  // Recycle the blossom number.
  m_label[a_b] = m_labelEnd[a_b] = -1;
  m_blossomChildren[a_b].clear();
#ifdef DEBUG
  printf("clear blossomchilds[%d]\n", a_b);
#endif
  m_blossomEndPts[a_b].clear();
  m_blossomBase[a_b] = -1;
  m_blossomBestEdges[a_b].clear();
  m_bestEdge[a_b] = -1;
  m_unusedBlossoms.push_back(a_b);
} // MeWeightMatcherImpl::expandBlossom
//------------------------------------------------------------------------------
/// \brief Swap matched/unmatched edges over an alternating path through
/// blossom b between vertex v and the base vertex. Keep blossom bookkeeping
/// consistent.
/// \param[in] a_b The blossom.
/// \param[in] a_v The vertex.
//------------------------------------------------------------------------------
void MeWeightMatcherImpl::augmentBlossom(int a_b, int a_v)
{
#ifdef DEBUG
  printf("augmentBlossom(%d,%d)\n", a_b, a_v);
#endif
  // Bubble up through the blossom tree from vertex v to an immediate
  // sub-blossom of b.
  int t = a_v;
  while (m_blossomParent[t] != a_b)
  {
    t = m_blossomParent[t];
  }
  // Recursively deal with the first sub-blossom.
  if (t >= m_nVertex)
  {
    augmentBlossom(t, a_v);
  }
  // Decide in which direction we will go round the blossom.
  int i = m_blossomChildren[a_b].index(t);
  int j = i;
  // Assume start index is even; go backward.
  int jstep = -1;
  int endptrick = 1;
  if (i & 1)
  {
    // Start index is odd; go forward and wrap.
    j -= (int)m_blossomChildren[a_b].size();
    jstep = 1;
    endptrick = 0;
  }
  // Move along the blossom until we get to the base.
  while (j != 0)
  {
    // Step to the next sub-blossom and augment it recursively.
    j += jstep;
    t = m_blossomChildren[a_b][j];
    int p = m_blossomEndPts[a_b][j - endptrick] ^ endptrick;
    if (t >= m_nVertex)
    {
#ifdef DEBUG
      printf("augmentBlossom 1\n");
#endif
      augmentBlossom(t, m_endPoint[p]);
    }
    // Step to the next sub-blossom and augment it recursively.
    j += jstep;
    t = m_blossomChildren[a_b][j];
    if (t >= m_nVertex)
    {
#ifdef DEBUG
      printf("augmentBlossom 2\n");
#endif
      augmentBlossom(t, m_endPoint[p ^ 1]);
    }
    // Match the edge connecting those sub-blossoms.
    m_mate[m_endPoint[p]] = p ^ 1;
    m_mate[m_endPoint[p ^ 1]] = p;
#ifdef DEBUG
    printf("PAIR %d %d (k=%d)\n", m_endPoint[p], m_endPoint[p ^ 1], p / 2);
#endif
  }
  // Rotate the list of sub-blossoms to put the new base at the front.
  m_blossomChildren[a_b].rotate(i);
#ifdef DEBUG
  printf("rotate blossomchilds[%d]=%s\n", a_b, TS_AS_STRING(m_blossomChildren[a_b].ToVecInt()));
#endif
  m_blossomEndPts[a_b].rotate(i);
  m_blossomBase[a_b] = m_blossomBase[m_blossomChildren[a_b][0]];
  XM_ASSERT(m_blossomBase[a_b] == a_v);
} // MeWeightMatcherImpl::augmentBlossom
//------------------------------------------------------------------------------
/// \brief Swap matched/unmatched edges over an alternating path between
/// two single vertices. The augmenting path runs through edge k which connects
/// a pair of S vertices
/// \param[in] a_k an MeEdge
//------------------------------------------------------------------------------
void MeWeightMatcherImpl::augmentMatching(int a_k)
{
  const MeEdge& edge = GetEdge(a_k);
  int v = edge.m_f0;
  int w = edge.m_f1;
#ifdef DEBUG
  printf("augmentMatching(%d) (v=%d w=%d)\n", a_k, v, w);
  printf("PAIR %d %d (k=%d)\n", v, w, a_k);
#endif
  VecInt2d vwkTemp = {{v, 2 * a_k + 1}, {w, 2 * a_k}};
  VecInt2dPy vwk(vwkTemp);
  for (int idx = 0; idx < vwk.size(); ++idx)
  {
    auto& sp = vwk[idx];
    int s = sp[0];
    int p = sp[1];

    // Match vertex s to remote endpoint p. Then trace back from s
    // until we find a single vertex, swapping matched and unmatched
    // edges as we go.
    while (true)
    {
      int bs = m_inBlossom[s];
      XM_ASSERT(m_label[bs] == 1);
      XM_ASSERT(m_labelEnd[bs] == m_mate[m_blossomBase[bs]]);

      // Augment through the S-blossom from s to base.
      if (bs >= m_nVertex)
      {
        augmentBlossom(bs, s);
      }
      // Update m_mate[s]
      m_mate[s] = p;
      // Trace one step back.
      if (m_labelEnd[bs] == -1)
      {
        // Reached single vertex; stop.
        break;
      }
      int t = m_endPoint[m_labelEnd[bs]];
      int bt = m_inBlossom[t];
      XM_ASSERT(m_label[bt] == 2);
      // Trace one step back.
      XM_ASSERT(m_labelEnd[bt] >= 0);
      s = m_endPoint[m_labelEnd[bt]];
      int j = m_endPoint[m_labelEnd[bt] ^ 1];
      // Augment through the T-blossom from j to base.
      XM_ASSERT(m_blossomBase[bt] == t);
      if (bt >= m_nVertex)
      {
        augmentBlossom(bt, j);
      }
      // Update m_mate[j]
      m_mate[j] = m_labelEnd[bt];
      // Keep the opposite endpoint;
      // it will be assigned to m_mate[s] in the next step.
      p = m_labelEnd[bt] ^ 1;
#ifdef DEBUG
      printf("PAIR %d %d (k=%d)\n", s, t, p / 2);
#endif
    }
  }
} // MeWeightMatcherImpl::augmentMatching
//------------------------------------------------------------------------------
/// \brief Verify that the optimum solution has been reached.
/// \param[in] a_maxCardinality The maximum cardinality.
//------------------------------------------------------------------------------
void MeWeightMatcherImpl::verifyOptimum(bool a_maxCardinality)
{
  int vdualoffset = 0;
  if (a_maxCardinality)
  {
    // Vertices may have negative dual;
    // find a constant non-negative number to add to all vertex duals.
    vdualoffset = std::max(0, -m_dualVar.min(0, m_nVertex));
#if _DEBUG
    XM_ASSERT(m_dualVar.min(0, m_nVertex) + vdualoffset >= 0);
    XM_ASSERT(m_dualVar.min(m_nVertex, m_dualVar.size()) >= 0);
#endif
  }
  // 0. all dual variables are non-negative
  // 0. all edges have non-negative slack and
  // 1. all matched edges have zero slack;
  for (int k = 0; k < m_edges.size(); ++k)
  {
    const MeEdge& edge = GetEdge(k);
    int i = edge.m_f0;
    int j = edge.m_f1;
    int wt = edge.m_weight;
    int s = m_dualVar[i] + m_dualVar[j] - 2 * wt;
    VecIntPy iblossoms(1, i);
    VecIntPy jblossoms(1, j);
    while (m_blossomParent[iblossoms.back()] != -1)
    {
      iblossoms.push_back(m_blossomParent[iblossoms.back()]);
    }
    while (m_blossomParent[jblossoms.back()] != -1)
    {
      jblossoms.push_back(m_blossomParent[jblossoms.back()]);
    }
    iblossoms.reverse();
    jblossoms.reverse();
    int ijmin = (int)std::min(iblossoms.size(), jblossoms.size());
    for (int ij = 0; ij < ijmin; ++ij)
    {
      int bi = iblossoms[ij];
      int bj = jblossoms[ij];
      if (bi != bj)
      {
        break;
      }
      s += 2 * m_dualVar[bi];
    }
    XM_ASSERT(s >= 0);
    int mate_i = m_mate[i];
    mate_i = mate_i == -1 ? -1 : mate_i / 2;
    int mate_j = m_mate[j];
    mate_j = mate_j == -1 ? -1 : mate_j / 2;
    if (mate_i == k || mate_j == k)
    {
      XM_ASSERT(mate_i == k && mate_j == k);
      XM_ASSERT(s == 0);
    }
  }
  // 2. all single vertices have zero dual value;
  for (int v = 0; v < m_nVertex; ++v)
  {
    XM_ASSERT(m_mate[v] >= 0 || m_dualVar[v] + vdualoffset == 0);
  }
  // 3. all blossoms with positive dual value are full.
  for (int b = m_nVertex; b < 2 * m_nVertex; ++b)
  {
#if _DEBUG
    if (m_blossomBase[b] >= 0 && m_dualVar[b] > 0)
    {
      VecIntPy& bendps = m_blossomEndPts[b];
      XM_ASSERT(bendps.size() % 2 == 1);
      for (int ip = 1; ip < bendps.size(); ip += 2)
      {
        int p = bendps[ip];
        XM_ASSERT(m_mate[m_endPoint[p]] == (p ^ 1));
        XM_ASSERT(m_mate[m_endPoint[p ^ 1]] == p);
      }
    }
#endif
  }
} // MeWeightMatcherImpl::verifyOptimum
//------------------------------------------------------------------------------
/// \brief Check optimized delta2 against a trivial computation.
//------------------------------------------------------------------------------
void MeWeightMatcherImpl::checkDelta2()
{
  for (int v = 0; v < m_nVertex; ++v)
  {
    if (m_label[m_inBlossom[v]] == 0)
    {
      int bd = 0xdeadc0de; // doesn't get used until after first time through loop
      int* pbd = NULL;
      int bk = -1;
      for (int idx = 0; idx < m_neighbEnd[v].size(); ++idx)
      {
        int p = m_neighbEnd[v][idx];
        int k = p / 2;
        int w = m_endPoint[p];
        if (m_label[m_inBlossom[w]] == 1)
        {
          int d = slack(k);
          if (bk == -1 || (pbd != NULL && d < *pbd))
          {
            bk = k;
            pbd = &bd;
            bd = d;
          }
        }
      }
#ifdef DEBUG
      int bev = m_bestEdge[v];
      int sbev = slack(bev);
      if ((bev != -1 || bk != -1) && (bev == -1 || (pbd != NULL && *pbd != sbev)))
      {
        printf("v=%d bk=%d bd=%s bestedge=%d slack=%d\n", v, bk, pbd ? TS_AS_STRING(*pbd) : "NULL",
               bev, sbev);
      }
      XM_ASSERT((bk == -1 && bev == -1) || (bev != -1 && bd == sbev));
#endif
    }
  }
} // MeWeightMatcherImpl::checkDelta2
//------------------------------------------------------------------------------
/// \brief Check optimized delta3 against a trivial computation.
//------------------------------------------------------------------------------
void MeWeightMatcherImpl::checkDelta3()
{
#ifdef DEBUG
  printf("checkDelta3\n");
#endif
  int bk = -1;
  int bd = 0xbaadf00d; // doesn't get used until after first time through loop
  int* pbd = NULL;
  int tbk = -1;
  int tbd = 0xbaadf00d; // doesn't get used until after first time through loop
  int* ptbd = NULL;
  for (int b = 0; b < 2 * m_nVertex; ++b)
  {
    if (m_blossomParent[b] == -1 && m_label[b] == 1)
    {
      VecIntPy bleaves = blossomLeaves(b);
      for (int idx = 0; idx < bleaves.size(); ++idx)
      {
        int v = bleaves[idx];
        for (int idx2 = 0; idx2 < m_neighbEnd[v].size(); ++idx2)
        {
          int p = m_neighbEnd[v][idx2];
          int k = p / 2;
          int w = m_endPoint[p];
          if (m_inBlossom[w] != b && m_label[m_inBlossom[w]] == 1)
          {
            int d = slack(k);
            if (bk == -1 || (pbd != NULL && d < *pbd))
            {
              bk = k;
              pbd = &bd;
              bd = d;
#ifdef DEBUG
              printf("changed bd b=%d,v=%d,p=%d,k=%d,bd=%d\n", b, v, p, k, bd);
#endif
            }
          }
        }
      }
      if (m_bestEdge[b] != -1)
      {
#if _DEBUG
        const MeEdge& edge = GetEdge(m_bestEdge[b]);
        int i = edge.m_f0;
        int j = edge.m_f1;
        XM_ASSERT(m_inBlossom[i] == b || m_inBlossom[j] == b);
        XM_ASSERT(m_inBlossom[i] != b || m_inBlossom[j] != b);
        XM_ASSERT(m_label[m_inBlossom[i]] == 1 && m_label[m_inBlossom[j]] == 1);
#endif
        if (tbk == -1 || (ptbd != NULL && slack(m_bestEdge[b]) < *ptbd))
        {
          tbk = m_bestEdge[b];
          ptbd = &tbd;
          tbd = slack(m_bestEdge[b]);
#ifdef DEBUG
          printf("changed tbd b=%d,tbk=%d,tbd=%d\n", b, tbk, tbd);
#endif
        }
      }
    }
  }
#ifdef DEBUG
  if (bd != tbd)
  {
    printf("bk=%d tbk=%d bd=%s tbd=%s\n", bk, tbk, pbd ? TS_AS_STRING(*pbd) : "NULL",
           ptbd ? TS_AS_STRING(*ptbd) : "NULL");
  }
#endif
  XM_ASSERT(bd == tbd);
} // MeWeightMatcherImpl::checkDelta3
//------------------------------------------------------------------------------
/// \brief Get an edge at a specified index location where the specified index
/// could be negative.
/// \param[in] a_idx The specified index. When negative the index retrieves a
/// value indexed from the end of the vector. The last item is at index
/// location -1.
/// \return The MeEdge at the indexed location.
//------------------------------------------------------------------------------
const MeEdge& MeWeightMatcherImpl::GetEdge(int a_idx)
{
  if (a_idx < 0)
  {
    a_idx = (int)m_edges.size() + a_idx;
  }
  XM_ASSERT(a_idx >= 0 && a_idx < (int)m_edges.size());
  return m_edges[a_idx];
} // MeWeightMatcherImpl::GetEdge

} // namespace

////////////////////////////////////////////////////////////////////////////////
/// \class MeWeightMatcher
/// \see MeWeightMatcherImpl
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief Create new weight matcher.
/// \return The new MeWeightMatcher.
//------------------------------------------------------------------------------
BSHP<MeWeightMatcher> MeWeightMatcher::New()
{
  BSHP<MeWeightMatcher> wm(new MeWeightMatcherImpl);
  return wm;
} // MeWeightMatcher::New
//------------------------------------------------------------------------------
/// \brief Constructor.
//------------------------------------------------------------------------------
MeWeightMatcher::MeWeightMatcher()
{
} // MeWeightMatcher::MeWeightMatcher
//------------------------------------------------------------------------------
/// \brief Destructor.
//------------------------------------------------------------------------------
MeWeightMatcher::~MeWeightMatcher()
{
} // MeWeightMatcher::~MeWeightMatcher

} // namespace xms

#if CXX_TEST
////////////////////////////////////////////////////////////////////////////////
// UNIT TESTS
////////////////////////////////////////////////////////////////////////////////

#include <xmsmesh/meshing/detail/MeWeightMatcher.t.h>

#include <xmscore/testing/TestTools.h>
#include <xmsinterp/triangulate/TrTriangulatorPoints.h>

//----- Namespace declaration --------------------------------------------------

using namespace xms;

namespace
{
/// Helper to run matcher tests.
#define TS_ASSERT_MATCH_WEIGHTS(a_expected, a_edges, a_cardinality) \
  _TS_ASSERT_MATCH_WEIGHTS(__FILE__, __LINE__, a_expected, a_edges, a_cardinality)
/// Helper to run matcher tests.
#define _TS_ASSERT_MATCH_WEIGHTS(a_file, a_line, a_expected, a_edges, a_cardinality) \
  iMatchWeights(a_file, a_line, a_expected, a_edges, a_cardinality)

//------------------------------------------------------------------------------
/// \brief Helper function to create a weight matcher and run a test.
/// \param[in] a_file The file for the test (this source file).
/// \param[in] a_line The line iMatchWeights was called from.
/// \param[in] a_expected The expected matcher output.
/// \param[in] a_edges The matcher edge input.
/// \param[in] a_useMaxCardinality Cardinality value passed to matcher.
//------------------------------------------------------------------------------
void iMatchWeights(const char* a_file,
                   int a_line,
                   const VecInt& a_expected,
                   const VecInt2d& a_edges,
                   bool a_useMaxCardinality = false)
{
  VecMeEdge edges;
  for (auto& edgeVec : a_edges)
  {
    MeEdge edge(edgeVec[0], edgeVec[1], edgeVec[2]);
    edges.push_back(edge);
  }

  MeWeightMatcherImpl weightMatcher;
  VecInt matchWeights = weightMatcher.MatchWeights(edges, a_useMaxCardinality);
  _TS_ASSERT_EQUALS(a_file, a_line, a_expected, matchWeights);
} // iTestMatchWeights

} // namespace
////////////////////////////////////////////////////////////////////////////////
/// \class MeWeightMatcherUnitTests
/// \brief Tests for MeWeightMatcher.
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief Test VecIntPy class.
//------------------------------------------------------------------------------
void MeWeightMatcherUnitTests::testPythonVectors()
{
  VecInt vals = {1, 5, -1, 1, 2, 3, 4};
  VecIntPy values(vals);
  TS_ASSERT_EQUALS(4, values.back());
  TS_ASSERT_EQUALS(3, values[-2]);

  auto min = values.min();
  TS_ASSERT_EQUALS(-1, min);

  min = values.min(1, 3);
  TS_ASSERT_EQUALS(-1, min);

  auto max = values.max();
  TS_ASSERT_EQUALS(5, max);

  max = values.max(2, 6);
  TS_ASSERT_EQUALS(3, max);

  values.rotate(4);
  VecInt expected = {2, 3, 4, 1, 5, -1, 1};
  TS_ASSERT_EQUALS(expected, values.ToVecInt());

  auto index = values.index(1);
  TS_ASSERT_EQUALS(3, index);

  values.resize_to_count(7);
  expected = {0, 1, 2, 3, 4, 5, 6};
  TS_ASSERT_EQUALS(expected, values.ToVecInt());

  values.resize_to_count(4, 5);
  expected = {5, 6, 7, 8};
  TS_ASSERT_EQUALS(expected, values.ToVecInt());

  values.reverse();
  expected = {8, 7, 6, 5};
  TS_ASSERT_EQUALS(expected, values.ToVecInt());
} // MeWeightMatcherUnitTests::testPythonVectors
//------------------------------------------------------------------------------
/// \brief Test empty input graph.
//------------------------------------------------------------------------------
void MeWeightMatcherUnitTests::test10_empty()
{
#ifdef DEBUG
  printf("test10_empty\n\n");
#endif
  TS_ASSERT_MATCH_WEIGHTS({}, {}, false);
} // MeWeightMatcherUnitTests::test10_empty
//------------------------------------------------------------------------------
/// \brief Test single edge.
//------------------------------------------------------------------------------
void MeWeightMatcherUnitTests::test11_singleedge()
{
#ifdef DEBUG
  printf("test11_singleedge\n\n");
#endif
  VecInt2d edges = {{0, 1, 1}};
  VecInt expected = {1, 0};
  TS_ASSERT_MATCH_WEIGHTS(expected, edges, false);
} // MeWeightMatcherUnitTests::test11_singleedge
//------------------------------------------------------------------------------
/// \brief Test with two edges.
//------------------------------------------------------------------------------
void MeWeightMatcherUnitTests::test12()
{
#ifdef DEBUG
  printf("test12\n\n");
#endif
  VecInt2d edges = {{1, 2, 10}, {2, 3, 11}};
  VecInt expected = {-1, -1, 3, 2};
  TS_ASSERT_MATCH_WEIGHTS(expected, edges, false);
} // MeWeightMatcherUnitTests::test12
//------------------------------------------------------------------------------
/// \brief Test with three edges.
//------------------------------------------------------------------------------
void MeWeightMatcherUnitTests::test13()
{
#ifdef DEBUG
  printf("test13\n\n");
#endif
  VecInt2d edges = {{1, 2, 5}, {2, 3, 11}, {3, 4, 5}};
  VecInt expected = {-1, -1, 3, 2, -1};
  TS_ASSERT_MATCH_WEIGHTS(expected, edges, false);
} // MeWeightMatcherUnitTests::test13
//------------------------------------------------------------------------------
/// \brief Test maximum cardinality.
//------------------------------------------------------------------------------
void MeWeightMatcherUnitTests::test14_maxcard()
{
#ifdef DEBUG
  printf("test14_maxcard\n\n");
#endif
  VecInt2d edges = {{1, 2, 5}, {2, 3, 11}, {3, 4, 5}};
  VecInt expected = {-1, 2, 1, 4, 3};
  TS_ASSERT_MATCH_WEIGHTS(expected, edges, true);
} // MeWeightMatcherUnitTests::test14_maxcard
//------------------------------------------------------------------------------
/// \brief Test negative weights.
//------------------------------------------------------------------------------
void MeWeightMatcherUnitTests::test16_negative()
{
#ifdef DEBUG
  printf("test16_negative\n\n");
#endif
  VecInt2d edges = {{1, 2, 2}, {1, 3, -2}, {2, 3, 1}, {2, 4, -1}, {3, 4, -6}};
  VecInt expected = {-1, 2, 1, -1, -1};
  TS_ASSERT_MATCH_WEIGHTS(expected, edges, false);
#ifdef DEBUG
  printf("test16_negative-true\n\n");
#endif
  expected = {-1, 3, 4, 1, 2};
  TS_ASSERT_MATCH_WEIGHTS(expected, edges, true);
} // MeWeightMatcherUnitTests::test16_negative
//------------------------------------------------------------------------------
/// \brief Test create S-blossom and use it for augmentation.
//------------------------------------------------------------------------------
void MeWeightMatcherUnitTests::test20_sblossom()
{
#ifdef DEBUG
  printf("test20_sblossom-1\n\n");
#endif
  VecInt2d edges = {{1, 2, 8}, {1, 3, 9}, {2, 3, 10}, {3, 4, 7}};
  VecInt expected = {-1, 2, 1, 4, 3};
  TS_ASSERT_MATCH_WEIGHTS(expected, edges, false);
#ifdef DEBUG
  printf("test20_sblossom-2\n\n");
#endif
  edges = {{1, 2, 8}, {1, 3, 9}, {2, 3, 10}, {3, 4, 7}, {1, 6, 5}, {4, 5, 6}};
  expected = {-1, 6, 3, 2, 5, 4, 1};
  TS_ASSERT_MATCH_WEIGHTS(expected, edges, false);
} // MeWeightMatcherUnitTests::test20_sblossom
//------------------------------------------------------------------------------
/// \brief Test create S-blossom, relabel as T-blossom, use for augmentation.
//------------------------------------------------------------------------------
void MeWeightMatcherUnitTests::test21_tblossom()
{
#ifdef DEBUG
  printf("test21_tblossom-1\n\n");
#endif
  VecInt2d edges = {{1, 2, 9}, {1, 3, 8}, {2, 3, 10}, {1, 4, 5}, {4, 5, 4}, {1, 6, 3}};
  VecInt expected = {-1, 6, 3, 2, 5, 4, 1};
  TS_ASSERT_MATCH_WEIGHTS(expected, edges, false);
#ifdef DEBUG
  printf("test21_tblossom-2\n\n");
#endif
  edges = {{1, 2, 9}, {1, 3, 8}, {2, 3, 10}, {1, 4, 5}, {4, 5, 3}, {1, 6, 4}};
  expected = {-1, 6, 3, 2, 5, 4, 1};
  TS_ASSERT_MATCH_WEIGHTS(expected, edges, false);
#ifdef DEBUG
  printf("test21_tblossom-3\n\n");
#endif
  edges = {{1, 2, 9}, {1, 3, 8}, {2, 3, 10}, {1, 4, 5}, {4, 5, 3}, {3, 6, 4}};
  expected = {-1, 2, 1, 6, 5, 4, 3};
  TS_ASSERT_MATCH_WEIGHTS(expected, edges, false);
} // MeWeightMatcherUnitTests::test21_tblossom
//------------------------------------------------------------------------------
/// \brief Test create nested S-blossom, use for augmentation.
//------------------------------------------------------------------------------
void MeWeightMatcherUnitTests::test22_s_nest()
{
#ifdef DEBUG
  printf("test22_s_nest\n\n");
#endif
  VecInt2d edges = {{1, 2, 9}, {1, 3, 9}, {2, 3, 10}, {2, 4, 8}, {3, 5, 8}, {4, 5, 10}, {5, 6, 6}};
  VecInt expected = {-1, 3, 4, 1, 2, 6, 5};
  TS_ASSERT_MATCH_WEIGHTS(expected, edges, false);
} // MeWeightMatcherUnitTests::test22_s_nest
//------------------------------------------------------------------------------
/// \brief Test create S-blossom, relabel as S, include in nested S-blossom.
//------------------------------------------------------------------------------
void MeWeightMatcherUnitTests::test23_s_relabel_nest()
{
#ifdef DEBUG
  printf("test23_s_relabel_nest\n\n");
#endif
  VecInt2d edges = {{1, 2, 10}, {1, 7, 10}, {2, 3, 12}, {3, 4, 20}, {3, 5, 20},
                    {4, 5, 25}, {5, 6, 10}, {6, 7, 10}, {7, 8, 8}};
  VecInt expected = {-1, 2, 1, 4, 3, 6, 5, 8, 7};
  TS_ASSERT_MATCH_WEIGHTS(expected, edges, false);
} // MeWeightMatcherUnitTests::test23_s_relabel_nest
//------------------------------------------------------------------------------
/// \brief Test create nested S-blossom, augment, expand recursively.
//------------------------------------------------------------------------------
void MeWeightMatcherUnitTests::test24_s_nest_expand()
{
#ifdef DEBUG
  printf("test24_s_nest_expand\n\n");
#endif
  VecInt2d edges = {{1, 2, 8},  {1, 3, 8},  {2, 3, 10}, {2, 4, 12}, {3, 5, 12},
                    {4, 5, 14}, {4, 6, 12}, {5, 7, 12}, {6, 7, 14}, {7, 8, 12}};
  VecInt expected = {-1, 2, 1, 5, 6, 3, 4, 8, 7};
  TS_ASSERT_MATCH_WEIGHTS(expected, edges, false);
} // MeWeightMatcherUnitTests::test24_s_nest_expand
//------------------------------------------------------------------------------
/// \brief Test create S-blossom, relabel as T, expand.
//------------------------------------------------------------------------------
void MeWeightMatcherUnitTests::test25_s_t_expand()
{
#ifdef DEBUG
  printf("test25_s_t_expand\n\n");
#endif
  VecInt2d edges = {{1, 2, 23}, {1, 5, 22}, {1, 6, 15}, {2, 3, 25},
                    {3, 4, 22}, {4, 5, 25}, {4, 8, 14}, {5, 7, 13}};
  VecInt expected = {-1, 6, 3, 2, 8, 7, 1, 5, 4};
  TS_ASSERT_MATCH_WEIGHTS(expected, edges, false);
} // MeWeightMatcherUnitTests::test25_s_t_expand
//------------------------------------------------------------------------------
/// \brief Test create nested S-blossom, relabel as T, expand.
//------------------------------------------------------------------------------
void MeWeightMatcherUnitTests::test26_s_nest_t_expand()
{
#ifdef DEBUG
  printf("test26_s_nest_t_expand\n\n");
#endif
  VecInt2d edges = {{1, 2, 19}, {1, 3, 20}, {1, 8, 8}, {2, 3, 25}, {2, 4, 18},
                    {3, 5, 18}, {4, 5, 13}, {4, 7, 7}, {5, 6, 7}};
  VecInt expected = {-1, 8, 3, 2, 7, 6, 5, 4, 1};
  TS_ASSERT_MATCH_WEIGHTS(expected, edges, false);
} // MeWeightMatcherUnitTests::test26_s_nest_t_expand
//------------------------------------------------------------------------------
/// \brief Test create blossom, relabel as T in more than one way, expand,
/// augment.
//------------------------------------------------------------------------------
void MeWeightMatcherUnitTests::test30_tnasty_expand()
{
#ifdef DEBUG
  printf("test30_tnasty_expand\n\n");
#endif
  VecInt2d edges = {{1, 2, 45}, {1, 5, 45}, {2, 3, 50}, {3, 4, 45}, {4, 5, 50},
                    {1, 6, 30}, {3, 9, 35}, {4, 8, 35}, {5, 7, 26}, {9, 10, 5}};
  VecInt expected = {-1, 6, 3, 2, 8, 7, 1, 5, 4, 10, 9};
  TS_ASSERT_MATCH_WEIGHTS(expected, edges, false);
} // MeWeightMatcherUnitTests::test30_tnasty_expand
//------------------------------------------------------------------------------
/// \brief Test again but slightly different.
//------------------------------------------------------------------------------
void MeWeightMatcherUnitTests::test31_tnasty2_expand()
{
#ifdef DEBUG
  printf("test31_tnasty2_expand\n\n");
#endif
  VecInt2d edges = {{1, 2, 45}, {1, 5, 45}, {2, 3, 50}, {3, 4, 45}, {4, 5, 50},
                    {1, 6, 30}, {3, 9, 35}, {4, 8, 26}, {5, 7, 40}, {9, 10, 5}};
  VecInt expected = {-1, 6, 3, 2, 8, 7, 1, 5, 4, 10, 9};
  TS_ASSERT_MATCH_WEIGHTS(expected, edges, false);
} // MeWeightMatcherUnitTests::test31_tnasty2_expand
//------------------------------------------------------------------------------
/// \brief Test create blossom, relabel as T, expand such that a new least-slack
/// S-to-free edge is produced, augment.
//------------------------------------------------------------------------------
void MeWeightMatcherUnitTests::test32_t_expand_leastslack()
{
#ifdef DEBUG
  printf("test32_t_expand_leastslack\n\n");
#endif
  VecInt2d edges = {{1, 2, 45}, {1, 5, 45}, {2, 3, 50}, {3, 4, 45}, {4, 5, 50},
                    {1, 6, 30}, {3, 9, 35}, {4, 8, 28}, {5, 7, 26}, {9, 10, 5}};
  VecInt expected = {-1, 6, 3, 2, 8, 7, 1, 5, 4, 10, 9};
  TS_ASSERT_MATCH_WEIGHTS(expected, edges, false);
} // MeWeightMatcherUnitTests::test32_t_expand_leastslack
//------------------------------------------------------------------------------
/// \brief Test create nested blossom, relabel as T in more than one way, expand
/// outer blossom such that inner blossom ends up on an augmenting path.
//------------------------------------------------------------------------------
void MeWeightMatcherUnitTests::test33_nest_tnasty_expand()
{
#ifdef DEBUG
  printf("test33_nest_tnasty_expand\n\n");
#endif
  VecInt2d edges = {{1, 2, 45}, {1, 7, 45},  {2, 3, 50}, {3, 4, 45}, {4, 5, 95},
                    {4, 6, 94}, {5, 6, 94},  {6, 7, 50}, {1, 8, 30}, {3, 11, 35},
                    {5, 9, 36}, {7, 10, 26}, {11, 12, 5}};
  VecInt expected = {-1, 8, 3, 2, 6, 9, 4, 10, 1, 5, 7, 12, 11};
  TS_ASSERT_MATCH_WEIGHTS(expected, edges, false);
} // MeWeightMatcherUnitTests::test33_nest_tnasty_expand
//------------------------------------------------------------------------------
/// \brief Test create nested S-blossom, relabel as S, expand recursively.
//------------------------------------------------------------------------------
void MeWeightMatcherUnitTests::test34_nest_relabel_expand()
{
#ifdef DEBUG
  printf("test34_nest_relabel_expand\n\n");
#endif
  VecInt2d edges = {{1, 2, 40}, {1, 3, 40}, {2, 3, 60}, {2, 4, 55},  {3, 5, 55}, {4, 5, 50},
                    {1, 8, 15}, {5, 7, 30}, {7, 6, 10}, {8, 10, 10}, {4, 9, 30}};
  VecInt expected = {-1, 2, 1, 5, 9, 3, 7, 6, 10, 4, 8};

  TS_ASSERT_MATCH_WEIGHTS(expected, edges, false);
} // MeWeightMatcherUnitTests::test34_nest_relabel_expand
//------------------------------------------------------------------------------
/// \brief Test simple triangle with three divisions on each side.
//------------------------------------------------------------------------------
void MeWeightMatcherUnitTests::testSimpleTriangle()
{
#ifdef DEBUG
  printf("testSimpleTriangle\n\n");
#endif
  // VecPt3d points = { {0, 0, 0}, {10, 0, 0}, {20, 0, 0}, {30, 0, 0},
  //                   {0, 10, 0}, {10, 10, 0}, {20, 10, 0},
  //                   {0, 20, 0}, {10, 20, 0},
  //                   {0, 30, 0}};
  // VecInt2d triangles = {{0, 1, 4}, {1, 5, 4}, {1, 2, 5}, {2, 6, 5}, {2, 3, 6},
  //                      {4, 5, 7}, {5, 8, 7}, {5, 6, 8}, {7, 8, 9}};
  // clang-format off
  VecInt2d edges = { {0, 1, 15}, {1, 2, 10}, {2, 3, 15}, {3, 4, 10},
                     {1, 5, 10}, {3, 7, 10}, {5, 6, 15}, {6, 7, 10},
                     {6, 8, 10},
                     {0, 2, -10}, {2, 4, -10}, {4, 7, -10}, {7, 8, -10}, {8, 5, -10}, {5, 0, -10}};
  // clang-format on
  VecInt expected = {1, 0, 3, 2, -1, 6, 5, -1, -1};
  //                 0, 1, 2, 3,  4, 5, 6,  7,  8

  TS_ASSERT_MATCH_WEIGHTS(expected, edges, false);
  expected = {1, 0, 3, 2, 7, 6, 5, 4, -1};
  //          0, 1, 2, 3, 4, 5, 6, 7,  8
  TS_ASSERT_MATCH_WEIGHTS(expected, edges, true);
} // MeWeightMatcherUnitTests::testSimpleTriangle
//------------------------------------------------------------------------------
/// \brief Test simple quad with three divisions on each side.
//------------------------------------------------------------------------------
void MeWeightMatcherUnitTests::testSimpleQuad()
{
#ifdef DEBUG
  printf("testSimpleQuad\n\n");
#endif
  // clang-format off
  //VecPt3d points = {
  //  {0, 0, 0}, {10, 0, 0}, {20, 0, 0}, {30, 0, 0},
  //  {0, 10, 0}, {10, 10, 0}, {20, 10, 0}, {30, 10, 0},
  //  {0, 20, 0}, {10, 20, 0}, {20, 20, 0}, {30, 20, 0},
  //  {0, 30, 0}, {10, 30, 0}, {20, 30, 0}, {30, 30, 0}};
  //VecInt2d triangles = {
  //  {0, 1, 4}, {1, 5, 4}, {1, 2, 5}, {2, 6, 5}, {2, 3, 6}, {3, 7, 6},
  //  {4, 9, 8}, {4, 5, 9}, {5, 10, 9}, {5, 6, 10}, {6, 7, 11}, {6, 11, 10},
  //  {8, 9, 12}, {9, 13, 12}, {9, 10, 13}, {10, 14, 13}, {10, 11, 14},
  //  {11, 15, 14}};
  VecInt2d edges = {
    {0, 1, 15}, {1, 2, 10}, {2, 3, 15}, {3, 4, 10}, {4, 5, 15},           // inside row 1
    {6, 7, 15}, {7, 8, 10}, {8, 9, 15}, {9, 10, 10}, {10, 11, 15},        // inside row 2
    {12, 13, 15}, {13, 14, 10}, {14, 15, 15}, {15, 16, 10}, {16, 17, 15}, // inside row 3
    {0, 2, -10}, {2, 4, -10}, {4, 5, -5}, {5, 11, -10}, {11, 17, -10},    // outside
    {17, 15, -10}, {15, 13, -10}, {12, 13, -5}, {12, 6, -10}, {6, 0, -10}
  };
  // clang-format on

  VecInt expected = {1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14, 17, 16};
  //                 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17

  TS_ASSERT_MATCH_WEIGHTS(expected, edges, false);
  TS_ASSERT_MATCH_WEIGHTS(expected, edges, true);
} // MeWeightMatcherUnitTests::testSimpleQuad
//------------------------------------------------------------------------------
/// \brief Test complex quad with three divisions on two opposite sides and
//  two on the horizontal bottom row increasing to six divisions on top row.
//------------------------------------------------------------------------------
void MeWeightMatcherUnitTests::testComplexQuad()
{
#ifdef DEBUG
  printf("testComplexQuad\n\n");
#endif
  // clang-format off
  //VecPt3d points = {
  //  {-10, 0, 0}, {0, 0, 0}, {10, 0,0},
  //  {-15, 10, 0}, {-5, 10, 0}, {5, 10, 0}, {15, 10, 0},
  //  {-20, 20, 0}, {-10, 20, 0}, {0, 20, 0}, {10, 20, 0}, {20, 20, 0},
  //  {-25, 30, 0}, {-15, 30, 0}, {-5, 30, 0}, {5, 30, 0}, {15, 30, 0}, {25, 30, 0},
  //  {-30, 40, 0}, {-20, 40, 0}, {-10, 40, 0}, {0, 40, 0}, {10, 40, 0}, {20, 40, 0}, {30, 40, 0}};
  //VecInt2d triangles = {
  //  {0, 4, 3}, {0, 1, 4}, {1, 5, 4}, {0, 2, 5}, {2, 6, 5},
  //  {3, 8, 7}, {3, 4, 8}, {4, 9, 8}, {4, 5, 9}, {5, 10, 9}, {5, 11, 10}, {5, 6, 11},
  //  {7, 13, 12}, {7, 8, 13}, {8, 14, 13}, {8, 9, 14}, {9, 15, 14}, {9, 10, 15}, {10, 16, 15}, {10, 11, 16}, {11, 17, 16},
  //  {12, 19, 18}, {12, 13, 19}, {13, 20, 19}, {13, 14, 20}, {14, 21, 20}, {14, 15, 21}, {15, 22, 21}, {15, 23, 22}, {15, 16, 23}, {16, 17, 23}, {17, 24, 23}};
  VecInt2d edges = {
    // vertical triangle edges bottom row
    {0, 1, 10}, {1, 2, 10}, {2, 3, 10}, {3, 4, 10},
    // horizontal triangle edges bottom row
    {0, 6, 10}, {2, 8, 10}, {4, 11, 10},
    // vertical triangles second row
    {5, 6, 10}, {6, 7, 10}, {7, 8, 10}, {8, 9, 10}, {9, 10, 10}, {10, 11, 10},
    // horizontal triangles second row
    {5, 13, 10}, {7, 15, 10}, {9, 17, 10}, {10, 19, 10},
    // vertical triangles third row
    {12, 13, 10}, {13, 14, 10}, {14, 15, 10}, {15, 16, 10}, {16, 17, 10},
    {17, 18, 10}, {18, 19, 10}, {19, 20, 10},
    // horizontal triangles third row
    {12, 22, 10}, {14, 24, 10}, {16, 26, 10}, {18, 29, 10}, {20, 30, 10},
    {21, 22, 10}, {22, 23, 10}, {23, 24, 10}, {24, 25, 10}, {25, 26, 10},
    {26, 27, 10}, {27, 28, 10}, {28, 29, 10}, {29, 30, 10}, {30, 31, 10},
    // boundary edges
    {0, 1, -5}, {1, 3, -10}, {3, 4, -5}, {4, 11, -5}, {11, 20, -15}, {20, 31, -10},
    {31, 28, -15}, {28, 27, -5}, {27, 25, -10}, {25, 23, -10}, {23, 21, -10}, {21, 12, -10},
    {12, 5, -10}, {5, 0, -10}
  };
  
  VecInt expected = { 1, 0, 3, 2, 11, 6, 5, 8, 7, 10,  9,  4, 13, 12, 15, 14, 17, 16, 29, 20, 19, 22, 21, 24, 23, 26, 25, 28, 27, 18, 31, 30 };
  //                  0, 1, 2, 3,  4, 5, 6, 7, 8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31
  // clang-format on

  TS_ASSERT_MATCH_WEIGHTS(expected, edges, false);
  TS_ASSERT_MATCH_WEIGHTS(expected, edges, true);
} // MeWeightMatcherUnitTests::testComplexQuad

#endif // CXX_TEST
