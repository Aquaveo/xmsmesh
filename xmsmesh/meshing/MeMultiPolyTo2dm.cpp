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
#include <xmsmesh/meshing/MeMultiPolyTo2dm.h>

// 3. Standard library headers
#include <fstream>
#include <sstream>

// 4. External library headers
#include <boost/format.hpp>

// 5. Shared code headers
#include <xmscore/misc/StringUtil.h>
#include <xmscore/misc/XmError.h>
#include <xmscore/misc/XmLog.h>
#include <xmscore/points/pt.h> // Pt3d
#include <xmscore/stl/vector.h>
#include <xmsmesh/meshing/MeMultiPolyMesher.h>
#include <xmsmesh/meshing/MeMultiPolyMesherIo.h>
#include <xmscore/misc/carray.h>

// 6. Non-shared code headers

//----- Forward declarations ---------------------------------------------------

//----- External globals -------------------------------------------------------

//----- Namespace declaration --------------------------------------------------

namespace xms
{
//----- Constants / Enumerations -----------------------------------------------

//----- Classes / Structs ------------------------------------------------------
class MeMultiPolyTo2dmImpl : public MeMultiPolyTo2dm
{
public:
  MeMultiPolyTo2dmImpl() {}

  bool Generate2dm(MeMultiPolyMesherIo& a_input, const std::string& a_outFileName) override;
  bool Generate2dm(MeMultiPolyMesherIo& a_input, std::ostream& a_os) override;

  void Write2dm(MeMultiPolyMesherIo& a_input, std::ostream& a_os);

  //------------------------------------------------------------------------------
  /// \brief sets the observer class to give feedback on the grid generation process
  /// \param a_: The observer.
  //------------------------------------------------------------------------------
  virtual void SetObserver(BSHP<Observer> a_) override { m_prog = a_; }

  BSHP<Observer> m_prog; ///< observer class to give feedback on grid generation process
};

//----- Internal functions -----------------------------------------------------

//----- Class / Function definitions -------------------------------------------
//////////////////////////////////////////////////////////////////////////////
/// \class MeMultiPolyTo2dm
//////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief Creates a class
/// \return MeMultiPolyTo2dm.
//------------------------------------------------------------------------------
BSHP<MeMultiPolyTo2dm> MeMultiPolyTo2dm::New()
{
  BSHP<MeMultiPolyTo2dm> ret(new MeMultiPolyTo2dmImpl);
  return ret;
} // MeMultiPolyTo2dm::New

////////////////////////////////////////////////////////////////////////////////
/// \class MeMultiPolyTo2dmImpl
/// \brief Creates a VTK Unstructured Grid from polygons.
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief Creates a 2dm file from polygons
/// \param[in] a_io: Input/output of polygons and options for generating a mesh.
/// \param[in] a_outFileName: output filename
/// \return true if the mesh was generated.
//------------------------------------------------------------------------------
bool MeMultiPolyTo2dmImpl::Generate2dm(MeMultiPolyMesherIo& a_io, const std::string& a_outFileName)
{
  std::fstream os(a_outFileName.c_str(), std::fstream::out);
  if (os.bad())
    return false;
  return Generate2dm(a_io, os);
} // MeMultiPolyTo2dmImpl::Generate2dm
//------------------------------------------------------------------------------
/// \brief Creates a 2dm file from polygons by meshing.
/// \param[in] a_io: Input/output of polygons and options for generating a mesh.
/// \param[in] a_os: output stream to store the mesh
/// \return true if a mesh was created
//------------------------------------------------------------------------------
bool MeMultiPolyTo2dmImpl::Generate2dm(MeMultiPolyMesherIo& a_io, std::ostream& a_os)
{
  if (a_os.bad())
  {
    XM_LOG(xmlog::error, "Invalid output specified. Aborting");
    return false;
  }

  // Mesh the polygons
  BSHP<MeMultiPolyMesher> mp = MeMultiPolyMesher::New();
  mp->SetObserver(m_prog);
  if (!mp->MeshIt(a_io))
  {
    XM_LOG(xmlog::error, "Failed to generate mesh from polygons.");
    return false;
  }

  // write the 2dm
  Write2dm(a_io, a_os);
  return true;
} // MeMultiPolyTo2dmImpl::GenerateUGrid
//------------------------------------------------------------------------------
/// \brief Writes 2dm data from a point and cell stream from the
/// MeMultiPolyMesherIo class.
/// \param[in] a_io: Input/output of polygons and options for generating a mesh.
/// \param[in] a_os: output stream to store the mesh
//------------------------------------------------------------------------------
void MeMultiPolyTo2dmImpl::Write2dm(MeMultiPolyMesherIo& a_io, std::ostream& a_os)
{
  XM_ENSURE_TRUE_NO_ASSERT(a_io.m_points.size() && a_io.m_cells.size());

  a_os << "MESH2D\n";
  std::string tempStr("E3T"), format("%5d");
  int ival = static_cast<int>(log(a_io.m_points.size())) + 1;
  if (ival > 5)
  {
    std::stringstream ss;
    ss << "%" << ival << "d";
    format = ss.str();
  }

  int id(0);
  for (size_t i = 0; i < a_io.m_cells.size();)
  {
    ++id;
    VecInt& cells(a_io.m_cells);
    int celltype = cells[i + 0];
    tempStr = "E3T";
    if (celltype != 5) // 5 = VTK_TRIANGLE
      tempStr = "E4Q";
    a_os << tempStr;

    tempStr = (boost::format(format.c_str()) % id).str();
    a_os << " " << tempStr;

    int numPoints = cells[i + 1];
    for (int j = 0; j < numPoints; ++j)
    {
      tempStr = (boost::format(format.c_str()) % (cells[i + 2 + j] + 1)).str();
      a_os << " " << tempStr;
    }
    tempStr = (boost::format(format.c_str()) % 1).str();
    a_os << " " << tempStr << "\n";

    i = i + 2 + numPoints;
  }

  id = 0;
  for (size_t i = 0; i < a_io.m_points.size(); ++i)
  {
    ++id;
    tempStr = (boost::format(format.c_str()) % id).str();
    Pt3d& p(a_io.m_points[i]);
    a_os << "ND " << tempStr << " " << STRstd(p.x, -1, 15, STR_FULLWIDTH) << " "
         << STRstd(p.y, -1, 15, STR_FULLWIDTH) << " " << STRstd(p.z, -1, 15, STR_FULLWIDTH) << "\n";
  }
} //  MeMultiPolyTo2dmImpl::Write2dm

} // namespace xms

#if CXX_TEST
////////////////////////////////////////////////////////////////////////////////
// UNIT TESTS
////////////////////////////////////////////////////////////////////////////////
#include <xmsmesh/meshing/MeMultiPolyTo2dm.t.h>

#include <xmscore/testing/TestTools.h>
#include <xmsmesh/tutorial/TutMeshing.t.h>

//----- Namespace declaration --------------------------------------------------
using namespace xms;

namespace
{
//------------------------------------------------------------------------------
/// \brief
//------------------------------------------------------------------------------
static VecPt3d iArrayToVecPt3d(double* a_array, int a_size)
{
  VecPt3d v(a_size / 2);
  for (int i = 0; i < a_size; i += 2)
  {
    v[i / 2].x = a_array[i];
    v[i / 2].y = a_array[i + 1];
  }
  return v;
} // iArrayToVecPt3d
//------------------------------------------------------------------------------
/// \brief
/// \param a_filename: Path and filename of the polygon file.
//------------------------------------------------------------------------------
static void iReadPolysAndCreate2dm(std::string a_filename, std::ostream& a_os)
{
  std::string inPolyFile(a_filename);
  MeMultiPolyTo2dmImpl p2g;
  std::vector<std::vector<std::vector<Pt3d>>> inside;
  std::vector<std::vector<Pt3d>> outside;
  MeMultiPolyMesherIo io;
  tutReadMeshIoFromFile(inPolyFile, io);
  p2g.Generate2dm(io, a_os);
} // iReadPolysAndCreate2dm
//------------------------------------------------------------------------------
/// \brief
/// \param a_fileBase: Base of filename (prefix) of the polygon file.
//------------------------------------------------------------------------------
static void iTestFromPolyFile(std::string a_fileBase)
{
  const std::string path(ttGetXmsngTestPath() + "/meshing/");
  std::string outFile = path + a_fileBase + "_out.2dm";
  std::string baseFile = path + a_fileBase + "_base.2dm";
  {
    std::fstream os;
    os.open(outFile.c_str(), std::fstream::out);
    iReadPolysAndCreate2dm(path + a_fileBase + ".txt", os);
  }
  TS_ASSERT_TXT_FILES_EQUAL(baseFile, outFile);
  // ugExportVtkUGridAndCompare(ug, path, a_fileBase);
} // iTestFromPolyFile
//------------------------------------------------------------------------------
/// \brief
/// \verbatim
// Two adjacent outer polygons each with two inner polygons
//
// 60-   0-----1-----2----3|0----1-----2-----3-----4
//   |   |                 |                       |
// 50-  17     0-----3    4|23--22     0-----3     5
//   |   |     | inA1|     |     |     | inB1|     |
// 40-  16     1-----2    5|20--21     1-----2     6
//   |   |                 |                       |
// 30-  15      outA      6|19--18      outB       7
//   |   |                 |     |                 |
// 20-  14     0-----3    7|16--17     0-----3     8
//   |   |     | inA2|     |           | inB2|     |
// 10-  13     1-----2    8|15--14     1-----2     9
//   |   |                 |     |                 |
//  0-  12----11----10----9|    13----12----11----10
//
//       |-----|-----|-----|-----|-----|-----|-----|
//       0    10    20    30    40    50    60    70
/// \endverbatim
//------------------------------------------------------------------------------
static void iBuildTestCase4Polys(std::vector<VecPt3d>& a_outside,
                                 std::vector<std::vector<VecPt3d>>& a_inside,
                                 std::vector<double>& a_bias)
{
  a_outside.clear();
  a_inside.clear();
  a_bias.clear();

  // Set up polygon A
  {
    double outAa[] = {0,  60, 10, 60, 20, 60, 30, 60, 30, 50, 30, 40, 30, 30, 30, 20, 30, 10,
                      30, 0,  20, 0,  10, 0,  0,  0,  0,  10, 0,  20, 0,  30, 0,  40, 0,  50};
    double inA1a[] = {10, 50, 10, 40, 20, 40, 20, 50};
    double inA2a[] = {10, 20, 10, 10, 20, 10, 20, 20};
    a_outside.push_back(iArrayToVecPt3d(outAa, XM_COUNTOF(outAa)));
    std::vector<VecPt3d> insides;
    insides.push_back(iArrayToVecPt3d(inA1a, XM_COUNTOF(inA1a)));
    insides.push_back(iArrayToVecPt3d(inA2a, XM_COUNTOF(inA2a)));
    a_inside.push_back(insides);
    a_bias.push_back(1.0);
  }

  // Set up polygon B
  {
    double outBa[] = {30, 60, 40, 60, 50, 60, 60, 60, 70, 60, 70, 50, 70, 40, 70, 30,
                      70, 20, 70, 10, 70, 0,  60, 0,  50, 0,  40, 0,  40, 10, 30, 10,
                      30, 20, 40, 20, 40, 30, 30, 30, 30, 40, 40, 40, 40, 50, 30, 50};
    double inB1a[] = {50, 50, 50, 40, 60, 40, 60, 50};
    double inB2a[] = {50, 20, 50, 10, 60, 10, 60, 20};
    a_outside.push_back(iArrayToVecPt3d(outBa, XM_COUNTOF(outBa)));
    std::vector<VecPt3d> insides;
    insides.push_back(iArrayToVecPt3d(inB1a, XM_COUNTOF(inB1a)));
    insides.push_back(iArrayToVecPt3d(inB2a, XM_COUNTOF(inB2a)));
    a_inside.push_back(insides);
    a_bias.push_back(1.0);
  }
} // iBuildTestCase4Polys

} // namespace unnamed

////////////////////////////////////////////////////////////////////////////////
/// \class MeMultiPolyTo2dmUnitTests
/// \brief Tests for MeMultiPolyTo2dm.
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief    Defines the test group.
/// \return CxxTest::TestGroup reference.
//------------------------------------------------------------------------------
#ifndef CXXTEST4
//------------------------------------------------------------------------------
/// \brief
//------------------------------------------------------------------------------
const CxxTest::TestGroup& MeMultiPolyTo2dmUnitTests::group()
{
  return CxxTest::TestSuite::group();
} // MeMultiPolyToUGridIntermediateTests::group
#endif
//------------------------------------------------------------------------------
/// \brief tests creating the class
//------------------------------------------------------------------------------
void MeMultiPolyTo2dmUnitTests::testCreateClass()
{
  BSHP<MeMultiPolyTo2dm> m = MeMultiPolyTo2dm::New();
  TS_ASSERT(m);
} // MeMultiPolyToUGridIntermediateTests::testCreateClass
//------------------------------------------------------------------------------
/// \brief tests meshing 2 adajacent polygons with holes
//------------------------------------------------------------------------------
void MeMultiPolyTo2dmUnitTests::testCase4()
{
  // Create polygons
  std::vector<VecPt3d> outside;
  std::vector<std::vector<VecPt3d>> inside;
  std::vector<double> bias;
  iBuildTestCase4Polys(outside, inside, bias);
  MePolyInput poly;
  poly.m_bias = bias[0];
  poly.m_outPoly = outside[0];
  poly.m_insidePolys = inside[0];
  MeMultiPolyMesherIo io;
  io.m_polys.push_back(poly);
  poly.m_bias = bias[1];
  poly.m_outPoly = outside[1];
  poly.m_insidePolys = inside[1];
  io.m_polys.push_back(poly);

  // Mesh the polys
  std::stringstream ss;
  MeMultiPolyTo2dmImpl p2g;
  p2g.Generate2dm(io, ss);

  // Compare
  std::string base =
    "MESH2D\n"
    "E3T     1     2     1    18     1\n"
    "E3T     2     3     2    24     1\n"
    "E3T     3    23     3    24     1\n"
    "E3T     4    24    18    17     1\n"
    "E3T     5     4    23    27     1\n"
    "E3T     6     4     3    23     1\n"
    "E3T     7    24     2    18     1\n"
    "E3T     8     6     5    20     1\n"
    "E3T     9     7     6    19     1\n"
    "E3T    10    20    19     6     1\n"
    "E3T    11    28    20    27     1\n"
    "E3T    12     8     7    19     1\n"
    "E3T    13     8    19    22     1\n"
    "E3T    14     5    27    20     1\n"
    "E3T    15     5     4    27     1\n"
    "E3T    16    27    23    29     1\n"
    "E3T    17    25    24    17     1\n"
    "E3T    18    25    17    16     1\n"
    "E3T    19    23    26    29     1\n"
    "E3T    20    15    26    25     1\n"
    "E3T    21    29    26    28     1\n"
    "E3T    22    14    26    15     1\n"
    "E3T    23    15    25    16     1\n"
    "E3T    24    13    26    14     1\n"
    "E3T    25    21    28    12     1\n"
    "E3T    26     8    22     9     1\n"
    "E3T    27     9    22    11     1\n"
    "E3T    28    21    12    22     1\n"
    "E3T    29    28    21    20     1\n"
    "E3T    30    13    12    28     1\n"
    "E3T    31    11    10     9     1\n"
    "E3T    32    12    11    22     1\n"
    "E3T    33    28    26    13     1\n"
    "E3T    34    29    28    27     1\n"
    "E3T    35    14    47    30     1\n"
    "E3T    36    47    14    15     1\n"
    "E3T    37    46    45    47     1\n"
    "E3T    38    53    30    47     1\n"
    "E3T    39    53    45    44     1\n"
    "E3T    40    30    52    31     1\n"
    "E3T    41    52    30    53     1\n"
    "E3T    42    53    47    45     1\n"
    "E3T    43    12    13    31     1\n"
    "E3T    44    10    11    33     1\n"
    "E3T    45    31    56    32     1\n"
    "E3T    46    56    31    52     1\n"
    "E3T    47    49    33    32     1\n"
    "E3T    48    48    34    33     1\n"
    "E3T    49    33    34    10     1\n"
    "E3T    50    49    32    56     1\n"
    "E3T    51    31    32    12     1\n"
    "E3T    52    52    58    56     1\n"
    "E3T    53    54    44    43     1\n"
    "E3T    54    57    49    56     1\n"
    "E3T    55    55    57    58     1\n"
    "E3T    56    52    55    58     1\n"
    "E3T    57    58    57    56     1\n"
    "E3T    58    42    55    54     1\n"
    "E3T    59    55    41    40     1\n"
    "E3T    60    41    55    42     1\n"
    "E3T    61    55    40    57     1\n"
    "E3T    62    42    54    43     1\n"
    "E3T    63    57    40    39     1\n"
    "E3T    64    50    49    57     1\n"
    "E3T    65    48    35    34     1\n"
    "E3T    66    51    35    48     1\n"
    "E3T    67    49    48    33     1\n"
    "E3T    68    51    39    38     1\n"
    "E3T    69    39    51    50     1\n"
    "E3T    70    36    38    37     1\n"
    "E3T    71    38    36    51     1\n"
    "E3T    72    51    36    35     1\n"
    "E3T    73    50    57    39     1\n"
    "E3T    74    54    53    44     1\n"
    "ND     1             0.0             0.0             0.0\n"
    "ND     2             0.0            10.0             0.0\n"
    "ND     3             0.0            20.0             0.0\n"
    "ND     4             0.0            30.0             0.0\n"
    "ND     5             0.0            40.0             0.0\n"
    "ND     6             0.0            50.0             0.0\n"
    "ND     7             0.0            60.0             0.0\n"
    "ND     8            10.0            60.0             0.0\n"
    "ND     9            20.0            60.0             0.0\n"
    "ND    10            30.0            60.0             0.0\n"
    "ND    11            30.0            50.0             0.0\n"
    "ND    12            30.0            40.0             0.0\n"
    "ND    13            30.0            30.0             0.0\n"
    "ND    14            30.0            20.0             0.0\n"
    "ND    15            30.0            10.0             0.0\n"
    "ND    16            30.0             0.0             0.0\n"
    "ND    17            20.0             0.0             0.0\n"
    "ND    18            10.0             0.0             0.0\n"
    "ND    19            10.0            50.0             0.0\n"
    "ND    20            10.0            40.0             0.0\n"
    "ND    21            20.0            40.0             0.0\n"
    "ND    22            20.0            50.0             0.0\n"
    "ND    23            10.0            20.0             0.0\n"
    "ND    24            10.0            10.0             0.0\n"
    "ND    25            20.0            10.0             0.0\n"
    "ND    26            20.0            20.0             0.0\n"
    "ND    27 8.5108202542291 31.757958032906             0.0\n"
    "ND    28 19.918901892813 32.426918311488             0.0\n"
    "ND    29 14.906323131347 26.659241919112             0.0\n"
    "ND    30            40.0            20.0             0.0\n"
    "ND    31            40.0            30.0             0.0\n"
    "ND    32            40.0            40.0             0.0\n"
    "ND    33            40.0            50.0             0.0\n"
    "ND    34            40.0            60.0             0.0\n"
    "ND    35            50.0            60.0             0.0\n"
    "ND    36            60.0            60.0             0.0\n"
    "ND    37            70.0            60.0             0.0\n"
    "ND    38            70.0            50.0             0.0\n"
    "ND    39            70.0            40.0             0.0\n"
    "ND    40            70.0            30.0             0.0\n"
    "ND    41            70.0            20.0             0.0\n"
    "ND    42            70.0            10.0             0.0\n"
    "ND    43            70.0             0.0             0.0\n"
    "ND    44            60.0             0.0             0.0\n"
    "ND    45            50.0             0.0             0.0\n"
    "ND    46            40.0             0.0             0.0\n"
    "ND    47            40.0            10.0             0.0\n"
    "ND    48            50.0            50.0             0.0\n"
    "ND    49            50.0            40.0             0.0\n"
    "ND    50            60.0            40.0             0.0\n"
    "ND    51            60.0            50.0             0.0\n"
    "ND    52            50.0            20.0             0.0\n"
    "ND    53            50.0            10.0             0.0\n"
    "ND    54            60.0            10.0             0.0\n"
    "ND    55            60.0            20.0             0.0\n"
    "ND    56 48.371140004857 31.770815509095             0.0\n"
    "ND    57 59.828479553661 32.416248616621             0.0\n"
    "ND    58 54.803686217591 26.664267399073             0.0\n";
  TS_ASSERT_EQUALS(base, ss.str());
} // MeMultiPolyTo2dmUnitTests::testCase4
////////////////////////////////////////////////////////////////////////////////
/// \class MeMultiPolyTo2dmIntermediateTests
/// \brief Tests for MeMultiPolyTo2dm.
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief    Defines the test group.
/// \return CxxTest::TestGroup reference.
//------------------------------------------------------------------------------
#ifndef CXXTEST4
//------------------------------------------------------------------------------
/// \brief
//------------------------------------------------------------------------------
const CxxTest::TestGroup& MeMultiPolyTo2dmIntermediateTests::group()
{
  return *CxxTest::TestGroup::GetGroup(CxxTest::TG_INTERMEDIATE);
  // return CxxTest::TestSuite::group();
} // MeMultiPolyTo2dmIntermediateTests::group
#endif
//------------------------------------------------------------------------------
/// \brief tests meshing a square with a "c" shaped hole in it
//------------------------------------------------------------------------------
void MeMultiPolyTo2dmIntermediateTests::testCase2()
{
  iTestFromPolyFile("case2");
} // MeMultiPolyTo2dmIntermediateTests::testCase2
//------------------------------------------------------------------------------
/// \brief tests meshing a square with a "c" shaped hole in it
/// Also has refine points, size function, and elevation function
//------------------------------------------------------------------------------
void MeMultiPolyTo2dmIntermediateTests::testCase100()
{
  iTestFromPolyFile("case100");
} // MeMultiPolyTo2dmIntermediateTests::testCase100
//------------------------------------------------------------------------------
/// \brief tests meshing a square with a "c" shaped hole in it
/// Also has refine points, constant size function, and elevation function
//------------------------------------------------------------------------------
void MeMultiPolyTo2dmIntermediateTests::testCase101()
{
  iTestFromPolyFile("case101");
} // MeMultiPolyTo2dmIntermediateTests::testCase101
//------------------------------------------------------------------------------
/// \brief tests meshing a square with a "c" shaped hole in it
/// Also has refine points, constant size function, and elevation function
/// This input file has the outer polygon in reverse order.
//------------------------------------------------------------------------------
void MeMultiPolyTo2dmIntermediateTests::testCase102()
{
  iTestFromPolyFile("case102");
} // MeMultiPolyTo2dmIntermediateTests::testCase102
//------------------------------------------------------------------------------
/// \brief tests meshing a square with a "c" shaped hole in it
/// Also has refine points, constant size function, and elevation function
/// This input file has the outer polygon in reverse order and the inner.
//------------------------------------------------------------------------------
void MeMultiPolyTo2dmIntermediateTests::testCase103()
{
  iTestFromPolyFile("case103");
} // MeMultiPolyTo2dmIntermediateTests::testCase103
//------------------------------------------------------------------------------
/// \brief Test a paving bug.
//------------------------------------------------------------------------------
void MeMultiPolyTo2dmIntermediateTests::testCasePaveGeo()
{
  iTestFromPolyFile("CasePaveGeo");
} // MeMultiPolyTo2dmIntermediateTests::testCasePaveGeo
//------------------------------------------------------------------------------
/// \brief Test San Diego bay
//------------------------------------------------------------------------------
void MeMultiPolyTo2dmIntermediateTests::testCasePaveSanDiego()
{
  iTestFromPolyFile("CasePaveSanDiego");
} // MeMultiPolyTo2dmIntermediateTests::testCasePaveSanDiego
//------------------------------------------------------------------------------
/// \brief Tests two patched polys next to each other, and one paved poly next
///        to them.
//------------------------------------------------------------------------------
void MeMultiPolyTo2dmIntermediateTests::testCasePatch6()
{
  const std::string path(ttGetXmsngTestPath() + "/meshing/");
  const std::string fname(path + "CasePatch6.txt");
  std::string fbase("CasePatch6");
  std::vector<std::vector<std::vector<Pt3d>>> inside;
  std::vector<std::vector<Pt3d>> outside;
  tutReadPolygons(fname, outside, inside);
  MeMultiPolyMesherIo input;
  VecInt corner;
  input.m_polys.push_back(MePolyInput(outside[0], inside[0], 1.0, nullptr, corner));
  corner = {0, 5, 10};
  input.m_polys.push_back(MePolyInput(outside[1], inside[1], 1.0, nullptr, corner));
  corner = {5, 10, 15};
  input.m_polys.push_back(MePolyInput(outside[2], inside[2], 1.0, nullptr, corner));
  MeMultiPolyTo2dmImpl imp;

  {
    std::string outFile = path + "CasePatch6_out.2dm";
    std::fstream os(outFile.c_str(), std::fstream::out);
    imp.Generate2dm(input, os);
    os.close();
    std::string baseFile = path + "CasePatch6_base.2dm";
    TS_ASSERT_TXT_FILES_EQUAL(baseFile, outFile);
  }

  {
    std::string outFile = path + "CasePatch6a_out.2dm";
    std::fstream os(outFile.c_str(), std::fstream::out);
    input.m_polys[1].m_polyCorners[0] = 5;
    imp.Generate2dm(input, os);
    os.close();
    std::string baseFile = path + "CasePatch6a_base.2dm";
    TS_ASSERT_TXT_FILES_EQUAL(baseFile, outFile);
  }

  {
    std::string outFile = path + "CasePatch6b_out.2dm";
    std::fstream os(outFile.c_str(), std::fstream::out);
    input.m_polys[1].m_polyCorners[1] = 10;
    imp.Generate2dm(input, os);
    os.close();
    std::string baseFile = path + "CasePatch6b_base.2dm";
    TS_ASSERT_TXT_FILES_EQUAL(baseFile, outFile);
  }
} // MeMultiPolyTo2dmIntermediateTests::testCasePatch6
//------------------------------------------------------------------------------
/// \brief Test paving to a constant size with a transition factor.
//------------------------------------------------------------------------------
void MeMultiPolyTo2dmIntermediateTests::testCasePaveConstSizeTransition()
{
  const std::string path(ttGetXmsngTestPath() + "/meshing/");
  std::string fbase("CaseTransitionToConstSize");
  std::string fname(path + "CaseTransitionToConstSize.txt");
  std::string inPolyFile(fname);
  std::vector<std::vector<std::vector<Pt3d>>> inside;
  std::vector<std::vector<Pt3d>> outside;
  MeMultiPolyMesherIo ip;
  ip.m_polys.push_back(MePolyInput());
  MePolyInput& pp(ip.m_polys.back());
  VecPt3d2d outPoly;
  VecPt3d3d inPolys;
  tutReadPolygons(inPolyFile, outPoly, inPolys);
  pp.m_outPoly = outPoly[0];
  pp.m_insidePolys = inPolys[0];
  pp.m_constSizeBias = 0.15;
  pp.m_constSizeFunction = 100.0;

  BSHP<MeMultiPolyTo2dm> p2g = MeMultiPolyTo2dm::New();
  std::string outFile, baseFile;
  outFile = path + "CaseTransitionToConstSize_out.2dm";
  baseFile = path + "CaseTransitionToConstSize_base.2dm";
  p2g->Generate2dm(ip, outFile);
  TS_ASSERT_TXT_FILES_EQUAL(baseFile, outFile);
} // MeMultiPolyTo2dmIntermediateTests::testCasePaveConstSizeTransition

#endif // CXX_TEST
