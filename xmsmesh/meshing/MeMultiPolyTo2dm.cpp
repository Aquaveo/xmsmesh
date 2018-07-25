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
#include <cmath>

// 4. External library headers
#include <boost/format.hpp>

// 5. Shared code headers
#include <xmscore/misc/StringUtil.h>
#include <xmscore/misc/XmConst.h>
#include <xmscore/misc/XmError.h>
#include <xmscore/misc/XmLog.h>
#include <xmscore/points/pt.h> // Pt3d
#include <xmscore/stl/vector.h>
#include <xmsinterp/geometry/geoms.h>
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

  bool Generate2dm(MeMultiPolyMesherIo& a_input, const std::string& a_outFileName,
                   int a_precision = 14) override;
  bool Generate2dm(MeMultiPolyMesherIo& a_input, std::ostream& a_os,
                   int a_precision = 14) override;

  void Write2dm(MeMultiPolyMesherIo& a_input, std::ostream& a_os, int a_precision);

  //------------------------------------------------------------------------------
  /// \brief sets the flag for testing.
  /// \param a_: The flag.
  //------------------------------------------------------------------------------
  void SetTranslatePointsForTesting(bool a_) { m_translatePointsForTesting = true; }

  //------------------------------------------------------------------------------
  /// \brief sets the observer class to give feedback on the grid generation process
  /// \param a_: The observer.
  //------------------------------------------------------------------------------
  virtual void SetObserver(BSHP<Observer> a_) override { m_prog = a_; }

  BSHP<Observer> m_prog; ///< observer class to give feedback on grid generation process

  /// translates input points to midpoint of input bounds. This removes a lot of insignificant
  /// digits when working in projected coordinate system (UTM).
  bool m_translatePointsForTesting = false; 
};

//----- Internal functions -----------------------------------------------------
//------------------------------------------------------------------------------
/// \brief finds the bounds of the input polygons and computes the middle. Then
/// the input points are translated by subtracting the middle. This can save
/// precision when working in projected coordinates (UTM).
/// \param[in,out] a_io MeMultiPolyMesherIo class that is modified by this
/// function
/// \param[out] a_ptMiddle The computed middle point of the input polygons
//------------------------------------------------------------------------------
static void iTransformInputs(MeMultiPolyMesherIo& a_io, Pt3d& a_ptMiddle)
{
    Pt3d overallMin(XM_DBL_HIGHEST), overallMax(XM_DBL_LOWEST);
    // get bounds of polygons
    for (auto& polyio : a_io.m_polys)
    {
      Pt3d polyioMin, polyioMax;
      gmExtents2D(polyio.m_outPoly, polyioMin, polyioMax);
      for (auto& poly: polyio.m_insidePolys)
      {
        Pt3d ptMin, ptMax;
        gmExtents2D(poly, ptMin, ptMax);
        gmAddToExtents(ptMin, polyioMin, polyioMax);
      }
      gmAddToExtents(polyioMin, overallMin, overallMax);
      gmAddToExtents(polyioMax, overallMin, overallMax);
    }
    // compute the middle of the polygon bounds
    a_ptMiddle = (overallMin + overallMax)/ 2.0;
    // subtract the bounds middle from all polygon points
    for (auto& polyio : a_io.m_polys)
    {
      for (auto& p : polyio.m_outPoly)
      {
        p -= a_ptMiddle;
      }
      for (auto& poly: polyio.m_insidePolys)
      {
        for (auto& p : poly)
        {
          p -= a_ptMiddle;
        }
      }
    }
    // subtract the bounds middle from all refine points
    for (auto& rpt : a_io.m_refPts)
    {
      rpt.m_pt -= a_ptMiddle;
    }
} // iTransformInputs
//------------------------------------------------------------------------------
/// \brief Translates the output points from the poly mesher based on 
/// a_ptMiddle. See iTransformInputs
/// \param[in] a_tmpIo MeMultiPolyMesherIo class that has an output mesh
/// \param[in] a_ptMiddle The computed middle point of the input polygons
/// \param[in] a_transform Flag indicating if the data needs to be transformed
/// back to the original coordinate space.
/// \param[out] a_io MeMultiPolyMesherIo class that will have the final mesh
/// transformed back to original coordinate space.
//------------------------------------------------------------------------------
static void iTranslateOutputs(const MeMultiPolyMesherIo& a_tmpIo, const Pt3d& a_ptMiddle,
                              bool a_transform, MeMultiPolyMesherIo& a_io)
{
  a_io.m_points = a_tmpIo.m_points;
  a_io.m_cells = a_tmpIo.m_cells;
  a_io.m_cellPolygons = a_tmpIo.m_cellPolygons;
  if (a_transform)
  {
    // translate the mesh outputs back to original coords
    for (auto& p : a_io.m_points)
    {
      p += a_ptMiddle;
    }
  }
} // iTranslateOutputs

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
/// \param[in] a_precision: pinted digits of precision of points in file
/// \return true if the mesh was generated.
//------------------------------------------------------------------------------
bool MeMultiPolyTo2dmImpl::Generate2dm(MeMultiPolyMesherIo& a_io, const std::string& a_outFileName,
                                       int a_precision)
{
  std::fstream os(a_outFileName.c_str(), std::fstream::out);
  if (os.bad())
    return false;
  return Generate2dm(a_io, os, a_precision);
} // MeMultiPolyTo2dmImpl::Generate2dm
//------------------------------------------------------------------------------
/// \brief Creates a 2dm file from polygons by meshing.
/// \param[in] a_io: Input/output of polygons and options for generating a mesh.
/// \param[in] a_os: output stream to store the mesh
/// \param[in] a_precision: pinted digits of precision of points in file
/// \return true if a mesh was created
//------------------------------------------------------------------------------
bool MeMultiPolyTo2dmImpl::Generate2dm(MeMultiPolyMesherIo& a_io, std::ostream& a_os,
                                       int a_precision)
{
  if (a_os.bad())
  {
    XM_LOG(xmlog::error, "Invalid output specified. Aborting");
    return false;
  }

  // translate
  Pt3d ptMiddle;
  MeMultiPolyMesherIo io = a_io;
  if (m_translatePointsForTesting)
  {
    iTransformInputs(io, ptMiddle);
  }

  // Mesh the polygons
  BSHP<MeMultiPolyMesher> mp = MeMultiPolyMesher::New();
  mp->SetObserver(m_prog);
  if (!mp->MeshIt(io))
  {
    XM_LOG(xmlog::error, "Failed to generate mesh from polygons.");
    return false;
  }

  iTranslateOutputs(io, ptMiddle, m_translatePointsForTesting, a_io);

  // write the 2dm
  Write2dm(a_io, a_os, a_precision);
  return true;
} // MeMultiPolyTo2dmImpl::GenerateUGrid
//------------------------------------------------------------------------------
/// \brief Writes 2dm data from a point and cell stream from the
/// MeMultiPolyMesherIo class.
/// \param[in] a_io: Input/output of polygons and options for generating a mesh.
/// \param[in] a_os: output stream to store the mesh
/// \param[in] a_precision: pinted digits of precision of points in file
//------------------------------------------------------------------------------
void MeMultiPolyTo2dmImpl::Write2dm(MeMultiPolyMesherIo& a_io, std::ostream& a_os,
                                    int a_precision)
{
  XM_ENSURE_TRUE_NO_ASSERT(a_io.m_points.size() && a_io.m_cells.size());

  a_os << "MESH2D\n";
  std::string tempStr("E3T"), format("%5d");
  int ival = static_cast<int>(std::log10(a_io.m_points.size())) + 1;
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

    // order the points so that the one with the lowest number is first
    int numPoints = cells[i + 1];
    int* points = &cells[i + 2]; // points for cell start at i + 2
    int lowPointNumber = points[0];
    int lowPointIndex = 0;
    for (int j = 1; j < numPoints; ++j)
    {
      if (points[j] < lowPointNumber)
      {
        lowPointNumber = points[j];
        lowPointIndex = j;
      }
    }
    for (int j = 0; j < numPoints; ++j)
    {
      int pointIndex = (j + lowPointIndex) % numPoints;
      // add 1 to change from 0 based to 1 based point numbers.
      tempStr = (boost::format(format.c_str()) % (points[pointIndex] + 1)).str();
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

    int test_precision = -1;
#ifdef TEST_PRECISION
    test_precision = TEST_PRECISION;
#endif
    a_os << "ND " << tempStr << " "
         << STRstd(p.x, test_precision, a_precision, STR_FULLWIDTH) << " "
         << STRstd(p.y, test_precision, a_precision, STR_FULLWIDTH) << " "
         << STRstd(p.z, test_precision, a_precision, STR_FULLWIDTH) << "\n";
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
/// \param a_precision: Precision used to write 2dm file points.
//------------------------------------------------------------------------------
static void iReadPolysAndCreate2dm(std::string a_filename, std::ostream& a_os,
                                   int a_precision)
{
  std::string inPolyFile(a_filename);
  MeMultiPolyTo2dmImpl p2g;
  std::vector<std::vector<std::vector<Pt3d>>> inside;
  std::vector<std::vector<Pt3d>> outside;
  MeMultiPolyMesherIo io;
  tutReadMeshIoFromFile(inPolyFile, io);
  p2g.Generate2dm(io, a_os, a_precision);
} // iReadPolysAndCreate2dm
//------------------------------------------------------------------------------
/// \brief
/// \param a_fileBase: Base of filename (prefix) of the polygon file.
/// \param a_precision: Precision used to write 2dm file points.
//------------------------------------------------------------------------------
static void iTestFromPolyFile(std::string a_fileBase, int a_precision = 14)
{
  const std::string path(std::string(XMSMESH_TEST_PATH) + "meshing/");
  // not using ttGetXmsngTestPath() because the XMSMESH_TEST_PATH is not defined
  // in xmscore.
  //const std::string path(ttGetXmsngTestPath() + "meshing/");
  std::string outFile;
  std::string baseFile;
  ttGetTestFilePaths(path, a_fileBase, ".2dm", baseFile, outFile);
  {
    std::fstream os;
    os.open(outFile.c_str(), std::fstream::out);
    iReadPolysAndCreate2dm(path + a_fileBase + ".txt", os, a_precision);
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
    "E3T     1     1    18     2     1\n"
    "E3T     2     2    24     3     1\n"
    "E3T     3     3    24    23     1\n"
    "E3T     4    17    24    18     1\n"
    "E3T     5     4    23    27     1\n"
    "E3T     6     3    23     4     1\n"
    "E3T     7     2    18    24     1\n"
    "E3T     8     5    20     6     1\n"
    "E3T     9     6    19     7     1\n"
    "E3T    10     6    20    19     1\n"
    "E3T    11    20    27    28     1\n"
    "E3T    12     7    19     8     1\n"
    "E3T    13     8    19    22     1\n"
    "E3T    14     5    27    20     1\n"
    "E3T    15     4    27     5     1\n"
    "E3T    16    23    29    27     1\n"
    "E3T    17    17    25    24     1\n"
    "E3T    18    16    25    17     1\n"
    "E3T    19    23    26    29     1\n"
    "E3T    20    15    26    25     1\n"
    "E3T    21    26    28    29     1\n"
    "E3T    22    14    26    15     1\n"
    "E3T    23    15    25    16     1\n"
    "E3T    24    13    26    14     1\n"
    "E3T    25    12    21    28     1\n"
    "E3T    26     8    22     9     1\n"
    "E3T    27     9    22    11     1\n"
    "E3T    28    12    22    21     1\n"
    "E3T    29    20    28    21     1\n"
    "E3T    30    12    28    13     1\n"
    "E3T    31     9    11    10     1\n"
    "E3T    32    11    22    12     1\n"
    "E3T    33    13    28    26     1\n"
    "E3T    34    27    29    28     1\n"
    "E3T    35    14    47    30     1\n"
    "E3T    36    14    15    47     1\n"
    "E3T    37    45    47    46     1\n"
    "E3T    38    30    47    53     1\n"
    "E3T    39    44    53    45     1\n"
    "E3T    40    30    52    31     1\n"
    "E3T    41    30    53    52     1\n"
    "E3T    42    45    53    47     1\n"
    "E3T    43    12    13    31     1\n"
    "E3T    44    10    11    33     1\n"
    "E3T    45    31    56    32     1\n"
    "E3T    46    31    52    56     1\n"
    "E3T    47    32    49    33     1\n"
    "E3T    48    33    48    34     1\n"
    "E3T    49    10    33    34     1\n"
    "E3T    50    32    56    49     1\n"
    "E3T    51    12    31    32     1\n"
    "E3T    52    52    58    56     1\n"
    "E3T    53    43    54    44     1\n"
    "E3T    54    49    56    57     1\n"
    "E3T    55    55    57    58     1\n"
    "E3T    56    52    55    58     1\n"
    "E3T    57    56    58    57     1\n"
    "E3T    58    42    55    54     1\n"
    "E3T    59    40    55    41     1\n"
    "E3T    60    41    55    42     1\n"
    "E3T    61    40    57    55     1\n"
    "E3T    62    42    54    43     1\n"
    "E3T    63    39    57    40     1\n"
    "E3T    64    49    57    50     1\n"
    "E3T    65    34    48    35     1\n"
    "E3T    66    35    48    51     1\n"
    "E3T    67    33    49    48     1\n"
    "E3T    68    38    51    39     1\n"
    "E3T    69    39    51    50     1\n"
    "E3T    70    36    38    37     1\n"
    "E3T    71    36    51    38     1\n"
    "E3T    72    35    51    36     1\n"
    "E3T    73    39    50    57     1\n"
    "E3T    74    44    54    53     1\n"
    "ND     1            0.0            0.0            0.0\n"
    "ND     2            0.0           10.0            0.0\n"
    "ND     3            0.0           20.0            0.0\n"
    "ND     4            0.0           30.0            0.0\n"
    "ND     5            0.0           40.0            0.0\n"
    "ND     6            0.0           50.0            0.0\n"
    "ND     7            0.0           60.0            0.0\n"
    "ND     8           10.0           60.0            0.0\n"
    "ND     9           20.0           60.0            0.0\n"
    "ND    10           30.0           60.0            0.0\n"
    "ND    11           30.0           50.0            0.0\n"
    "ND    12           30.0           40.0            0.0\n"
    "ND    13           30.0           30.0            0.0\n"
    "ND    14           30.0           20.0            0.0\n"
    "ND    15           30.0           10.0            0.0\n"
    "ND    16           30.0            0.0            0.0\n"
    "ND    17           20.0            0.0            0.0\n"
    "ND    18           10.0            0.0            0.0\n"
    "ND    19           10.0           50.0            0.0\n"
    "ND    20           10.0           40.0            0.0\n"
    "ND    21           20.0           40.0            0.0\n"
    "ND    22           20.0           50.0            0.0\n"
    "ND    23           10.0           20.0            0.0\n"
    "ND    24           10.0           10.0            0.0\n"
    "ND    25           20.0           10.0            0.0\n"
    "ND    26           20.0           20.0            0.0\n"
    "ND    27 8.510820254229 31.75795803291            0.0\n"
    "ND    28 19.91890189281 32.42691831149            0.0\n"
    "ND    29 14.90632313135 26.65924191911            0.0\n"
    "ND    30           40.0           20.0            0.0\n"
    "ND    31           40.0           30.0            0.0\n"
    "ND    32           40.0           40.0            0.0\n"
    "ND    33           40.0           50.0            0.0\n"
    "ND    34           40.0           60.0            0.0\n"
    "ND    35           50.0           60.0            0.0\n"
    "ND    36           60.0           60.0            0.0\n"
    "ND    37           70.0           60.0            0.0\n"
    "ND    38           70.0           50.0            0.0\n"
    "ND    39           70.0           40.0            0.0\n"
    "ND    40           70.0           30.0            0.0\n"
    "ND    41           70.0           20.0            0.0\n"
    "ND    42           70.0           10.0            0.0\n"
    "ND    43           70.0            0.0            0.0\n"
    "ND    44           60.0            0.0            0.0\n"
    "ND    45           50.0            0.0            0.0\n"
    "ND    46           40.0            0.0            0.0\n"
    "ND    47           40.0           10.0            0.0\n"
    "ND    48           50.0           50.0            0.0\n"
    "ND    49           50.0           40.0            0.0\n"
    "ND    50           60.0           40.0            0.0\n"
    "ND    51           60.0           50.0            0.0\n"
    "ND    52           50.0           20.0            0.0\n"
    "ND    53           50.0           10.0            0.0\n"
    "ND    54           60.0           10.0            0.0\n"
    "ND    55           60.0           20.0            0.0\n"
    "ND    56 48.37114000486  31.7708155091            0.0\n"
    "ND    57 59.82847955366 32.41624861662            0.0\n"
    "ND    58 54.80368621759 26.66426739907            0.0\n";
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
  iTestFromPolyFile("CasePaveSanDiego", 10);
} // MeMultiPolyTo2dmIntermediateTests::testCasePaveSanDiego
//------------------------------------------------------------------------------
/// \brief Tests two patched polys next to each other, and one paved poly next
///        to them.
//------------------------------------------------------------------------------
void MeMultiPolyTo2dmIntermediateTests::testCasePatch6()
{
  const std::string path(std::string(XMSMESH_TEST_PATH) + "meshing/");
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
  const std::string path(std::string(XMSMESH_TEST_PATH) + "meshing/");
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
  if (!tutReadPolygons(inPolyFile, outPoly, inPolys))
  {
    TS_ASSERT(std::string(fname + " does not exist.").empty());
    return;
  }
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
