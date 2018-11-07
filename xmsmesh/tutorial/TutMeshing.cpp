//------------------------------------------------------------------------------
/// \file
/// \ingroup tutorial
/// \copyright (C) Copyright Aquaveo 2018. Distributed under the xmsng
///  Software License, Version 1. (See accompanying file
///  LICENSE_1_0.txt or copy at http://www.aquaveo.com/xmsng/LICENSE_1_0.txt)
//------------------------------------------------------------------------------

//----- Included files ---------------------------------------------------------

// 1. Precompiled header

// 2. My own header

// 3. Standard library headers

// 4. External library headers

// 5. Shared code headers

// 6. Non-shared code headers

//----- Forward declarations ---------------------------------------------------

//----- External globals -------------------------------------------------------

//----- Namespace declaration --------------------------------------------------

//----- Constants / Enumerations -----------------------------------------------

//----- Classes / Structs ------------------------------------------------------

//----- Internal functions -----------------------------------------------------

//----- Class / Function definitions -------------------------------------------

#if CXX_TEST
////////////////////////////////////////////////////////////////////////////////
// UNIT TESTS
////////////////////////////////////////////////////////////////////////////////
#include <xmsmesh/tutorial/TutMeshing.t.h>

#include <fstream>

#include <xmscore/testing/TestTools.h>
#include <xmsinterp/geometry/geoms.h>
#include <xmsinterp/interpolate/InterpIdw.h>
#include <xmsinterp/interpolate/InterpLinear.h>
#include <xmsinterp/triangulate/TrTin.h>
#include <xmsinterp/triangulate/TrTriangulatorPoints.h>
#include <xmsmesh/meshing/MeMeshUtils.h>
#include <xmsmesh/meshing/MeMultiPolyMesherIo.h>
#include <xmsmesh/meshing/MeMultiPolyTo2dm.h>
#include <xmsmesh/meshing/MePolyRedistributePts.h>

namespace xms
{
//------------------------------------------------------------------------------
/// \brief  helper function to read MeMultiPolyMesherIo
/// \param a_fname File name where the polygon is stored
/// \param a_io MeMultiPolyMesherIo class that is filled by this function
/// \return true if the file was read with no errors
//------------------------------------------------------------------------------
bool tutReadMeshIoFromFile(const std::string& a_fname, MeMultiPolyMesherIo& a_io)
{
  a_io = MeMultiPolyMesherIo();
  std::fstream os(a_fname, std::fstream::in);
  if (!os.is_open())
    return false;

  std::string card;
  size_t numpts;
  Pt3d pt;
  MePolyInput* p(nullptr);
  while (os.good())
  {
    os >> card;
    if ("END_POLYGON" == card)
    { // do nothing
      p = nullptr;
    }
    else if ("BEGIN_POLYGON" == card)
    {
      a_io.m_polys.push_back(MePolyInput());
      p = &a_io.m_polys.back();
    }
    else if ("OUTSIDE" == card && p)
    {
      os >> numpts;
      p->m_outPoly.reserve(numpts);
      for (size_t i = 0; i < numpts; ++i)
      {
        os >> pt.x >> pt.y;
        p->m_outPoly.push_back(pt);
      }

      VecPt3d& vPoly(p->m_outPoly);
      double area = gmPolygonArea(&vPoly[0], vPoly.size());
      if (area > 0.0)
      {
        std::reverse(vPoly.begin(), vPoly.end());
        area = gmPolygonArea(&vPoly[0], vPoly.size());
      }
    }
    else if ("INSIDE" == card && p)
    {
      p->m_insidePolys.push_back(VecPt3d());
      std::vector<Pt3d>& in(p->m_insidePolys.back());
      os >> numpts;
      in.reserve(numpts);
      for (size_t i = 0; i < numpts; ++i)
      {
        os >> pt.x >> pt.y;
        in.push_back(pt);
      }

      VecPt3d& vPoly(in);
      double area = gmPolygonArea(&vPoly[0], vPoly.size());
      if (area < 0.0)
      {
        std::reverse(vPoly.begin(), vPoly.end());
        area = gmPolygonArea(&vPoly[0], vPoly.size());
      }
    }
    else if ("BIAS" == card && p)
    {
      os >> p->m_bias;
    }
    else if (("SIZE_FUNCTION" == card || "ELEVATION_FUNCTION" == card) && p)
    {
      BSHP<InterpBase> interp;
      std::string interpType;
      os >> interpType; // LINEAR or IDW
      if ("LINEAR" == interpType)
      {
        interp = InterpLinear::New();
      }
      else if ("IDW" == interpType)
      {
        interp = InterpIdw::New();
      }

      os >> numpts;
      BSHP<VecPt3d> vpts(new VecPt3d());
      VecPt3d& pts(*vpts);
      for (int i = 0; i < (int)numpts; ++i)
      {
        os >> pt.x >> pt.y >> pt.z;
        pts.push_back(pt);
      }

      if (interp)
      {
        BSHP<VecInt> tris(new VecInt());
        TrTriangulatorPoints tri(pts, *tris);
        tri.Triangulate();
        interp->SetPtsTris(vpts, tris);
        if ("SIZE_FUNCTION" == card)
          p->m_sizeFunction = interp;
        else
          p->m_elevFunction = interp;
      }
    }
    else if ("CONST_SIZE_FUNCTION" == card && p)
    {
      os >> p->m_constSizeFunction >> p->m_constSizeBias;
    }
    else if ("PATCH_CORNERS" == card && p)
    {
      p->m_polyCorners.assign(3, -1);
      os >> p->m_polyCorners[0] >> p->m_polyCorners[1] >> p->m_polyCorners[2];
    }
    else if ("CHECK_TOPOLOGY" == card)
    {
      a_io.m_checkTopology = true;
    }
    else if ("RETURN_CELL_POLYGONS" == card)
    {
      a_io.m_returnCellPolygons = true;
    }
    else if ("REFINE_POINTS" == card)
    {
      os >> numpts;
      a_io.m_refPts.reserve(numpts);
      MeRefinePoint rpt(Pt3d(), -1, false);
      for (int i = 0; i < (int)numpts; ++i)
      {
        os >> rpt.m_pt.x >> rpt.m_pt.y >> rpt.m_size >> rpt.m_createMeshPoint;
        a_io.m_refPts.push_back(rpt);
      }
    }
    else if ("SEED_POINTS" == card)
    {
      os >> numpts;
      VecPt3d& pts = p->m_seedPoints;
      pts.assign(numpts, Pt3d());
      for (int i = 0; i < numpts; ++i)
      {
        os >> pts[i].x >> pts[i].y;
      }
    }
    else if ("RELAXATION_METHOD" == card)
    {
      std::string method;
      os >> method;
      p->m_relaxationMethod = method;
    }
    card = "";
  }
  return true;
} // tutReadMeshIoFromFile
//------------------------------------------------------------------------------
/// \brief  helper function to read polygons from a text file
/// \param a_fname File name where the polygon is stored
/// \param a_outside Vector of locations for the definition of the outside of
/// the polygon.
/// \param a_inside 2d vector of locations to define multiple inside polygons
/// \return true if the file was read without any errors
//------------------------------------------------------------------------------
bool tutReadPolygons(const std::string& a_fname, VecPt3d2d& a_outside, VecPt3d3d& a_inside)
{
  a_outside.resize(0);
  a_inside.resize(0);

  MeMultiPolyMesherIo io;
  if (!tutReadMeshIoFromFile(a_fname, io))
    return false;
  for (size_t i = 0; i < io.m_polys.size(); ++i)
  {
    a_outside.push_back(io.m_polys[i].m_outPoly);
    a_inside.push_back(io.m_polys[i].m_insidePolys);
  }
  return true;
} // tutReadPolygons
//------------------------------------------------------------------------------
/// \brief  helper function to generate 2dm file and compare to a baseline
/// \param a_io mesher input
/// \param a_fileBase base file name for output file and base line file
//------------------------------------------------------------------------------
void tutGenerateAndCompare2dm(MeMultiPolyMesherIo& a_io, const std::string& a_fileBase)
{
  const std::string path(std::string(XMS_TEST_PATH) + "Tutorial_Meshing/");
  std::string outFile;
  std::string baseFile;
  ttGetTestFilePaths(path, a_fileBase, ".2dm", baseFile, outFile);
  {
    std::fstream os(outFile.c_str(), std::fstream::out);
    TS_ASSERT(!os.bad());
    if (os.bad())
      return;
    BSHP<MeMultiPolyTo2dm> pm = MeMultiPolyTo2dm::New();
    pm->Generate2dm(a_io, os, 10);
  }
  TS_ASSERT_TXT_FILES_EQUAL(baseFile, outFile);
} // tutGenerateAndCompare2dm

} // namespace xms

////////////////////////////////////////////////////////////////////////////////
/// \class TutMeshingIntermediateTests
/// \brief Tests for meshing tutorial.
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
/// \brief    Defines the test group.
/// \return CxxTest::TestGroup reference.
//------------------------------------------------------------------------------
#ifndef CXXTEST4
//------------------------------------------------------------------------------
/// \brief
//------------------------------------------------------------------------------
const CxxTest::TestGroup& TutMeshingIntermediateTests::group()
{
  return *CxxTest::TestGroup::GetGroup(CxxTest::TG_INTERMEDIATE);
  // return CxxTest::TestSuite::group();
} // TutMeshingIntermediateTests::group
#endif
//------------------------------------------------------------------------------
/// \brief Example for meshing a simple square polygon
//------------------------------------------------------------------------------
//! [snip_test_Example_SimplePolygon]
void TutMeshingIntermediateTests::test_Example_SimplePolygon()
{
  // The MePolyInput class defines a single polygon to be meshed and which
  // consists of one outer polygon and, optionally, inner polygons and some
  // other options.
  xms::MePolyInput inputPoly;

  // Outer polygon points are in clockwise order. Do not repeat the first point.
  inputPoly.m_outPoly = {
    {0, 10},   {0, 20},   {0, 30},   {0, 40},    {0, 50},   {0, 60},   {0, 70},   {0, 80},
    {0, 90},   {0, 100},  {10, 100}, {20, 100},  {30, 100}, {40, 100}, {50, 100}, {60, 100},
    {70, 100}, {80, 100}, {90, 100}, {100, 100}, {100, 90}, {100, 80}, {100, 70}, {100, 60},
    {100, 50}, {100, 40}, {100, 30}, {100, 20},  {100, 10}, {100, 0},  {90, 0},   {80, 0},
    {70, 0},   {60, 0},   {50, 0},   {40, 0},    {30, 0},   {20, 0},   {10, 0},   {0, 0}};

  // The MeMultiPolyMesherIo class holds all the options for generating
  // UGrids from polygon data
  xms::MeMultiPolyMesherIo input;

  // The m_polys vector holds the polygons that will be filled with a UGrid.
  // This case has only 1 polygon.
  input.m_polys.push_back(inputPoly);

  // generate the mesh and check the base line
  const std::string baseFile = "Example_SimplePolygon";
  tutGenerateAndCompare2dm(input, baseFile);
} // TutMeshingIntermediateTests::test_Example_SimplePolygon
//! [snip_test_Example_SimplePolygon]
//------------------------------------------------------------------------------
/// \brief Example for meshing a complex polygon
//------------------------------------------------------------------------------
//! [snip_test_Example_ComplexPolygon]
void TutMeshingIntermediateTests::test_Example_ComplexPolygon()
{
  // The MePolyInput class defines a single polygon to be meshed and which
  // consists of one outer polygon and, optionally, inner polygons and some
  // other options.
  xms::MePolyInput inputPoly;

  // Outer polygon points are in clockwise order. Do not repeat the first point.
  inputPoly.m_outPoly = {
    {0, 10},    {0, 20},    {0, 30},    {0, 40},    {0, 50},    {0, 60},    {0, 70},    {0, 80},
    {0, 90},    {0, 100},   {0, 110},   {0, 120},   {0, 130},   {0, 140},   {0, 150},   {0, 160},
    {0, 170},   {0, 180},   {0, 190},   {0, 200},   {10, 200},  {20, 200},  {30, 200},  {40, 200},
    {50, 200},  {60, 200},  {70, 200},  {80, 200},  {90, 200},  {100, 200}, {110, 200}, {120, 200},
    {130, 200}, {140, 200}, {150, 200}, {160, 200}, {170, 200}, {180, 200}, {190, 200}, {200, 200},
    {200, 190}, {200, 180}, {200, 170}, {200, 160}, {200, 150}, {200, 140}, {200, 130}, {200, 120},
    {200, 110}, {200, 100}, {200, 90},  {200, 80},  {200, 70},  {200, 60},  {200, 50},  {200, 40},
    {200, 30},  {200, 20},  {200, 10},  {200, 0},   {190, 0},   {180, 0},   {170, 0},   {160, 0},
    {150, 0},   {140, 0},   {130, 0},   {120, 0},   {110, 0},   {110, 10},  {110, 20},  {110, 30},
    {110, 40},  {120, 40},  {130, 40},  {140, 40},  {150, 40},  {150, 50},  {150, 60},  {150, 70},
    {150, 80},  {150, 90},  {150, 100}, {150, 110}, {150, 120}, {150, 130}, {150, 140}, {150, 150},
    {140, 150}, {130, 150}, {120, 150}, {110, 150}, {100, 150}, {90, 150},  {80, 150},  {70, 150},
    {60, 150},  {50, 150},  {50, 140},  {50, 130},  {50, 120},  {50, 110},  {50, 100},  {50, 90},
    {50, 80},   {50, 70},   {50, 60},   {50, 50},   {50, 40},   {60, 40},   {70, 40},   {80, 40},
    {90, 40},   {90, 30},   {90, 20},   {90, 10},   {90, 0},    {80, 0},    {70, 0},    {60, 0},
    {50, 0},    {40, 0},    {30, 0},    {20, 0},    {10, 0},    {0, 0}};

  // The MeMultiPolyMesherIo class holds all the options for generating
  // UGrids from polygon data
  xms::MeMultiPolyMesherIo input;

  // The m_polys vector holds the polygons that will be filled with a UGrid.
  // This case has only 1 polygon.
  input.m_polys.push_back(inputPoly);

  // generate the mesh and check the base line
  const std::string baseFile = "Example_ComplexPolygon";
  tutGenerateAndCompare2dm(input, baseFile);
} // TutMeshingIntermediateTests::test_Example_ComplexPolygon
//! [snip_test_Example_ComplexPolygon]
//------------------------------------------------------------------------------
/// \brief Example for meshing a simple polygon with a hole
//------------------------------------------------------------------------------
//! [snip_test_Example_SimplePolygonWithHole]
void TutMeshingIntermediateTests::test_Example_SimplePolygonWithHole()
{
  // The MePolyInput class defines a single polygon to be meshed and which
  // consists of one outer polygon and, optionally, inner polygons and some
  // other options.
  xms::MePolyInput inputPoly;

  // Outer polygon points are in clockwise order. Do not repeat the first point.
  inputPoly.m_outPoly = {
    {0, 10},   {0, 20},   {0, 30},   {0, 40},    {0, 50},   {0, 60},   {0, 70},   {0, 80},
    {0, 90},   {0, 100},  {10, 100}, {20, 100},  {30, 100}, {40, 100}, {50, 100}, {60, 100},
    {70, 100}, {80, 100}, {90, 100}, {100, 100}, {100, 90}, {100, 80}, {100, 70}, {100, 60},
    {100, 50}, {100, 40}, {100, 30}, {100, 20},  {100, 10}, {100, 0},  {90, 0},   {80, 0},
    {70, 0},   {60, 0},   {50, 0},   {40, 0},    {30, 0},   {20, 0},   {10, 0},   {0, 0}};

  // Inner polygons are in counter clockwise order. Do not repeat the first
  // point.
  inputPoly.m_insidePolys.push_back(xms::VecPt3d());
  inputPoly.m_insidePolys[0] = {{40, 40}, {50, 40}, {60, 40}, {60, 50},
                                {60, 60}, {50, 60}, {40, 60}, {40, 50}};

  // The MeMultiPolyMesherIo class holds all the options for generating
  // UGrids from polygon data
  xms::MeMultiPolyMesherIo input;

  // The m_polys vector holds the polygons that will be filled with a UGrid.
  // This case has only 1 polygon.
  input.m_polys.push_back(inputPoly);

  // generate the mesh and check the base line
  const std::string baseFile = "Example_SimplePolygonWithHole";
  tutGenerateAndCompare2dm(input, baseFile);
} // TutMeshingIntermediateTests::test_Example_SimplePolygonWithHole
//! [snip_test_Example_SimplePolygonWithHole]
//------------------------------------------------------------------------------
/// \brief Example for a simple polygon with a breakline
//------------------------------------------------------------------------------
//! [snip_test_Example_Breakline]
void TutMeshingIntermediateTests::test_Example_Breakline()
{
  // The MePolyInput class defines a single polygon to be meshed and which
  // consists of one outer polygon and, optionally, inner polygons and some
  // other options.
  xms::MePolyInput inputPoly;

  // Outer polygon points are in clockwise order. Do not repeat the first point.
  inputPoly.m_outPoly = {
    {0, 10},   {0, 20},   {0, 30},   {0, 40},    {0, 50},   {0, 60},   {0, 70},   {0, 80},
    {0, 90},   {0, 100},  {10, 100}, {20, 100},  {30, 100}, {40, 100}, {50, 100}, {60, 100},
    {70, 100}, {80, 100}, {90, 100}, {100, 100}, {100, 90}, {100, 80}, {100, 70}, {100, 60},
    {100, 50}, {100, 40}, {100, 30}, {100, 20},  {100, 10}, {100, 0},  {90, 0},   {80, 0},
    {70, 0},   {60, 0},   {50, 0},   {40, 0},    {30, 0},   {20, 0},   {10, 0},   {0, 0}};

  // Inner polygons are in counter clockwise order. Do not repeat the first
  // point.
  inputPoly.m_insidePolys.push_back(xms::VecPt3d());
  inputPoly.m_insidePolys[0] = {{50, 0}, {50, 10}, {50, 20}, {50, 10}};

  // The MeMultiPolyMesherIo class holds all the options for generating
  // UGrids from polygon data
  xms::MeMultiPolyMesherIo input;

  // The m_polys vector holds the polygons that will be filled with a UGrid.
  // This case has only 1 polygon.
  input.m_polys.push_back(inputPoly);

  // generate the mesh and check the base line
  const std::string baseFile = "Example_Breakline";
  tutGenerateAndCompare2dm(input, baseFile);
} // TutMeshingIntermediateTests::test_Example_Breakline
//! [snip_test_Example_Breakline]
//------------------------------------------------------------------------------
/// \brief Example for a simple polygon with refine points
//------------------------------------------------------------------------------
//! [snip_test_Example_RefinePoints]
void TutMeshingIntermediateTests::test_Example_RefinePoints()
{
  // The MePolyInput class defines a single polygon to be meshed and which
  // consists of one outer polygon and, optionally, inner polygons and some
  // other options.
  xms::MePolyInput inputPoly;

  // Outer polygon points are in clockwise order. Do not repeat the first point.
  inputPoly.m_outPoly = {
    {0, 10},   {0, 20},   {0, 30},   {0, 40},    {0, 50},   {0, 60},   {0, 70},   {0, 80},
    {0, 90},   {0, 100},  {10, 100}, {20, 100},  {30, 100}, {40, 100}, {50, 100}, {60, 100},
    {70, 100}, {80, 100}, {90, 100}, {100, 100}, {100, 90}, {100, 80}, {100, 70}, {100, 60},
    {100, 50}, {100, 40}, {100, 30}, {100, 20},  {100, 10}, {100, 0},  {90, 0},   {80, 0},
    {70, 0},   {60, 0},   {50, 0},   {40, 0},    {30, 0},   {20, 0},   {10, 0},   {0, 0}};

  // The MeMultiPolyMesherIo class holds all the options for generating
  // UGrids from polygon data
  xms::MeMultiPolyMesherIo input;

  // Refine points are specified independent of polygons. The mesher will
  // determine the polygon that contains the point.
  input.m_refPts.assign(3, xms::MeRefinePoint(xms::Pt3d(), -1, true));

  // specify a refine point where the point is at the cell center with
  // a desired size of 2
  input.m_refPts[0].m_pt = xms::Pt3d(20, 20, 0);
  input.m_refPts[0].m_createMeshPoint = false;
  input.m_refPts[0].m_size = 2.0;

  // specify a refine point where the point is at a mesh node with a
  // desired size of 7
  input.m_refPts[1].m_pt = xms::Pt3d(20, 80, 0);
  input.m_refPts[1].m_createMeshPoint = true;
  input.m_refPts[1].m_size = 5.0;

  // specify a "hard point"
  input.m_refPts[2].m_pt = xms::Pt3d(80, 20, 0);
  input.m_refPts[2].m_size = -1;

  // The m_polys vector holds the polygons that will be filled with a UGrid.
  // This case has only 1 polygon.
  input.m_polys.push_back(inputPoly);

  // generate the mesh and check the base line
  const std::string baseFile = "Example_RefinePoints";
  tutGenerateAndCompare2dm(input, baseFile);
} // TutMeshingIntermediateTests::test_Example_RefinePoints
//! [snip_test_Example_RefinePoints]
//------------------------------------------------------------------------------
/// \brief Example for multiple polygons with variable spacing, holes,
/// breaklines, and refine points.
//------------------------------------------------------------------------------
//! [snip_test_Example_MultiPolygon]
void TutMeshingIntermediateTests::test_Example_MultiplePolygons()
{
  const std::string path(std::string(XMS_TEST_PATH) + "Tutorial_Meshing/");
  const std::string fname(path + "Example_MultiPolys.txt");
  xms::VecPt3d3d inside;
  xms::VecPt3d2d outside;
  // Read the polygons from a text file to avoid typing out all of the
  // coordinates
  xms::tutReadPolygons(fname, outside, inside);

  // put the polygons into the meshing input class
  xms::MeMultiPolyMesherIo input;
  for (size_t i = 0; i < outside.size(); ++i)
  {
    input.m_polys.push_back(xms::MePolyInput(outside[i], inside[i]));
  }
  // add refine points to the meshing input
  xms::Pt3d p0(80, 80, 0), p1(125, 125, 0);
  input.m_refPts.push_back(xms::MeRefinePoint(p0, -1.0, true));
  input.m_refPts.push_back(xms::MeRefinePoint(p1, 1.0, true));

  // generate the mesh and check the base line
  const std::string baseFile = "Example_MultiPolys";
  tutGenerateAndCompare2dm(input, baseFile);
} // TutMeshingIntermediateTests::test_Example_MultiplePolygons
//! [snip_test_Example_MultiPolygon]
//------------------------------------------------------------------------------
/// \brief Example for scalar paving.
//------------------------------------------------------------------------------
//! [snip_test_Example_ScalarPaving]
void TutMeshingIntermediateTests::test_Example_ScalarPaving()
{
  // The MePolyInput class defines a single polygon to be meshed and which
  // consists of one outer polygon and, optionally, inner polygons and some
  // other options.
  xms::MePolyInput inputPoly;

  // Outer polygon points are in clockwise order. Do not repeat the first point.
  inputPoly.m_outPoly = {
    {0, 10},   {0, 20},   {0, 30},   {0, 40},    {0, 50},   {0, 60},   {0, 70},   {0, 80},
    {0, 90},   {0, 100},  {10, 100}, {20, 100},  {30, 100}, {40, 100}, {50, 100}, {60, 100},
    {70, 100}, {80, 100}, {90, 100}, {100, 100}, {100, 90}, {100, 80}, {100, 70}, {100, 60},
    {100, 50}, {100, 40}, {100, 30}, {100, 20},  {100, 10}, {100, 0},  {90, 0},   {80, 0},
    {70, 0},   {60, 0},   {50, 0},   {40, 0},    {30, 0},   {20, 0},   {10, 0},   {0, 0}};

  // The MeMultiPolyMesherIo class holds all the options for generating
  // UGrids from polygon data
  xms::MeMultiPolyMesherIo input;

  // The m_polys vector holds the polygons that will be filled with a UGrid.
  // This case has only 1 polygon.
  input.m_polys.push_back(inputPoly);

  // create a size function interpolator

  // These are the interpolator point locations. The size is specified by the
  // "z" component of the points.
  BSHP<xms::VecPt3d> sPts(new xms::VecPt3d());
  *sPts = {{-10, -10, 10}, {-10, 110, 10}, {110, 110, 10}, {110, -10, 10}, {60, 70, 1}};
  // These are the interpolator triangles.
  // Triangles are specified as point indexes in ccw order.  You can see the
  // 4 triangles in the vector below.
  BSHP<xms::VecInt> sTris(new xms::VecInt());
  *sTris = {0, 4, 1, 1, 4, 2, 2, 4, 3, 3, 4, 0};
  // create a linear interpolator for the size function
  BSHP<xms::InterpBase> linear(xms::InterpLinear::New());
  // sets the points and the triangles for the interpolator
  linear->SetPtsTris(sPts, sTris);
  // now set the size function on the mesh generator input class
  input.m_polys[0].m_sizeFunction = linear;

  // generate the mesh and check the base line
  const std::string baseFile = "Example_ScalarPaving";
  tutGenerateAndCompare2dm(input, baseFile);
} // TutMeshingIntermediateTests::test_Example_ScalarPaving
//! [snip_test_Example_ScalarPaving]
//------------------------------------------------------------------------------
/// \brief Example for generating a UGrid with the patch algorithm.
//------------------------------------------------------------------------------
//! [snip_test_Example_Patch]
void TutMeshingIntermediateTests::test_Example_Patch()
{
  // The MePolyInput class defines a single polygon to be meshed and which
  // consists of one outer polygon and, optionally, inner polygons and some
  // other options.
  xms::MePolyInput inputPoly;

  // Outer polygon points are in clockwise order. Do not repeat the first point.
  inputPoly.m_outPoly = {{0, 0},    {0, 10},    {0, 20},   {0, 30},   {0, 40},   {0, 60},
                         {0, 70},   {0, 80},    {0, 90},   {0, 100},  {20, 100}, {40, 100},
                         {60, 100}, {100, 100}, {100, 90}, {100, 80}, {100, 70}, {100, 60},
                         {100, 50}, {100, 40},  {100, 30}, {100, 20}, {100, 10}, {100, 0},
                         {85, 0},   {70, 0},    {55, 0},   {40, 0},   {25, 0},   {10, 0}};

  // Specify the polygon corners. It is assumed that the first point above is
  // one of the corners so we specify the other 3 corners.
  inputPoly.m_polyCorners = {9, 13, 23};

  // The MeMultiPolyMesherIo class holds all the options for generating
  // UGrids from polygon data
  xms::MeMultiPolyMesherIo input;

  // The m_polys vector holds the polygons that will be filled with a UGrid.
  // This case has only 1 polygon.
  input.m_polys.push_back(inputPoly);

  // generate the mesh and check the base line
  const std::string baseFile = "Example_Patch";
  tutGenerateAndCompare2dm(input, baseFile);
} // TutMeshingIntermediateTests::test_Example_Patch
//! [snip_test_Example_Patch]
//------------------------------------------------------------------------------
/// \brief Example for smooth transition with a constant size function
//------------------------------------------------------------------------------
//! [snip_test_Example_ConstantSmooth]
void TutMeshingIntermediateTests::test_Example_ConstantSmooth()
{
  // The MePolyInput class defines a single polygon to be meshed and which
  // consists of one outer polygon and, optionally, inner polygons and some
  // other options.
  xms::MePolyInput inputPoly;

  // Outer polygon points are in clockwise order. Do not repeat the first point.
  inputPoly.m_outPoly = {
    {0, 10},   {0, 20},   {0, 30},   {0, 40},    {0, 50},   {0, 60},   {0, 70},   {0, 80},
    {0, 90},   {0, 100},  {10, 100}, {20, 100},  {30, 100}, {40, 100}, {50, 100}, {60, 100},
    {70, 100}, {80, 100}, {90, 100}, {100, 100}, {100, 90}, {100, 80}, {100, 70}, {100, 60},
    {100, 50}, {100, 40}, {100, 30}, {100, 20},  {100, 10}, {100, 0},  {90, 0},   {80, 0},
    {70, 0},   {60, 0},   {50, 0},   {40, 0},    {30, 0},   {20, 0},   {10, 0},   {0, 0}};

  // specify a constant size function
  inputPoly.m_constSizeFunction = 1.0;
  inputPoly.m_constSizeBias = 0.2;

  // The MeMultiPolyMesherIo class holds all the options for generating
  // UGrids from polygon data
  xms::MeMultiPolyMesherIo input;

  // The m_polys vector holds the polygons that will be filled with a UGrid.
  // This case has only 1 polygon.
  input.m_polys.push_back(inputPoly);

  // generate the mesh and check the base line
  const std::string baseFile = "Example_ConstantSmooth";
  tutGenerateAndCompare2dm(input, baseFile);

  // now make the polygon grow to a bigger size element
  input.m_polys.front().m_constSizeFunction = 50;

  const std::string baseFile2 = "Example_ConstantSmooth2";
  tutGenerateAndCompare2dm(input, baseFile2);
} // TutMeshingIntermediateTests::test_Example_ConstantSmooth
//! [snip_test_Example_ConstantSmooth]
//------------------------------------------------------------------------------
/// \brief Example for creating a size function from a set of points with
/// depth scalars
//------------------------------------------------------------------------------
//! [snip_test_Example_SizeFuncFromDepth]
void TutMeshingIntermediateTests::test_Example_SizeFuncFromDepth()
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

  // now create an interpolator to pass along to a mesher

  // init the locations of the points used to interpolate
  BSHP<xms::VecPt3d> pts(new xms::VecPt3d());
  *pts = {{10, 10}, {25, 10}, {10, 25}, {50, 10}, {50, 25}, {50, 50}, {25, 50}};
  BSHP<xms::VecInt> tris(new xms::VecInt());

  // convert the sizes to a float array for the interpolator
  xms::VecFlt sizeFlt(elemSize.size(), 0);
  int i(0);
  for (auto& d : elemSize)
    sizeFlt[i++] = (float)d;

  // create an IDW interpolator
  BSHP<xms::InterpBase> interp = xms::InterpIdw::New();
  interp->SetPtsTris(pts, tris);

  // now the interpolator could be used with a mesher
  xms::MePolyInput poly;
  poly.m_sizeFunction = interp;

} // TutMeshingIntermediateTests::test_Example_SizeFuncFromDepth
//! [snip_test_Example_SizeFuncFromDepth]
//------------------------------------------------------------------------------
/// \brief Example for smoothing a size function
/// \verbatim
///   *     - point location
///   []    - array index of point
/// 1 / 100 - size function value
///
/// Input
///
///       100     100     100      100
///  20- * [8]   * [9]    * [10]   * [11]
///    |
///    |
///    |  1       100      100      100
///  10- * [4]    * [5]    * [6]    * [7]
///    |
///    |
///    |  100     100      100      100
///   0- * [0]    * [1]    * [2]    * [3]
///      0-------10-------20-------30
///
/// Output
///
///       4.4     7.9     11.3      14.8
///  20- * [8]   * [9]     * [10]   * [11]
///    |
///    |
///    |  1       4.4      7.9      11.3
///  10- * [4]    * [5]    * [6]    * [7]
///    |
///    |
///    |  4.4     5.9      9.3       12.8
///   0- * [0]    * [1]    * [2]    * [3]
///      0-------10-------20-------30
/// \endverbatim
//------------------------------------------------------------------------------
//! [snip_test_Example_SmoothSizeFunc]
void TutMeshingIntermediateTests::test_Example_SmoothSizeFunc()
{
  // create a grid of points
  BSHP<xms::VecPt3d> pts(new xms::VecPt3d());
  *pts = {{0, 0},   {10, 0},  {20, 0}, {30, 0},  {0, 10},  {10, 10},
          {20, 10}, {30, 10}, {0, 20}, {10, 20}, {20, 20}, {30, 20}};
  // assign all points a size of 100
  xms::VecFlt sizes(pts->size(), 100);
  // change the size to 1.0 of the point at 0, 10
  sizes[4] = 1;
  // create a TrTin class from the array of points
  BSHP<xms::VecInt> tris(new xms::VecInt());
  BSHP<xms::VecInt2d> adjTris(new xms::VecInt2d());
  xms::TrTriangulatorPoints tr(*pts, *tris, &*adjTris);
  tr.Triangulate();
  BSHP<xms::TrTin> tin = xms::TrTin::New();
  tin->SetGeometry(pts, tris, adjTris);

  // smooth the size function. The size ratio is set to 0.5. The min element
  // size is 1.0. The "anchor type" is 0 (meaning min). This means the minimum
  // size will be honored and the other values smoothed from the min.
  xms::VecFlt vSmooth;
  xms::DynBitset ptFlgs;
  xms::meSmoothSizeFunction(tin, sizes, 0.5, 1.0, 0, ptFlgs, vSmooth);
  xms::VecFlt baseSmooth = {4.46f, 5.90f,  9.36f, 12.83f, 1.0f,   4.46f,
                            7.93f, 11.39f, 4.46f, 7.93f,  11.39f, 14.86f};
  // ensure the results are as expected
  TS_ASSERT_DELTA_VEC(baseSmooth, vSmooth, .1);
} // TutMeshingIntermediateTests::test_Example_SmoothSizeFunc
  //! [snip_test_Example_SmoothSizeFunc]
//------------------------------------------------------------------------------
/// \brief Example for meshing a simple polygon with a hole and using the
/// spring relaxation method
//------------------------------------------------------------------------------
//! [snip_test_Example_SpringRelax]
void TutMeshingIntermediateTests::test_Example_SpringRelax()
{
  // The MePolyInput class defines a single polygon to be meshed and which
  // consists of one outer polygon and, optionally, inner polygons and some
  // other options.
  xms::MePolyInput inputPoly;

  // Outer polygon points are in clockwise order. Do not repeat the first point.
  inputPoly.m_outPoly = {
    {0, 10},   {0, 20},   {0, 30},   {0, 40},    {0, 50},   {0, 60},   {0, 70},   {0, 80},
    {0, 90},   {0, 100},  {10, 100}, {20, 100},  {30, 100}, {40, 100}, {50, 100}, {60, 100},
    {70, 100}, {80, 100}, {90, 100}, {100, 100}, {100, 90}, {100, 80}, {100, 70}, {100, 60},
    {100, 50}, {100, 40}, {100, 30}, {100, 20},  {100, 10}, {100, 0},  {90, 0},   {80, 0},
    {70, 0},   {60, 0},   {50, 0},   {40, 0},    {30, 0},   {20, 0},   {10, 0},   {0, 0}};

  // Inner polygons are in counter clockwise order. Do not repeat the first
  // point.
  inputPoly.m_insidePolys.push_back(xms::VecPt3d());
  inputPoly.m_insidePolys[0] = {{40, 40}, {50, 40}, {60, 40}, {60, 50},
                                {60, 60}, {50, 60}, {40, 60}, {40, 50}};
  inputPoly.m_relaxationMethod = "spring_relaxation";

  // The MeMultiPolyMesherIo class holds all the options for generating
  // UGrids from polygon data
  xms::MeMultiPolyMesherIo input;

  // The m_polys vector holds the polygons that will be filled with a UGrid.
  // This case has only 1 polygon.
  input.m_polys.push_back(inputPoly);

  // generate the mesh and check the base line
  const std::string baseFile = "Example_SpringRelax";
  tutGenerateAndCompare2dm(input, baseFile);
} // TutMeshingIntermediateTests::test_Example_SpringRelax
//! [snip_test_Example_SpringRelax]

//------------------------------------------------------------------------------
/// \brief Example for redistributing points on a polygon boundary to a constant
/// spacing.
//------------------------------------------------------------------------------
//! [snip_test_example_SimplePolygon_Redistribute]
void TutRedistributionIntermediateTests::test_Example_SimplePolygon_Redistribute()
{
  // Outer polygon points are in clockwise order
  xms::VecPt3d polygon = {{0, 0},    {0, 10},   {0, 20},    {0, 30},   {0, 40},   {0, 50},
                          {0, 60},   {0, 70},   {0, 80},    {0, 90},   {0, 100},  {10, 100},
                          {20, 100}, {30, 100}, {40, 100},  {50, 100}, {60, 100}, {70, 100},
                          {80, 100}, {90, 100}, {100, 100}, {100, 90}, {100, 80}, {100, 70},
                          {100, 60}, {100, 50}, {100, 40},  {100, 30}, {100, 20}, {100, 10},
                          {100, 0},  {90, 0},   {80, 0},    {70, 0},   {60, 0},   {50, 0},
                          {40, 0},   {30, 0},   {20, 0},    {10, 0},   {0, 0}};
  // create the redistribution class
  BSHP<xms::MePolyRedistributePts> redist = xms::MePolyRedistributePts::New();
  // set the redistribution class to use a constant spacing
  redist->SetConstantSizeFunc(20.0);
  // redistribute the points
  xms::VecPt3d outPts = redist->Redistribute(polygon);
  {
    xms::VecPt3d expectedPts = {
      {0, 0, 0},    {0, 20, 0},   {0, 40, 0},   {0, 60, 0},   {0, 80, 0},    {0, 100, 0},
      {20, 100, 0}, {40, 100, 0}, {60, 100, 0}, {80, 100, 0}, {100, 100, 0}, {100, 80, 0},
      {100, 60, 0}, {100, 40, 0}, {100, 20, 0}, {100, 0, 0},  {80, 0, 0},    {60, 0, 0},
      {40, 0, 0},   {20, 0, 0},   {0, 0, 0}};
    TS_ASSERT_DELTA_VECPT3D(expectedPts, outPts, 1e-3);
  }
} // TutRedistributionIntermediateTests::test_Example_SimplePolygon_Redistribute
//! [snip_test_example_SimplePolygon_Redistribute]
//------------------------------------------------------------------------------
/// \brief Example for redistributing a polygon using scalar paving.
//------------------------------------------------------------------------------
//! [snip_test_Example_Redistribute_SizeFunction]
void TutRedistributionIntermediateTests::test_Example_Redistribute_SizeFunction()
{
  // Outer polygon points are in clockwise order
  xms::VecPt3d polygon = {{0, 0},    {0, 10},   {0, 20},    {0, 30},   {0, 40},   {0, 50},
                          {0, 60},   {0, 70},   {0, 80},    {0, 90},   {0, 100},  {10, 100},
                          {20, 100}, {30, 100}, {40, 100},  {50, 100}, {60, 100}, {70, 100},
                          {80, 100}, {90, 100}, {100, 100}, {100, 90}, {100, 80}, {100, 70},
                          {100, 60}, {100, 50}, {100, 40},  {100, 30}, {100, 20}, {100, 10},
                          {100, 0},  {90, 0},   {80, 0},    {70, 0},   {60, 0},   {50, 0},
                          {40, 0},   {30, 0},   {20, 0},    {10, 0},   {0, 0}};

  // create a size function interpolator

  // These are the interpolator point locations. The size is specified by the
  // "z" component of the points.
  BSHP<xms::VecPt3d> sPts(new xms::VecPt3d());
  *sPts = {{-10, -10, 10}, {-10, 110, 10}, {110, 110, 1}, {110, -10, 10}, {60, 70, 1}};
  // These are the interpolator triangles.
  // Triangles are specified as point indexes in ccw order.  You can see the
  // 4 triangles in the vector below.
  BSHP<xms::VecInt> sTris(new xms::VecInt());
  *sTris = {0, 4, 1, 1, 4, 2, 2, 4, 3, 3, 4, 0};
  // create a linear interpolator for the size function
  BSHP<xms::InterpBase> linear(xms::InterpLinear::New());
  // sets the points and the triangles for the interpolator
  linear->SetPtsTris(sPts, sTris);
  // create the redistribution class
  BSHP<xms::MePolyRedistributePts> redist = xms::MePolyRedistributePts::New();
  // set the redistribution class to use a constant spacing
  redist->SetSizeFunc(linear);
  // redistribute the points
  xms::VecPt3d outPts = redist->Redistribute(polygon);
  {
    xms::VecPt3d expectedPts = {
      {0.000, 0.000, 0},    {0.000, 8.794, 0},    {0.000, 17.574, 0},   {0.000, 26.355, 0},
      {0.000, 35.135, 0},   {0.000, 43.916, 0},   {0.000, 52.697, 0},   {0.000, 61.477, 0},
      {0.000, 70.258, 0},   {0.000, 79.038, 0},   {0.000, 87.819, 0},   {0.000, 96.599, 0},
      {5.230, 100.000, 0},  {13.015, 100.000, 0}, {20.214, 100.000, 0}, {26.898, 100.000, 0},
      {33.102, 100.000, 0}, {38.864, 100.000, 0}, {44.212, 100.000, 0}, {49.181, 100.000, 0},
      {53.794, 100.000, 0}, {58.079, 100.000, 0}, {62.059, 100.000, 0}, {65.757, 100.000, 0},
      {69.192, 100.000, 0}, {72.382, 100.000, 0}, {75.348, 100.000, 0}, {78.104, 100.000, 0},
      {80.666, 100.000, 0}, {83.046, 100.000, 0}, {85.260, 100.000, 0}, {87.319, 100.000, 0},
      {89.234, 100.000, 0}, {91.015, 100.000, 0}, {92.672, 100.000, 0}, {94.214, 100.000, 0},
      {95.649, 100.000, 0}, {96.986, 100.000, 0}, {98.230, 100.000, 0}, {99.390, 100.000, 0},
      {100.000, 99.523, 0}, {100.000, 98.406, 0}, {100.000, 97.223, 0}, {100.000, 95.963, 0},
      {100.000, 94.609, 0}, {100.000, 93.156, 0}, {100.000, 91.594, 0}, {100.000, 89.915, 0},
      {100.000, 88.111, 0}, {100.000, 86.172, 0}, {100.000, 84.087, 0}, {100.000, 81.844, 0},
      {100.000, 79.431, 0}, {100.000, 76.836, 0}, {100.000, 74.044, 0}, {100.000, 71.040, 0},
      {100.000, 67.806, 0}, {100.000, 64.327, 0}, {100.000, 60.580, 0}, {100.000, 56.547, 0},
      {100.000, 52.206, 0}, {100.000, 47.530, 0}, {100.000, 42.497, 0}, {100.000, 37.076, 0},
      {100.000, 31.239, 0}, {100.000, 24.951, 0}, {100.000, 18.178, 0}, {100.000, 10.884, 0},
      {100.000, 3.028, 0},  {94.462, 0.000, 0},   {85.520, 0.000, 0},   {76.579, 0.000, 0},
      {67.638, 0.000, 0},   {58.696, 0.000, 0},   {49.755, 0.000, 0},   {40.814, 0.000, 0},
      {31.873, 0.000, 0},   {22.931, 0.000, 0},   {13.990, 0.000, 0},   {5.049, 0.000, 0},
      {0.000, 0.000, 0}};
    TS_ASSERT_DELTA_VECPT3D(expectedPts, outPts, 1e-3);
  }

} // TutMeshingIntermediateTests::test_Example_ScalarPaving
//! [snip_test_Example_Redistribute_SizeFunction]
//------------------------------------------------------------------------------
/// \brief Example for redistributing a polyline using curvature.
//------------------------------------------------------------------------------
//! [test_Example_Redistribute_Curvature]
void TutRedistributionIntermediateTests::test_Example_Redistribute_Curvature()
{
  xms::VecPt3d polyline = {{0, 0, 0},   {5, 5, 0},   {10, 10, 0}, {15, 5, 0},
                           {20, 10, 0}, {21, 10, 0}, {25, 0, 0}};
  // create the redistribution class
  BSHP<xms::MePolyRedistributePts> redist = xms::MePolyRedistributePts::New();
  // set the redistribution class to curvature
  double featureSize(3.0), meanSpacing(0.5), minimumCurvature(.0001);
  bool smooth(false);
  redist->SetUseCurvatureRedistribution(featureSize, meanSpacing, minimumCurvature, smooth);
  // redistribute the points
  xms::VecPt3d outPts = redist->Redistribute(polyline);
  {
    xms::VecPt3d expectedPts = {
      {0.000, 0.000, 0},   {9.526, 9.526, 0},   {9.582, 9.582, 0},   {9.639, 9.639, 0},
      {9.695, 9.695, 0},   {9.751, 9.751, 0},   {9.808, 9.808, 0},   {9.864, 9.864, 0},
      {9.921, 9.921, 0},   {9.977, 9.977, 0},   {10.034, 9.966, 0},  {10.090, 9.910, 0},
      {10.146, 9.854, 0},  {10.203, 9.797, 0},  {10.259, 9.741, 0},  {10.316, 9.684, 0},
      {10.372, 9.628, 0},  {10.429, 9.571, 0},  {10.485, 9.515, 0},  {14.481, 5.519, 0},
      {14.537, 5.463, 0},  {14.594, 5.406, 0},  {14.650, 5.350, 0},  {14.707, 5.293, 0},
      {14.763, 5.237, 0},  {14.819, 5.181, 0},  {14.876, 5.124, 0},  {14.932, 5.068, 0},
      {14.989, 5.011, 0},  {15.045, 5.045, 0},  {15.102, 5.102, 0},  {15.158, 5.158, 0},
      {15.215, 5.215, 0},  {15.271, 5.271, 0},  {15.327, 5.327, 0},  {15.384, 5.384, 0},
      {15.440, 5.440, 0},  {15.497, 5.497, 0},  {19.484, 9.484, 0},  {19.518, 9.518, 0},
      {19.552, 9.552, 0},  {19.587, 9.587, 0},  {19.621, 9.621, 0},  {19.655, 9.655, 0},
      {19.690, 9.690, 0},  {19.724, 9.724, 0},  {19.759, 9.759, 0},  {19.793, 9.793, 0},
      {19.827, 9.827, 0},  {19.862, 9.862, 0},  {19.896, 9.896, 0},  {19.930, 9.930, 0},
      {19.965, 9.965, 0},  {19.999, 9.999, 0},  {20.047, 10.000, 0}, {20.096, 10.000, 0},
      {20.144, 10.000, 0}, {20.193, 10.000, 0}, {20.242, 10.000, 0}, {20.790, 10.000, 0},
      {20.838, 10.000, 0}, {20.886, 10.000, 0}, {20.934, 10.000, 0}, {20.982, 10.000, 0},
      {21.011, 9.972, 0},  {21.029, 9.928, 0},  {21.047, 9.883, 0},  {21.065, 9.839, 0},
      {21.082, 9.794, 0},  {21.100, 9.749, 0},  {21.118, 9.705, 0},  {21.136, 9.660, 0},
      {21.154, 9.616, 0},  {21.172, 9.571, 0},  {21.189, 9.527, 0},  {21.207, 9.482, 0},
      {21.225, 9.437, 0},  {21.243, 9.393, 0},  {21.261, 9.348, 0},  {25.000, 0.000, 0}};
    TS_ASSERT_DELTA_VECPT3D(expectedPts, outPts, 1e-3);
  }
} // TutRedistributionIntermediateTests::test_Example_Redistribute_Curvature
//! [test_Example_Redistribute_Curvature]
//------------------------------------------------------------------------------
/// \brief Example for redistributing a polyline using curvature.
//------------------------------------------------------------------------------
//! [test_Example_Redistribute_Polygon_Curvature]
void TutRedistributionIntermediateTests::test_Example_Redistribute_Polygon_Curvature()
{
  xms::VecPt3d polyline = {{0, 0, 0},   {5, 5, 0},   {10, 10, 0}, {15, 5, 0},
                           {20, 10, 0}, {21, 10, 0}, {25, 0, 0},  {0, 0, 0}};
  // create the redistribution class
  BSHP<xms::MePolyRedistributePts> redist = xms::MePolyRedistributePts::New();
  // set the redistribution class to curvature
  double featureSize(3.0), meanSpacing(0.5), minimumCurvature(.0001);
  bool smooth(false);
  redist->SetUseCurvatureRedistribution(featureSize, meanSpacing, minimumCurvature, smooth);
  // redistribute the points
  xms::VecPt3d outPts = redist->Redistribute(polyline);
  {
    xms::VecPt3d expectedPts = {
      {0.000, 0.000, 0},   {0.042, 0.042, 0},   {0.084, 0.084, 0},   {0.126, 0.126, 0},
      {0.168, 0.168, 0},   {0.210, 0.210, 0},   {0.252, 0.252, 0},   {0.294, 0.294, 0},
      {0.336, 0.336, 0},   {0.378, 0.378, 0},   {0.420, 0.420, 0},   {0.462, 0.462, 0},
      {0.505, 0.505, 0},   {9.491, 9.491, 0},   {9.546, 9.546, 0},   {9.601, 9.601, 0},
      {9.655, 9.655, 0},   {9.710, 9.710, 0},   {9.765, 9.765, 0},   {9.820, 9.820, 0},
      {9.875, 9.875, 0},   {9.930, 9.930, 0},   {9.985, 9.985, 0},   {10.040, 9.960, 0},
      {10.095, 9.905, 0},  {10.150, 9.850, 0},  {10.205, 9.795, 0},  {10.260, 9.740, 0},
      {10.315, 9.685, 0},  {10.370, 9.630, 0},  {10.424, 9.576, 0},  {10.479, 9.521, 0},
      {14.474, 5.526, 0},  {14.529, 5.471, 0},  {14.584, 5.416, 0},  {14.639, 5.361, 0},
      {14.693, 5.307, 0},  {14.748, 5.252, 0},  {14.803, 5.197, 0},  {14.858, 5.142, 0},
      {14.913, 5.087, 0},  {14.968, 5.032, 0},  {15.023, 5.023, 0},  {15.078, 5.078, 0},
      {15.133, 5.133, 0},  {15.188, 5.188, 0},  {15.243, 5.243, 0},  {15.298, 5.298, 0},
      {15.353, 5.353, 0},  {15.408, 5.408, 0},  {15.462, 5.462, 0},  {15.517, 5.517, 0},
      {19.495, 9.495, 0},  {19.529, 9.529, 0},  {19.562, 9.562, 0},  {19.596, 9.596, 0},
      {19.629, 9.629, 0},  {19.662, 9.662, 0},  {19.696, 9.696, 0},  {19.729, 9.729, 0},
      {19.763, 9.763, 0},  {19.796, 9.796, 0},  {19.830, 9.830, 0},  {19.863, 9.863, 0},
      {19.897, 9.897, 0},  {19.930, 9.930, 0},  {19.964, 9.964, 0},  {19.997, 9.997, 0},
      {20.043, 10.000, 0}, {20.090, 10.000, 0}, {20.138, 10.000, 0}, {20.185, 10.000, 0},
      {20.232, 10.000, 0}, {20.779, 10.000, 0}, {20.826, 10.000, 0}, {20.873, 10.000, 0},
      {20.919, 10.000, 0}, {20.966, 10.000, 0}, {21.005, 9.988, 0},  {21.022, 9.945, 0},
      {21.039, 9.901, 0},  {21.057, 9.858, 0},  {21.074, 9.815, 0},  {21.092, 9.771, 0},
      {21.109, 9.728, 0},  {21.126, 9.684, 0},  {21.144, 9.641, 0},  {21.161, 9.598, 0},
      {21.178, 9.554, 0},  {21.196, 9.511, 0},  {21.213, 9.467, 0},  {21.230, 9.424, 0},
      {21.248, 9.381, 0},  {21.265, 9.337, 0},  {24.727, 0.682, 0},  {24.752, 0.621, 0},
      {24.776, 0.559, 0},  {24.801, 0.498, 0},  {24.826, 0.436, 0},  {24.850, 0.375, 0},
      {24.875, 0.313, 0},  {24.899, 0.251, 0},  {24.924, 0.190, 0},  {24.949, 0.128, 0},
      {24.973, 0.067, 0},  {24.998, 0.005, 0},  {24.939, 0.000, 0},  {24.873, 0.000, 0},
      {24.806, 0.000, 0},  {24.740, 0.000, 0},  {24.674, 0.000, 0},  {24.607, 0.000, 0},
      {24.541, 0.000, 0},  {24.475, 0.000, 0},  {24.408, 0.000, 0},  {24.342, 0.000, 0},
      {24.276, 0.000, 0},  {0.713, 0.000, 0},   {0.654, 0.000, 0},   {0.595, 0.000, 0},
      {0.535, 0.000, 0},   {0.476, 0.000, 0},   {0.416, 0.000, 0},   {0.357, 0.000, 0},
      {0.297, 0.000, 0},   {0.238, 0.000, 0},   {0.178, 0.000, 0},   {0.119, 0.000, 0},
      {0.059, 0.000, 0},   {0.000, 0.000, 0}};
    TS_ASSERT_DELTA_VECPT3D(expectedPts, outPts, 1e-3);
  }
} // TutRedistributionIntermediateTests::test_Example_Redistribute_Polygon_Curvature
//! [test_Example_Redistribute_Polygon_Curvature]

//------------------------------------------------------------------------------
/// \brief Example for using the MeQuadBlossom class to convert a triangle mesh
/// into a mesh of quads and then to use MeBadQuadRemover to remove narrow,
/// ill-formed quads while still preserving the basic sizes of cells.
//------------------------------------------------------------------------------
//! [snip_test_Example_QuadBlossom_BadQuadRemover]
#include <xmsgrid/ugrid/XmUGrid.h>
#include <xmsgrid/ugrid/XmUGridUtils.h>
#include <xmsmesh/meshing/MeQuadBlossom.h>
#include <xmsmesh/meshing/MeBadQuadRemover.h>

void TutMeshingIntermediateTests::test_Example_QuadBlossom_BadQuadRemover()
{
  // read a UGrid from a file.
  const std::string path(std::string(XMS_TEST_PATH) + 
                         "Tutorial_Meshing/");
  const std::string inputFilePath = path + "Example_QuadBlossom_triangleUGridInput.txt";
  BSHP<xms::XmUGrid> ugrid = xms::XmReadUGridFromAsciiFile(inputFilePath); // read from file.
  BSHP<xms::MeQuadBlossom> quadBlossom = xms::MeQuadBlossom::New(ugrid);

  // PreMakeQuads returns the number of boundary edges of the triangular mesh.
  // An even number of boundary edges insures no triangles in the output UGrid. 
  int nBoundaryEdges = quadBlossom->PreMakeQuads();
  TS_ASSERT((nBoundaryEdges & 0x1) == 0); 

  // The MeQuadBlossom basic algorithm is O(N^3). Check the estimated minutes.  If too large,
  // Then split the mesh into smaller sub-UGrids, then call MakeQuads on each one.
  double minutes = xms::MeQuadBlossom::EstimatedRunTimeInMinutes(ugrid->GetPointCount());
  TS_ASSERT(minutes < 2.0);

  bool splitVertices = true;
  bool useAngle = false;
  BSHP<xms::XmUGrid> quadUGrid = quadBlossom->MakeQuads(splitVertices, useAngle);

  // Test the quad UGrid generated by the Quad Blossom algorithm
  const std::string outFile = path + "Example_QuadBlossom_quadUGrid_out.txt";
  xms::XmWriteUGridToAsciiFile(quadUGrid, outFile);
  const std::string baseFile = path + "Example_QuadBlossom_quadUGrid_base.txt";
  TS_ASSERT_TXT_FILES_EQUAL(baseFile, outFile);

  // Use MeBadQuadRemover to remove bad quads from the quad UGrid.
  double maxAspect = 0.7;
  BSHP<xms::MeBadQuadRemover> badQuadRemover = xms::MeBadQuadRemover::New(quadUGrid);
  BSHP<xms::XmUGrid> quadUGridImproved = badQuadRemover->RemoveBadQuads(maxAspect);

  // Test the improved quad UGrid returned from the Bad Quad Remover algorithm
  const std::string outFile2 = path + "Example_QuadBlossom_quadUGridImproved_out.txt";
  xms::XmWriteUGridToAsciiFile(quadUGridImproved, outFile2);
  const std::string baseFile2 = path + "Example_QuadBlossom_quadUGridImproved_base.txt";
  TS_ASSERT_TXT_FILES_EQUAL(baseFile2, outFile2);

  // Use MeBadQuadRemover to remove bad quads from the quad UGrid a second time.
  maxAspect = 0.7;
  badQuadRemover = xms::MeBadQuadRemover::New(quadUGridImproved);
  BSHP<xms::XmUGrid> quadUGridImproved2 = badQuadRemover->RemoveBadQuads(maxAspect);

  // Test the improved quad UGrid returned from the Bad Quad Remover algorithm
  const std::string outFile3 = path + "Example_QuadBlossom_quadUGridImproved2_out.txt";
  xms::XmWriteUGridToAsciiFile(quadUGridImproved2, outFile3);
  const std::string baseFile3 = path + "Example_QuadBlossom_quadUGridImproved2_base.txt";
  TS_ASSERT_TXT_FILES_EQUAL(baseFile3, outFile3);
}  // TutMeshingIntermediateTests::test_Example_QuadBlossom_BadQuadRemover
//! [snip_test_Example_QuadBlossom_BadQuadRemover]
#endif
