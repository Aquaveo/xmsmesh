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
#include <xmsmesh/meshing/MeMeshUtils.h>
#include <xmsmesh/meshing/MeMultiPolyMesherIo.h>
#include <xmsmesh/meshing/MeMultiPolyTo2dm.h>
#include <xmsmesh/meshing/MePolyRedistributePts.h>
#include <xmsinterp/triangulate/TrTin.h>
#include <xmsinterp/triangulate/TrTriangulatorPoints.h>

namespace xms
{
//------------------------------------------------------------------------------
/// \brief  helper function to read MeMultiPolyMesherIo
/// \param a_fname File name where the polygon is stored
/// \param a_io MeMultiPolyMesherIo class that is filled by this function
//------------------------------------------------------------------------------
void tutReadMeshIoFromFile(const std::string& a_fname, MeMultiPolyMesherIo& a_io)
{
  a_io = MeMultiPolyMesherIo();
  std::fstream os(a_fname, std::fstream::in);
  if (!os.is_open())
    return;

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
    card = "";
  }
} // tutReadMeshIoFromFile
//------------------------------------------------------------------------------
/// \brief  helper function to read polygons from a text file
/// \param a_fname File name where the polygon is stored
/// \param a_outside Vector of locations for the definition of the outside of
/// the polygon.
/// \param a_inside 2d vector of locations to define multiple inside polygons
//------------------------------------------------------------------------------
void tutReadPolygons(const std::string& a_fname, VecPt3d2d& a_outside, VecPt3d3d& a_inside)
{
  a_outside.resize(0);
  a_inside.resize(0);

  MeMultiPolyMesherIo io;
  tutReadMeshIoFromFile(a_fname, io);
  for (size_t i = 0; i < io.m_polys.size(); ++i)
  {
    a_outside.push_back(io.m_polys[i].m_outPoly);
    a_inside.push_back(io.m_polys[i].m_insidePolys);
  }
} // tutReadPolygons
//------------------------------------------------------------------------------
/// \brief  helper function to generate 2dm file and compare to a baseline
/// \param a_input poly mesher input
/// \param a_baseFilenameWithPath base file name for output file and base line file
//------------------------------------------------------------------------------
void tutGenerateAndCompare2dm(MeMultiPolyMesherIo& a_io, const std::string& a_baseFilenameWithPath)
{
  std::string outFile = a_baseFilenameWithPath + "_out.2dm";
  std::string baseFile = a_baseFilenameWithPath + "_base.2dm";
  {
    std::fstream os(outFile.c_str(), std::fstream::out);
    TS_ASSERT(!os.bad());
    if (os.bad())
      return;
    BSHP<MeMultiPolyTo2dm> pm = MeMultiPolyTo2dm::New();
    pm->Generate2dm(a_io, os);
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
  const std::string path(xms::ttGetXmsngTestPath() + "/Tutorial_Meshing/");
  const std::string baseFile = path + "Example_SimplePolygon";
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
  const std::string path(xms::ttGetXmsngTestPath() + "/Tutorial_Meshing/");
  const std::string baseFile = path + "Example_ComplexPolygon";
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

  // Inner polygons are in counter clockwise order. Do not repeat the first point.
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
  const std::string path(xms::ttGetXmsngTestPath() + "/Tutorial_Meshing/");
  const std::string baseFile = path + "Example_SimplePolygonWithHole";
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

  // Inner polygons are in counter clockwise order. Do not repeat the first point.
  inputPoly.m_insidePolys.push_back(xms::VecPt3d());
  inputPoly.m_insidePolys[0] = {{50, 0}, {50, 10}, {50, 20}, {50, 10}};

  // The MeMultiPolyMesherIo class holds all the options for generating
  // UGrids from polygon data
  xms::MeMultiPolyMesherIo input;

  // The m_polys vector holds the polygons that will be filled with a UGrid.
  // This case has only 1 polygon.
  input.m_polys.push_back(inputPoly);

  // generate the mesh and check the base line
  const std::string path(xms::ttGetXmsngTestPath() + "/Tutorial_Meshing/");
  const std::string baseFile = path + "Example_Breakline";
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
  const std::string path(xms::ttGetXmsngTestPath() + "/Tutorial_Meshing/");
  const std::string baseFile = path + "Example_RefinePoints";
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
  const std::string path(xms::ttGetXmsngTestPath() + "/Tutorial_Meshing/");
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
  const std::string baseFile = path + "Example_MultiPolys";
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
  const std::string path(xms::ttGetXmsngTestPath() + "/Tutorial_Meshing/");
  const std::string baseFile = path + "Example_ScalarPaving";
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
  const std::string path(xms::ttGetXmsngTestPath() + "/Tutorial_Meshing/");
  const std::string baseFile = path + "Example_Patch";
  tutGenerateAndCompare2dm(input, baseFile);
} // TutMeshingIntermediateTests::test_Example_Patch
//! [snip_test_Example_Patch]
//------------------------------------------------------------------------------
/// \brief Example for redistributing points on a polyline
//------------------------------------------------------------------------------
//! [snip_test_Example_Redist]
void TutMeshingIntermediateTests::test_Example_PolyLineRedist()
{
  // create a polyline
  xms::VecPt3d pl = {{0, 0}, {100, 0}, {100, 100}, {0, 100}, {0, 0}};
  // create a point redistributor class
  BSHP<xms::MePolyRedistributePts> r = xms::MePolyRedistributePts::New();
  // set a constant size function
  r->SetConstantSizeFunc(25.0);

  // redistribute the poly line
  xms::VecPt3d p1 = r->Redistribute(pl);
  // verify the output
  xms::VecPt3d base = {{0, 0},    {25, 0},   {50, 0},    {75, 0},   {100, 0},  {100, 25},
                       {100, 50}, {100, 75}, {100, 100}, {75, 100}, {50, 100}, {25, 100},
                       {0, 100},  {0, 75},   {0, 50},    {0, 25},   {0, 0}};
  TS_ASSERT_DELTA_VEC_MP3(base, p1, 1e-9);
} // TutMeshingIntermediateTests::test_Example_PolyLineRedist
//! [snip_test_Example_Redist]
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
  const std::string path(xms::ttGetXmsngTestPath() + "/Tutorial_Meshing/");
  const std::string baseFile = path + "Example_ConstantSmooth";
  tutGenerateAndCompare2dm(input, baseFile);

  // now make the polygon grow to a bigger size element
  input.m_polys.front().m_constSizeFunction = 50;

  const std::string baseFile2 = path + "Example_ConstantSmooth2";
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

#endif