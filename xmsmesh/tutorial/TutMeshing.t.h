#pragma once
#ifdef CXX_TEST
//------------------------------------------------------------------------------
/// \file
/// \ingroup meshing
/// \copyright (C) Copyright Aquaveo 2018. Distributed under the xmsng
///  Software License, Version 1.0. (See accompanying file
///  LICENSE_1_0.txt or copy at http://www.aquaveo.com/xmsng/LICENSE_1_0.txt)
//------------------------------------------------------------------------------

// 3. Standard Library Headers

// 4. External Library Headers
#include <cxxtest/TestSuite.h>

// 5. Shared Headers
#include <xmscore/stl/vector.h>

// 6. Non-shared Headers

//----- Namespace declaration --------------------------------------------------

namespace xms
{
class MeMultiPolyMesherIo;

bool tutReadMeshIoFromFile(const std::string& a_fname, MeMultiPolyMesherIo& a_io);
bool tutReadPolygons(const std::string& a_fname, VecPt3d2d& a_outside, VecPt3d3d& a_inside);
} // namespace xms

/// \brief Class for testing meshing functionality
class TutMeshingIntermediateTests : public CxxTest::TestSuite
{
public:
#ifndef CXXTEST4
  virtual const CxxTest::TestGroup& group();
#endif
  void test_Example_SimplePolygon();
  void test_Example_ComplexPolygon();
  void test_Example_SimplePolygonWithHole();
  void test_Example_Breakline();
  void test_Example_RefinePoints();
  void test_Example_MultiplePolygons();
  void test_Example_ScalarPaving();
  void test_Example_Patch();
  void test_Example_ConstantSmooth();
  void test_Example_SizeFuncFromDepth();
  void test_Example_SmoothSizeFunc();
  void test_Example_SpringRelax();
};

/// \brief Class for testing polyline/polygon point redistribution
class TutRedistributionIntermediateTests : public CxxTest::TestSuite
{
public:
  void test_Example_SimplePolygon_Redistribute();
  void test_Example_Redistribute_SizeFunction();
  void test_Example_Redistribute_Curvature();
  void test_Example_Redistribute_Polygon_Curvature();
};
#endif
