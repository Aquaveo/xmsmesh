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

// 6. Non-shared Headers

//----- Namespace declaration --------------------------------------------------

// namespace xms {

class MeMultiPolyTo2dmUnitTests : public CxxTest::TestSuite
{
public:
#ifndef CXXTEST4
  virtual const CxxTest::TestGroup& group();
#endif
  void testCreateClass();
  void testCase4();
};

class MeMultiPolyTo2dmIntermediateTests : public CxxTest::TestSuite
{
public:
#ifndef CXXTEST4
  virtual const CxxTest::TestGroup& group();
#endif
  void testCase2();
  void testCase100();
  void testCase101();
  void testCase102();
  void testCase103();
  void testCasePaveGeo();
  void testCasePaveSanDiego();
  void testCasePaveSanDiego_SpringRelaxation();
  void testCasePatch6();
  void testCasePaveConstSizeTransition();
  void testSeedPoints();
  void testSeedPoints_PolygonWithHole();
};

//} // namespace xms
#endif
