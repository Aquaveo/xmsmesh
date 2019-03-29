#pragma once
#ifdef CXX_TEST
//------------------------------------------------------------------------------
/// \file
/// \ingroup meshing
/// \copyright (C) Copyright Aquaveo 2018. Distributed under FreeBSD License
/// (See accompanying file LICENSE or https://aqaveo.com/bsd/license.txt)
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
  void testCreateClass();
  void testCase4();
};

class MeMultiPolyTo2dmIntermediateTests : public CxxTest::TestSuite
{
public:
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
  void testbug11299();
  void testInternalFeaturesCase0();
  void testInternalFeaturesCase1();
  void testInternalFeaturesCase2();
  void testbug11646();
};

//} // namespace xms
#endif
