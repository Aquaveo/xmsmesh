//------------------------------------------------------------------------------
/// \file
/// \ingroup meshing
/// \copyright (C) Copyright Aquaveo 2018. Distributed under FreeBSD License
/// (See accompanying file LICENSE or https://aqaveo.com/bsd/license.txt)
//------------------------------------------------------------------------------
#pragma once
#ifdef CXX_TEST

//----- Included files ---------------------------------------------------------
// 3. Standard Library Headers

// 4. External Library Headers
#include <cxxtest/TestSuite.h>

// 5. Shared Headers

// 6. Non-shared Headers

//----- Forward declarations ---------------------------------------------------

//----- Namespace declaration --------------------------------------------------
// namespace xms {

//----- Constants / Enumerations -----------------------------------------------

//----- Structs / Classes ------------------------------------------------------
/// \brief Class for testing MePolyRedistributePtsCurvature
class MePolyRedistributePtsCurvatureUnitTests : public CxxTest::TestSuite
{
public:
  void testCreateClass();
  void testSetup();
  void testGetSignificantPoints();
  void testCalculateCurvature();
  void testShiftAndAggregateOpen();
  void testShiftAndAggregateClosed();
  void testDoSmoothing();
  void testNewPointsFromParamCurvs();
};

/// \brief Class for testing MePolyRedistributePtsCurvature
class MePolyRedistributePtsCurvatureIntermediateTests : public CxxTest::TestSuite
{
public:
  void testCoastline();
  void testIsland();
};

//----- Function prototypes ----------------------------------------------------

//} // namespace xms
#endif
