//------------------------------------------------------------------------------
/// \file
/// \ingroup meshing
/// \copyright (C) Copyright Aquaveo 2018. Distributed under the xmsng
///  Software License, Version 1.0. (See accompanying file
///  LICENSE_1_0.txt or copy at http://www.aquaveo.com/xmsng/LICENSE_1_0.txt)
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
