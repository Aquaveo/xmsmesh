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
class MePolyRedistributePtsUnitTests : public CxxTest::TestSuite
{
public:
  void testCreateClass();
  void testInterpEdgeLengths();
  void testInterpEdgeLengths1();
  void testInterpEdgeLengths2();
  void testInterpEdgeLengths3();
  void testInterpEdgeLengths4();
  void testRedistPts();
  void testRedistPts1();
  void testRedistPts2();
  void testRedistPts3();
  void testRedistPts4();
  void testRedistPts5();
  void testIntersectWithTris();
  void testRedistPolyLine();
  void testRedistPolyLine1();
  void testRedistPolyLine2();
};
//----- Function prototypes ----------------------------------------------------

//} // namespace xms
#endif
