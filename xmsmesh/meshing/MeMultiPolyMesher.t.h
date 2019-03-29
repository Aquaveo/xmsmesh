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

class MeMultiPolyMesherUnitTests : public CxxTest::TestSuite
{
public:
  void testCreateClass();
  void testCheckForIntersections1();
  void testCheckForIntersections2();
  void testCheckForIntersections3();
  void testCheckForIntersections4();
  void testCheckForIntersections5();
  void testEmptyOuterPolygon();
};

//} // namespace xms
#endif
