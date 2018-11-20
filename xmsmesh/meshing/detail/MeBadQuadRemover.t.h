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

class MeBadQuadRemoverUnitTests : public CxxTest::TestSuite
{
public:
  void testGetAdjacentPointCounts();
  void testReplacePoints();
  void testCollapse();
  void testCollapseQuadTri();
};

//} // namespace xms
#endif
