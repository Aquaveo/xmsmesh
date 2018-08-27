#pragma once
#ifdef CXX_TEST
//------------------------------------------------------------------------------
/// \file
/// \ingroup meshing_detail
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

////////////////////////////////////////////////////////////////////////////////
class MeRelaxerUnitTests : public CxxTest::TestSuite
{
public:
  void testRelaxWhileMeshing();
  void testSpringRelaxSetup();
  void testSpringRelaxSetup2();
  void testSpringRelaxSinglePoint();
  void testSpringRelaxSinglePoint2();
  void testSpringRelaxSinglePoint3();
  void testNewLocationIsValid();
  void testAllTrianglesHavePositiveArea();
};

//} // namespace xms
#endif
