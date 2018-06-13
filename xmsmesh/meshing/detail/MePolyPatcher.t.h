//------------------------------------------------------------------------------
/// \file
/// \ingroup meshing_detail
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

//----- Constants / Enumerations -----------------------------------------------

//----- Structs / Classes ------------------------------------------------------
/// \brief tests MePolyPatcher
/// \see MePolyPatcher
class MePolyPatcherUnitTests : public CxxTest::TestSuite
{
public:
  void testCreateClass();
  void testPolyPtsToSides();
  void testPolyPtsToSides1();
  void testPatch00();
  void testQuadPatchErrors();
  void testBug9226();
};
//----- Function prototypes ----------------------------------------------------

#endif
