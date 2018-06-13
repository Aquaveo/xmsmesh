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
// namespace xms {

//----- Constants / Enumerations -----------------------------------------------

//----- Structs / Classes ------------------------------------------------------
class MePolyCleanerUnitTests : public CxxTest::TestSuite
{
public:
  void testCreateClass();
  void testCase0();
  void testCase1();
  void testCase2();
  void testCase2a();
  void testCase3();
  void testCase3a();
  void testCase4();
  void testCase5();
  void testCase5a();
  void testCase6();
  void testCase6a();
  void testCase6b();
  void testCase6c();
  void testCase7();
  void testCase7a();
  void testCase7b();
  void testCase7c();
  void testCase8();
  void testCase8a();
  void testCase8b();
  void testCase8c();
  void testCase9();

  void testCleanIn0();
  void testCleanIn1();
  void testCleanIn1a();
  void testCleanIn1b();
  void testCleanIn2();
  void testCleanIn2a();
  void testCleanIn2b();
  void testCleanIn3();
  void testCleanIn3a();
  void testCleanIn4();
  void testCleanIn4a();
  void testCleanIn4b();
  void testCleanIn4c();
  void testCleanIn4d();
  void testCleanIn5();
  void testCleanIn6();
  void testCleanIn7();
  void testCleanIn8();
  void testCleanIn9();
  void testCleanIn10();
  void testCleanIn11();

  void testCleanInOut0();
  void testCleanInOut0a();
  void testCleanInOut0b();
  void testCleanInOut1();
  void testCleanInOut1a();
  void testCleanInOut1b();
};
//----- Function prototypes ----------------------------------------------------

//} // namespace xms
#endif // CXX_TEST
