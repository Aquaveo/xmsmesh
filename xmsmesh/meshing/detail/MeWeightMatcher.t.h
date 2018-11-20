#pragma once
#ifdef CXX_TEST
//------------------------------------------------------------------------------
/// \file
/// \ingroup meshing_detail
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

class MeWeightMatcherUnitTests : public CxxTest::TestSuite
{
public:
  void testPythonVectors();
  void test10_empty();
  void test11_singleedge();
  void test12();
  void test13();
  void test14_maxcard();
  void test16_negative();
  void test20_sblossom();
  void test21_tblossom();
  void test22_s_nest();
  void test23_s_relabel_nest();
  void test24_s_nest_expand();
  void test25_s_t_expand();
  void test26_s_nest_t_expand();
  void test30_tnasty_expand();
  void test31_tnasty2_expand();
  void test32_t_expand_leastslack();
  void test33_nest_tnasty_expand();
  void test34_nest_relabel_expand();
  void testSimpleTriangle();
  void testSimpleQuad();
  void testComplexQuad();
};

//} // namespace xms
#endif
