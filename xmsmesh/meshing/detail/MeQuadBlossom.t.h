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

class MeQuadBlossomUnitTests : public CxxTest::TestSuite
{
public:
  void testGetInteriorEdges();
  void testGetBoundaryEdges();
  void testProcessVertexChains();
  void testEta();
  void testGetAngle();
  void testExtractCellsAdjacentToPoint();
  void testGetEdges();
  void testEliminateEdges();
  void testSimpleTriangle();
  void testSimpleQuad();
  void testComplexQuad();
  void testSplitToQuads();
  void testEstimatedRunTime();
  void testPreMakeQuads();
};

//} // namespace xms
#endif
