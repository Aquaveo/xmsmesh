#pragma once
//------------------------------------------------------------------------------
/// \file
/// \brief initializer functions for members of meshing python module.
/// \copyright (C) Copyright Aquaveo 2018. Distributed under the xmsng
///  Software License, Version 1.0. (See accompanying file
///  LICENSE_1_0.txt or copy at http://www.aquaveo.com/xmsng/LICENSE_1_0.txt)
//------------------------------------------------------------------------------

//----- Included files ---------------------------------------------------------
#include <pybind11/pybind11.h>

//----- Namespace declaration --------------------------------------------------
namespace py = pybind11;

//----- Function declarations --------------------------------------------------
void initMeshing(py::module &);

void initMeBadQuadRemover(py::module &);
void initMeMeshUtils(py::module &);
void initMeMultiPolyMesher(py::module &);
void initMeMultiPolyMesherIo(py::module &);
void initMeMultiPolyTo2dm(py::module &);
void initMePolyInput(py::module &);
void initMePolyMesher(py::module &);
void initMePolyRedistributePts(py::module &);
void initMeQuadBlossom(py::module &);
void initMeRefinePoint(py::module &);
