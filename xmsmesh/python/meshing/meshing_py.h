#pragma once
//------------------------------------------------------------------------------
/// \file
/// \brief initializer functions for members of meshing python module.
/// \copyright (C) Copyright Aquaveo 2018. Distributed under FreeBSD License
/// (See accompanying file LICENSE or https://aqaveo.com/bsd/license.txt)
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
