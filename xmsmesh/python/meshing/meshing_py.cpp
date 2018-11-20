//------------------------------------------------------------------------------
/// \file
/// \brief
/// \copyright (C) Copyright Aquaveo 2018. Distributed under FreeBSD License
/// (See accompanying file LICENSE or https://aqaveo.com/bsd/license.txt)
//------------------------------------------------------------------------------

//----- Included files ---------------------------------------------------------
#include <pybind11/pybind11.h>
#include <xmsmesh/python/meshing/meshing_py.h>

//----- Namespace declaration --------------------------------------------------
namespace py = pybind11;

//----- Python Interface -------------------------------------------------------

void initMeshing(py::module &m) {
    initMeBadQuadRemover(m);
    initMeMeshUtils(m);
    initMeMultiPolyMesher(m);
    initMePolyInput(m);
    initMeRefinePoint(m);
    initMeMultiPolyMesherIo(m);
    initMeMultiPolyTo2dm(m);
    initMePolyMesher(m);
    initMePolyRedistributePts(m);
    initMeQuadBlossom(m);
}