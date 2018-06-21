//------------------------------------------------------------------------------
/// \file
/// \brief
/// \copyright (C) Copyright Aquaveo 2018. Distributed under the xmsng
///  Software License, Version 1.0. (See accompanying file
///  LICENSE_1_0.txt or copy at http://www.aquaveo.com/xmsng/LICENSE_1_0.txt)
//------------------------------------------------------------------------------

//----- Included files ---------------------------------------------------------
#include <pybind11/pybind11.h>
#include <xmsmesh/python/meshing/meshing_py.h>

//----- Namespace declaration --------------------------------------------------
namespace py = pybind11;

//----- Python Interface -------------------------------------------------------

void initMeshing(py::module &m) {
    initMeMeshUtils(m);
    initMeMultiPolyMesher(m);
    initMePolyInput(m);
    initMeRefinePoint(m);
    initMeMultiPolyMesherIo(m);
    initMeMultiPolyTo2dm(m);
    initMePolyMesher(m);
    initMePolyRedistributePts(m);
}