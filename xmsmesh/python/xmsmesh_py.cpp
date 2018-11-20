//------------------------------------------------------------------------------
/// \file
/// \brief root module for xmsmesh Python bindings.
/// \copyright (C) Copyright Aquaveo 2018. Distributed under FreeBSD License
/// (See accompanying file LICENSE or https://aqaveo.com/bsd/license.txt)
//------------------------------------------------------------------------------

//----- Included files ---------------------------------------------------------
#include <pybind11/pybind11.h>
#include <xmsmesh/python/meshing/meshing_py.h>

//----- Namespace declaration --------------------------------------------------
namespace py = pybind11;

//----- Python Interface -------------------------------------------------------
//PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

std::string version() {
    return "1.0.0";
}


//------ Primary Module --------------------------------------------------------
PYBIND11_MODULE(xmsmesh_py, m) {
    m.doc() = "Python bindings for xmsmesh"; // optional module docstring
    m.def("version", &version,
          "Get current version of xmsmesh Python bindings.");

    // Interpolate module
    py::module modMeshing = m.def_submodule("meshing");
    initMeshing(modMeshing);
}