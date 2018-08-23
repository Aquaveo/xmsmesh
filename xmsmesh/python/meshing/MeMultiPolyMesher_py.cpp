//------------------------------------------------------------------------------
/// \file
/// \brief
/// \copyright (C) Copyright Aquaveo 2018. Distributed under the xmsng
///  Software License, Version 1.0. (See accompanying file
///  LICENSE_1_0.txt or copy at http://www.aquaveo.com/xmsng/LICENSE_1_0.txt)
//------------------------------------------------------------------------------

//----- Included files ---------------------------------------------------------
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <boost/shared_ptr.hpp>
#include <xmscore/misc/XmError.h>
#include <xmscore/python/misc/PublicObserver.h>
#include <xmsmesh/meshing/MeMultiPolyMesherIo.h>
#include <xmsmesh/meshing/MeMultiPolyMesher.h>

//----- Namespace declaration --------------------------------------------------
namespace py = pybind11;

//----- Python Interface -------------------------------------------------------
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void initMeMultiPolyMesher(py::module &m) {
    py::class_<xms::MeMultiPolyMesher, boost::shared_ptr<xms::MeMultiPolyMesher>>(m, "MeMultiPolyMesher")
        .def(py::init(&xms::MeMultiPolyMesher::New))
        .def("mesh_it", [](xms::MeMultiPolyMesher &self,
                           xms::MeMultiPolyMesherIo &mesh_io) -> py::iterable {
            if (self.MeshIt(mesh_io)) {
              return py::make_tuple(true, "");
            } else {
              std::string errors = xms::XmLog::Instance().GetAndClearStackStr();
              return py::make_tuple(false, errors);
            }
        },"Creates a triangle mesh from the input polygons. The polygons can not overlap.",
          py::arg("mesh_io"))
        .def("set_observer", [](xms::MeMultiPolyMesher &self,
                                boost::shared_ptr<xms::PublicObserver> obs) {
            self.SetObserver(obs);
        },"Set the observer to use for feedback while processing.",
          py::arg("obs")
        )
        ;
}