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
#include <sstream>
#include <boost/shared_ptr.hpp>
#include <xmscore/python/misc/PublicObserver.h>
#include <xmsmesh/meshing/MeMultiPolyTo2dm.h>
#include <xmsmesh/meshing/MeMultiPolyMesherIo.h>

//----- Namespace declaration --------------------------------------------------
namespace py = pybind11;

//----- Python Interface -------------------------------------------------------
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void initMeMultiPolyTo2dm(py::module &m) {
    py::class_<xms::MeMultiPolyTo2dm, boost::shared_ptr<xms::MeMultiPolyTo2dm>>(m, "MeMultiPolyTo2dm")
        .def(py::init(&xms::MeMultiPolyTo2dm::New))
        .def("generate_2dm", [](xms::MeMultiPolyTo2dm &self,
                           xms::MeMultiPolyMesherIo &mesh_io,
                           std::string fname) -> py::tuple {
          if (fname.empty()) {
            std::stringstream ss;
            bool result = self.Generate2dm(mesh_io, ss);
            return py::make_tuple(result, ss.str());

          }
          return py::make_tuple(self.Generate2dm(mesh_io, fname), "");
        })
        .def("set_observer", [](xms::MeMultiPolyTo2dm &self,
                                boost::shared_ptr<xms::PublicObserver> obs) {
            self.SetObserver(obs);
        })
        ;
}