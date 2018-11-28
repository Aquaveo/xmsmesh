//------------------------------------------------------------------------------
/// \file
/// \brief
/// \copyright (C) Copyright Aquaveo 2018. Distributed under FreeBSD License
/// (See accompanying file LICENSE or https://aqaveo.com/bsd/license.txt)
//------------------------------------------------------------------------------

//----- Included files ---------------------------------------------------------
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <sstream>
#include <boost/shared_ptr.hpp>
#include <xmsmesh/meshing/MeMultiPolyTo2dm.h>
#include <xmsmesh/meshing/MeMultiPolyMesherIo.h>

//----- Namespace declaration --------------------------------------------------
namespace py = pybind11;

//----- Python Interface -------------------------------------------------------
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void initMeMultiPolyTo2dm(py::module &m) {
    py::class_<xms::MeMultiPolyTo2dm, boost::shared_ptr<xms::MeMultiPolyTo2dm>> poly2dm(m, "MultiPolyTo2dm");
    poly2dm.def(py::init(&xms::MeMultiPolyTo2dm::New));
    // -------------------------------------------------------------------------
    // function: generate_2dm
    // -------------------------------------------------------------------------
    const char* generate_2dm_doc = R"pydoc(
        Creates a 2dm file from polygons

        Args:
          mesh_io (MeMultiPolyMesherIo): Input/output of polygons and options 
            for generating a mesh.
          fname (str): output filename
          precision (:obj:`int`, optional) decimal point precision
        
        Returns:
          tuple: true if the mesh was generated., and resultant filename
    )pydoc";
    poly2dm.def("generate_2dm", [](xms::MeMultiPolyTo2dm &self,
                           xms::MeMultiPolyMesherIo &mesh_io,
                           std::string fname, int precision) -> py::tuple {
          if (fname.empty()) {
            std::stringstream ss;
            bool result = self.Generate2dm(mesh_io, ss, precision);
            return py::make_tuple(result, ss.str());

          }
          return py::make_tuple(self.Generate2dm(mesh_io, fname), "");
        },generate_2dm_doc,py::arg("mesh_io"),py::arg("fname"),py::arg("precision")=15);
}