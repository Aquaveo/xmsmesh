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
#include <xmscore/python/misc/PublicObserver.h>
#include <xmsmesh/meshing/MeMultiPolyMesherIo.h>
#include <xmsmesh/meshing/MePolyMesher.h>

//----- Namespace declaration --------------------------------------------------
namespace py = pybind11;

//----- Python Interface -------------------------------------------------------
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void initMePolyMesher(py::module &m) {
    py::class_<xms::MePolyMesher, boost::shared_ptr<xms::MePolyMesher>>(m, "MePolyMesher")
        .def(py::init(&xms::MePolyMesher::New))
        .def("mesh_it", [](xms::MePolyMesher &self,
                           xms::MeMultiPolyMesherIo &mesh_io,
                           size_t polyIdx) -> py::tuple {
            xms::VecPt3d vec_pts;
            xms::VecInt vec_tris, vec_cell;
            bool result = self.MeshIt(mesh_io, polyIdx, vec_pts, vec_tris, vec_cell);

            py::array_t<double, py::array::c_style> ret_points({(int)vec_pts.size(), 3});
            auto r = ret_points.mutable_unchecked<2>();
            int i = 0;
            for (ssize_t i = 0; i < r.shape(0); i++) {
              r(i, 0) = vec_pts[i].x;
              r(i, 1) = vec_pts[i].y;
              r(i, 2) = vec_pts[i].z;
            }
            py::array ret_tris(vec_tris.size(), vec_tris.data());
            py::array ret_cell(vec_cell.size(), vec_cell.data());
            return py::make_tuple(result, ret_points, ret_tris, ret_cell);
        })
        .def("set_observer", [](xms::MePolyMesher &self,
                                boost::shared_ptr<xms::PublicObserver> obs) {
            self.SetObserver(obs);
        })
        .def("get_processed_refine_pts", [](xms::MePolyMesher &self) -> py::iterable {
            xms::VecPt3d vec_pts;

            self.GetProcessedRefinePts(vec_pts);

            py::array_t<double, py::array::c_style> ret_points({(int)vec_pts.size(), 3});
            auto r = ret_points.mutable_unchecked<2>();
            int i = 0;
            for (ssize_t i = 0; i < r.shape(0); i++) {
              r(i, 0) = vec_pts[i].x;
              r(i, 1) = vec_pts[i].y;
              r(i, 2) = vec_pts[i].z;
            }

            return ret_points;
        })
        ;
}