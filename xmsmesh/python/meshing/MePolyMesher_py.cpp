//------------------------------------------------------------------------------
/// \file
/// \brief
/// \copyright (C) Copyright Aquaveo 2018. Distributed under FreeBSD License
/// (See accompanying file LICENSE or https://aqaveo.com/bsd/license.txt)
//------------------------------------------------------------------------------

//----- Included files ---------------------------------------------------------
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <boost/shared_ptr.hpp>
#include <xmsmesh/meshing/MeMultiPolyMesherIo.h>
#include <xmsmesh/meshing/MePolyMesher.h>

//----- Namespace declaration --------------------------------------------------
namespace py = pybind11;

//----- Python Interface -------------------------------------------------------
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void initMePolyMesher(py::module &m) {
    py::class_<xms::MePolyMesher, boost::shared_ptr<xms::MePolyMesher>> polyMesher(m, "PolyMesher");


    polyMesher.def(py::init(&xms::MePolyMesher::New));
    // -------------------------------------------------------------------------
    // function: mesh_it
    // -------------------------------------------------------------------------
    const char* mesh_it_doc = R"pydoc(
        Perform MESH_PAVE, MESH_SPAVE, MESH_PATCH meshing on a polygon.

        Args:
          mesh_io (MeMultiPolyMesherIo):Meshing input: polygons and optional
            inputs
          poly_idx (int): Index to the polygon in a_input to mesh
        
        Returns:
          tuple: True if no errors encountered, computed mesh points, computed 
            mesh triangles from paving, computed mesh cells from patch.
    )pydoc";
    polyMesher.def("mesh_it", [](xms::MePolyMesher &self,
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
        },mesh_it_doc,py::arg("mesh_io"),py::arg("poly_idx"));
    // -------------------------------------------------------------------------
    // function: get_processed_refine_pts
    // -------------------------------------------------------------------------
    const char* get_processed_refine_pts_doc = R"pydoc(
        Gets the refine points that were inside the polygon, both points that 
        are included in the meshing process and those that were not.

        Returns:
            iterable: Locations of refine points used inside of this polygon.
    )pydoc";
    polyMesher.def("get_processed_refine_pts", [](xms::MePolyMesher &self) -> py::iterable {
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
        },get_processed_refine_pts_doc);
}