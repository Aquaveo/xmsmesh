//------------------------------------------------------------------------------
/// \file
/// \brief
/// \copyright (C) Copyright Aquaveo 2018. Distributed under the xmsng
///  Software License, Version 1.0. (See accompanying file
///  LICENSE_1_0.txt or copy at http://www.aquaveo.com/xmsng/LICENSE_1_0.txt)
//------------------------------------------------------------------------------

//----- Included files ---------------------------------------------------------
#include <pybind11/pybind11.h>
#include <iostream>
#include <sstream>
#include <pybind11/numpy.h>
#include <boost/shared_ptr.hpp>

#include <xmscore/python/misc/PyUtils.h>
#include <xmsinterp/interpolate/InterpLinear.h>
#include <xmsmesh/meshing/MeMultiPolyMesherIo.h>

//----- Namespace declaration --------------------------------------------------
namespace py = pybind11;

//----- Python Interface -------------------------------------------------------
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void initMeMultiPolyMesherIo(py::module &m) {
    py::class_<xms::MeMultiPolyMesherIo, boost::shared_ptr<xms::MeMultiPolyMesherIo>>(m, "MeMultiPolyMesherIo")
        .def(py::init<>())
        .def_readwrite("check_topology", &xms::MeMultiPolyMesherIo::m_checkTopology,
            "Optional. If true, checks polygon input topology for errors."
        )
        .def_readwrite("return_cell_polygons", &xms::MeMultiPolyMesherIo::m_returnCellPolygons,
            " If true, returns the polygon index of each cell."
        )
        .def_property("points",
            [](xms::MeMultiPolyMesherIo &self) -> py::iterable {
                xms::VecPt3d &vec_pts = self.m_points;
                py::array_t<double, py::array::c_style> ret_points({(int)vec_pts.size(), 3});
                auto r = ret_points.mutable_unchecked<2>();
                int i = 0;
                for (ssize_t i = 0; i < r.shape(0); i++) {
                  r(i, 0) = vec_pts[i].x;
                  r(i, 1) = vec_pts[i].y;
                  r(i, 2) = vec_pts[i].z;
                }
                return ret_points;
            },
            [](xms::MeMultiPolyMesherIo &self, py::iterable out_poly) {
                 self.m_points.clear();
                 self.m_points.reserve(py::len(out_poly));
                 for (auto item : out_poly) {
                   if(!py::isinstance<py::iterable>(item)) {
                     throw py::type_error("First arg must be a n-tuple of 3-tuples");
                   }
                   py::tuple tuple = item.cast<py::tuple>();
                   if (py::len(tuple) != 3) {
                     throw py::type_error("Input points must be 3-tuples");
                   } else {
                     xms::Pt3d point(tuple[0].cast<double>(), tuple[1].cast<double>(), tuple[2].cast<double>());
                     self.m_points.push_back(point);
                   }
                 }
            }
        )
        .def_property("cells",
            [](xms::MeMultiPolyMesherIo &self) -> py::iterable {
                return py::array(self.m_cells.size(), self.m_cells.data());
            },
            [](xms::MeMultiPolyMesherIo &self, py::iterable cells) {
                 xms::VecInt &vecInt = self.m_cells;
                 vecInt.clear();
                 vecInt.reserve(py::len(cells));
                 for (auto item : cells) {
                    vecInt.push_back(item.cast<int>());
                 }
            }
        )
        .def_property("cell_polygons",
            [](xms::MeMultiPolyMesherIo &self) -> py::iterable {
                return py::array(self.m_cellPolygons.size(), self.m_cellPolygons.data());
            },
            [](xms::MeMultiPolyMesherIo &self, py::iterable cell_polygons) {
                 xms::VecInt &vecInt = self.m_cellPolygons;
                 vecInt.clear();
                 vecInt.reserve(py::len(cell_polygons));
                 for (auto item : cell_polygons) {
                    vecInt.push_back(item.cast<int>());
                 }
            }
        )
        .def_property("poly_inputs",
            [](xms::MeMultiPolyMesherIo &self) -> py::iterable {
                py::tuple ret_tuple(self.m_polys.size());
                for (int i = 0; i < self.m_polys.size(); i++) {
                    ret_tuple[i] = self.m_polys[i];
                }
                return ret_tuple;
            },
            [](xms::MeMultiPolyMesherIo &self, py::iterable polys) {
                 std::vector<xms::MePolyInput> &vecPolys = self.m_polys;
                 vecPolys.clear();
                 vecPolys.reserve(py::len(polys));
                 for (auto item : polys) {
                    vecPolys.push_back(item.cast<xms::MePolyInput>());
                 }
            }
        )
        .def_property("refine_points",
            [](xms::MeMultiPolyMesherIo &self) -> py::iterable {
                py::tuple ret_tuple(self.m_refPts.size());
                for (int i = 0; i < self.m_refPts.size(); i++) {
                    ret_tuple[i] = self.m_refPts[i];
                }
                return ret_tuple;
            },
            [](xms::MeMultiPolyMesherIo &self, py::iterable refine_points) {
                 std::vector<xms::MeRefinePoint> &vecRefinePoints = self.m_refPts;
                 vecRefinePoints.clear();
                 vecRefinePoints.reserve(py::len(refine_points));
                 for (auto item : refine_points) {
                    vecRefinePoints.push_back(item.cast<xms::MeRefinePoint>());
                 }
            }
        )
        ;
}

void initMePolyInput(py::module &m) {
    py::class_<xms::MePolyInput, boost::shared_ptr<xms::MePolyInput>>(m, "MePolyInput")
        .def(py::init<>([](py::iterable out_poly, py::iterable inside_polys, double bias,
                           boost::shared_ptr<xms::InterpBase> &size_function,
                           py::iterable poly_corners, boost::shared_ptr<xms::InterpBase> &elev_function) {
            xms::VecPt3d vec_out_poly;
            for (auto item : out_poly) {
              if(!py::isinstance<py::iterable>(item)) {
                throw py::type_error("First arg (out_poly) must be a n-tuple of 3-tuples");
              }
              py::tuple tuple = item.cast<py::tuple>();
              if (py::len(tuple) != 3) {
                throw py::type_error("Input points must be 3-tuples");
              } else {
                xms::Pt3d point(tuple[0].cast<double>(), tuple[1].cast<double>(), tuple[2].cast<double>());
                vec_out_poly.push_back(point);
              }
            }
            xms::VecPt3d2d vec_inside_polys;
            vec_inside_polys.reserve(py::len(inside_polys));
            for (auto poly : inside_polys) {
                if (!py::isinstance<py::iterable>(poly)) {
                    throw py::type_error("Arg (inside_polys) must be an n-tuple of @n-tuples of 3-tuples");
                }
                xms::VecPt3d vec_poly;
                vec_poly.reserve(py::len(poly));
                for (auto p : poly) {
                    if (!py::isinstance<py::iterable>(poly)) {
                        throw py::type_error("Arg (inside_polys) must be an n-tuple of @n-tuples of 3-tuples");
                    }
                    py::tuple pt = p.cast<py::tuple>();
                    if (py::len(pt) != 3) {
                        throw py::type_error("Arg (inside_polys) must be an n-tuple of n-tuples of @3-tuples");
                    } else {
                        xms::Pt3d point(pt[0].cast<double>(), pt[1].cast<double>(), pt[2].cast<double>());
                        vec_poly.push_back(point);
                    }
                }
                vec_inside_polys.push_back(vec_poly);
            }
            xms::VecInt vec_poly_corners(py::len(poly_corners));
            int k = 0;
            for (auto item : poly_corners) {
              vec_poly_corners.at(k) = item.cast<int>();
              k++;
            }
            return new xms::MePolyInput(vec_out_poly, vec_inside_polys, bias, size_function,
                                        vec_poly_corners, elev_function);
        }))
        .def(py::init<>([]() {
            return new xms::MePolyInput();
        }))
        .def_property("outside_poly",
            [](xms::MePolyInput &self) -> py::iterable {
                xms::VecPt3d &vec_pts = self.m_outPoly;
                py::array_t<double, py::array::c_style> ret_points({(int)vec_pts.size(), 3});
                auto r = ret_points.mutable_unchecked<2>();
                int i = 0;
                for (ssize_t i = 0; i < r.shape(0); i++) {
                  r(i, 0) = vec_pts[i].x;
                  r(i, 1) = vec_pts[i].y;
                  r(i, 2) = vec_pts[i].z;
                }
                return ret_points;
            },
            [](xms::MePolyInput &self, py::iterable out_poly) {
                 self.m_outPoly.clear();
                 self.m_outPoly.reserve(py::len(out_poly));
                 for (auto item : out_poly) {
                   if(!py::isinstance<py::iterable>(item)) {
                     throw py::type_error("First arg must be a n-tuple of 3-tuples");
                   }
                   py::tuple tuple = item.cast<py::tuple>();
                   if (py::len(tuple) != 3) {
                     throw py::type_error("Input points must be 3-tuples");
                   } else {
                     xms::Pt3d point(tuple[0].cast<double>(), tuple[1].cast<double>(), tuple[2].cast<double>());
                     self.m_outPoly.push_back(point);
                   }
                 }
            }
        )
        .def_property("inside_polys",
            [](xms::MePolyInput &self) -> py::iterable {
                //xms::VecPt3d2d &vec_inside_polys = self.m_insidePolys;
                py::tuple py_inside_polys(self.m_insidePolys.size());
                for (int i = 0; i < self.m_insidePolys.size(); i++) {
                    auto inside_poly = self.m_insidePolys[i];
                    py::array_t<double, py::array::c_style> poly_points({(int)inside_poly.size(), 3});
                    auto r = poly_points.mutable_unchecked<2>();
                    for (ssize_t i = 0; i < r.shape(0); i++) {
                      r(i, 0) = inside_poly[i].x;
                      r(i, 1) = inside_poly[i].y;
                      r(i, 2) = inside_poly[i].z;
                    }
                    py_inside_polys[i] = poly_points;
                }
                return py_inside_polys;
            },
            [](xms::MePolyInput &self, py::iterable inside_polys) {
                self.m_insidePolys.clear();
                self.m_insidePolys.reserve(py::len(inside_polys));
                for (auto poly : inside_polys) {
                    if (!py::isinstance<py::iterable>(poly)) {
                        throw py::type_error("Arg (inside_polys) must be an n-tuple of @n-tuples of 3-tuples");
                    }
                    xms::VecPt3d vec_poly;
                    vec_poly.reserve(py::len(poly));
                    for (auto p : poly) {
                        if (!py::isinstance<py::iterable>(poly)) {
                            throw py::type_error("Arg (inside_polys) must be an n-tuple of @n-tuples of 3-tuples");
                        }
                        py::tuple pt = p.cast<py::tuple>();
                        if (py::len(pt) != 3) {
                            throw py::type_error("Arg (inside_polys) must be an n-tuple of n-tuples of @3-tuples");
                        } else {
                            xms::Pt3d point(pt[0].cast<double>(), pt[1].cast<double>(), pt[2].cast<double>());
                            vec_poly.push_back(point);
                        }
                    }
                    self.m_insidePolys.push_back(vec_poly);
                }
            }
        )
        .def_property("poly_corners",
            [](xms::MePolyInput &self) -> py::iterable {
                return py::array(self.m_polyCorners.size(), self.m_polyCorners.data());
            },
            [](xms::MePolyInput &self, py::iterable poly_corners) {
                 self.m_polyCorners.clear();
                 for (auto item : poly_corners) {
                    self.m_polyCorners.push_back(item.cast<int>());
                 }
            }
        )
        .def_property("bound_pts_to_remove",
            [](xms::MePolyInput &self) -> py::iterable {
                xms::VecPt3d &vec_pts = self.m_boundPtsToRemove;
                py::array_t<double, py::array::c_style> ret_points({(int)vec_pts.size(), 3});
                auto r = ret_points.mutable_unchecked<2>();
                int i = 0;
                for (ssize_t i = 0; i < r.shape(0); i++) {
                  r(i, 0) = vec_pts[i].x;
                  r(i, 1) = vec_pts[i].y;
                  r(i, 2) = vec_pts[i].z;
                }
                return ret_points;
            },
            [](xms::MePolyInput &self, py::iterable out_poly) {
                 self.m_boundPtsToRemove.clear();
                 self.m_boundPtsToRemove.reserve(py::len(out_poly));
                 for (auto item : out_poly) {
                   if(!py::isinstance<py::iterable>(item)) {
                     throw py::type_error("First arg must be a n-tuple of 3-tuples");
                   }
                   py::tuple tuple = item.cast<py::tuple>();
                   if (py::len(tuple) != 3) {
                     throw py::type_error("Input points must be 3-tuples");
                   } else {
                     xms::Pt3d point(tuple[0].cast<double>(), tuple[1].cast<double>(), tuple[2].cast<double>());
                     self.m_boundPtsToRemove.push_back(point);
                   }
                 }
            }
        )
        .def_readwrite("bias", &xms::MePolyInput::m_bias," Optional. Factor for transitioning between areas of"
            " high refinement to less refinement."
        )
        .def_readwrite("size_function", &xms::MePolyInput::m_sizeFunction,
            "Optional. Size function for scalar paving."
        )
        .def_readwrite("elev_function", &xms::MePolyInput::m_elevFunction,
            "Optional. Elevation function for interpolating z coordinate of mesh points."
        )
        .def_readwrite("const_size_function", &xms::MePolyInput::m_constSizeFunction,
            "Optional. Transition factor for constant size function."
        )
        .def_readwrite("const_size_bias", &xms::MePolyInput::m_constSizeBias,
            "Optional. Transition factor for constant size function."
        )
        .def_readwrite("remove_internal_four_triangle_pts", &xms::MePolyInput::m_removeInternalFourTrianglePts,
            "Optional. Remove internal points that are only connected to 4 cells. Used by the ugAutoCorrectCells class"
        )
        .def_readwrite("poly_id", &xms::MePolyInput::m_polyId,
            "Optional. Set when needed. Can be useful for classes who need an ID."
        )
        .def_property("seed_points",
            [](xms::MePolyInput &self) -> py::iterable {
                return xms::PyIterFromVecPt3d(self.m_seedPts);
            },
            [](xms::MePolyInput &self, py::iterable seed_points) {
                self.m_seedPts = *xms::VecPt3dFromPyIter(seed_points);
            }, "Optional array of seed points. If the user has some methodology for creating points inside the polygon"
            " then those points can be specified here. If these points are specified then the paving is not performed."
            " These points will not be used if the meshing option is patch."
        )
        .def("__str__", [](xms::MePolyInput &self) {
             std::string szf = self.m_sizeFunction == nullptr ? "none" : self.m_sizeFunction->ToString();
             std::string elevf = self.m_elevFunction == nullptr ? "none" : self.m_elevFunction->ToString();
             std::stringstream ss;
             ss <<  "bias: " << self.m_bias << std::endl <<
                    "size_func: " << szf << std::endl <<
                    "elev_func: " << elevf << std::endl <<
                    "const_size_function: " << self.m_constSizeFunction << std::endl <<
                    "const_size_bias: " << self.m_constSizeBias << std::endl <<
                    "outside_poly size: " << self.m_outPoly.size() << std::endl <<
                    "inside_polys size: " << self.m_insidePolys.size() << std::endl;
             return ss.str();
        },"outputs contents as string")
    ;
}

void initMeRefinePoint(py::module &m) {
    py::class_<xms::MeRefinePoint, boost::shared_ptr<xms::MeRefinePoint>>(m, "MeRefinePoint")
        .def(py::init<>([](py::tuple pt, double size, bool create_mesh_point) {
            if(py::len(pt) != 3) {
                throw py::type_error("Input point should be a 3-tuple");
            }
            xms::Pt3d point(pt[0].cast<double>(), pt[1].cast<double>(), pt[2].cast<double>());
            return new xms::MeRefinePoint(point, size, create_mesh_point);
        }))
        .def_property("point",
          [](xms::MeRefinePoint &self) -> py::tuple {
            return py::make_tuple(self.m_pt.x, self.m_pt.y, self.m_pt.z);
          },
          [](xms::MeRefinePoint &self, py::tuple pt) {
            if(py::len(pt) != 3) {
              throw py::type_error("Input point should be a 3-tuple");
            } else {
              xms::Pt3d point(pt[0].cast<double>(), pt[1].cast<double>(), pt[2].cast<double>());
              self.m_pt = point;
            }
          }
         )
        .def_readwrite("size", &xms::MeRefinePoint::m_size,"  /// Element size at the refine point."
            " A negative value indicates a hard point.")
        .def_readwrite("create_mesh_point", &xms::MeRefinePoint::m_createMeshPoint)
        .def("__str__", [](xms::MeRefinePoint &self) {
            std::stringstream ss;
            ss << "point: (" << self.m_pt << ")" <<  std::endl <<
                   "size: " << self.m_size <<  std::endl <<
                   "create_mesh_point: " << self.m_createMeshPoint << std::endl;
            return ss.str();
        }, "Returns contents as string")
        ;
}
