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
#include <xmsinterp/interpolate/InterpBase.h>
#include <xmsmesh/meshing/MePolyRedistributePts.h>

//----- Namespace declaration --------------------------------------------------
namespace py = pybind11;

//----- Python Interface -------------------------------------------------------
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void initMePolyRedistributePts(py::module &m) {
    py::class_<xms::MePolyRedistributePts, boost::shared_ptr<xms::MePolyRedistributePts>>(m, "MePolyRedistributePts")
        .def(py::init(&xms::MePolyRedistributePts::New))
        .def("set_size_func", [](xms::MePolyRedistributePts &self,
        boost::shared_ptr<xms::InterpBase> interp)
        {
            self.SetSizeFunc(interp);
        })
        .def("set_size_func_from_poly", [](xms::MePolyRedistributePts &self,
                                           py::iterable out_poly,
                                           py::iterable inside_polys,
                                           double size_bias) {

            xms::VecPt3d vec_out_poly;
            for (auto item : out_poly) {
              if(!py::isinstance<py::iterable>(item)) {
                throw py::type_error("First arg must be a n-tuple of 3-tuples");
              }
              py::tuple tuple = item.cast<py::tuple>();
              if (py::len(tuple) != 3) {
                throw py::type_error("Input points must be 3-tuples");
              } else {
                xms::Pt3d point(tuple[0].cast<double>(), tuple[1].cast<double>(), tuple[2].cast<double>());
                vec_out_poly.push_back(point);
              }
            }
            xms::VecPt3d2d vec_inside_polys(py::len(inside_polys));
            int i = 0;
            for (auto poly : inside_polys) {
                if (!py::isinstance<py::iterable>(poly)) {
                    throw py::type_error("Second arg must be an n-tuple of n-tuples of 3-tuples");
                }
                xms::VecPt3d vec_poly(py::len(poly));
                int j = 0;
                for (auto p : poly) {
                    if (!py::isinstance<py::iterable>(poly)) {
                        throw py::type_error("Second arg must be an n-tuple of n-tuples of 3-tuples");
                    }
                    py::tuple pt = p.cast<py::tuple>();
                    if (py::len(pt) != 3) {
                        throw py::type_error("Second arg must be an n-tuple of n-tuples of 3-tuples");
                    } else {
                        xms::Pt3d point(pt[0].cast<double>(), pt[1].cast<double>(), pt[2].cast<double>());
                        vec_poly.at(j) = point;
                    }
                }
                vec_inside_polys.at(i) = vec_poly;
            }
            self.SetSizeFuncFromPoly(vec_out_poly, vec_inside_polys, size_bias);
        })
        .def("set_constant_size_func", &xms::MePolyRedistributePts::SetConstantSizeFunc)
        .def("set_constant_size_bias", &xms::MePolyRedistributePts::SetConstantSizeBias)
        .def("redistribute", [](xms::MePolyRedistributePts &self,
                                py::iterable poly_line) -> py::iterable {
           xms::VecPt3d vec_poly_line;
            for (auto item : poly_line) {
              if(!py::isinstance<py::iterable>(item)) {
                throw py::type_error("First arg must be a n-tuple of 3-tuples");
              }
              py::tuple tuple = item.cast<py::tuple>();
              if (py::len(tuple) != 3) {
                throw py::type_error("Input points must be 3-tuples");
              } else {
                xms::Pt3d point(tuple[0].cast<double>(), tuple[1].cast<double>(), tuple[2].cast<double>());
                vec_poly_line.push_back(point);
              }
            }

            xms::VecPt3d vec_pts(self.Redistribute(vec_poly_line));

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