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
#include <xmsinterp/interpolate/InterpBase.h>
#include <xmsmesh/meshing/MePolyRedistributePts.h>

//----- Namespace declaration --------------------------------------------------
namespace py = pybind11;

//----- Python Interface -------------------------------------------------------
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void initMePolyRedistributePts(py::module &m) {
    py::class_<xms::MePolyRedistributePts,
                boost::shared_ptr<xms::MePolyRedistributePts>> polyRedistribute(m, "PolyRedistributePts");


    polyRedistribute.def(py::init(&xms::MePolyRedistributePts::New));
    // -------------------------------------------------------------------------
    // function: set_size_func
    // -------------------------------------------------------------------------
    const char* set_size_func_doc = R"pydoc(
        Sets the size function interpolator

        Args:
            interp (InterpBase): Size function interpolator class
    )pydoc";
    polyRedistribute.def("set_size_func", [](xms::MePolyRedistributePts &self,
        boost::shared_ptr<xms::InterpBase> interp)
        {
            self.SetSizeFunc(interp);
        },set_size_func_doc,py::arg("interp"));
    // -------------------------------------------------------------------------
    // function: __repr__
    // -------------------------------------------------------------------------
    polyRedistribute.def("__repr__", [](xms::MePolyRedistributePts &self) {
      return self.ToPyRepr();
    });
    // -------------------------------------------------------------------------
    // function: set_size_func_from_poly
    // -------------------------------------------------------------------------
    const char* set_size_func_from_poly_doc = R"pydoc(
        Creates an interpolator that uses the spacing on the input polygon as 
        its scalar

        Args:
            out_poly (iterable): The outside polygon
            inside_polys (iterable): Inside polygons that are inside of 
              a_outPoly
            size_bias (float): A factor used in transitioning the size
    )pydoc";
    polyRedistribute.def("set_size_func_from_poly", [](xms::MePolyRedistributePts &self,
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
        },set_size_func_from_poly_doc, py::arg("out_poly"),
          py::arg("inside_polys"),py::arg("size_bias"));
    // -------------------------------------------------------------------------
    // function: set_constant_size_func
    // -------------------------------------------------------------------------
    const char* set_constant_size_func_doc = R"pydoc(
        Sets the size function to a constant value

        Args:
            size (float): The element edge size.
    )pydoc";
    polyRedistribute.def("set_constant_size_func", &xms::MePolyRedistributePts::SetConstantSizeFunc,
          set_constant_size_func_doc,py::arg("size"));
    // -------------------------------------------------------------------------
    // function: set_constant_size_bias
    // -------------------------------------------------------------------------
    const char* set_constant_size_bias_doc = R"pydoc(
        Sets the bias for constant value size function

        Args:
            size_bias (float): Transition rate for size function
    )pydoc";
    polyRedistribute.def("set_constant_size_bias", &xms::MePolyRedistributePts::SetConstantSizeBias,
          set_constant_size_bias_doc,py::arg("size_bias"));
    // -------------------------------------------------------------------------
    // function: set_use_curvature_redistribution
    // -------------------------------------------------------------------------
    const char* set_use_curvature_redistribution_doc = R"pydoc(
        Specifies that curvature redistribution will be used.

        Args:
            feature_size (float): The size of the smallest feature in the 
              polyline to be detected.
            mean_spacing (float): The mean spacing between the distributed 
              points.
            minimum_curvature (float): The value of the curvature to be used 
              instead of 0 in staight lines. It limits the maximum spacing 
              between points. If not included, the default is 0.001.
            smooth (bool): Detemines if the curvatures are to be averaged by a 
              rolling 0.25-0.5-0.25 weighted rolling average.
    )pydoc";
    polyRedistribute.def("set_use_curvature_redistribution", &xms::MePolyRedistributePts::SetUseCurvatureRedistribution,
          set_use_curvature_redistribution_doc,py::arg("feature_size"),
          py::arg("mean_spacing"),py::arg("minimum_curvature"),py::arg("smooth")
        );
    // -------------------------------------------------------------------------
    // function: redistribute
    // -------------------------------------------------------------------------
    const char* redistribute_doc = R"pydoc(
        Redistributes points on closed loop polylines. The length of edges in 
        the redistribution comes from a size function that is interpolated to 
        the points that make up the polylines. By default this size function 
        comes from the edge lengths in the original polygon.

        Args:
            poly_line (iterable): Input closed loop polylines

        Returns:
          iterable: Redistributed closed loop polylines, Number of iterations 
            from the polygon boundary
    )pydoc";
    polyRedistribute.def("redistribute", [](xms::MePolyRedistributePts &self,
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
        },redistribute_doc,py::arg("poly_line"));
}