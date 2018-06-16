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
#include <xmscore/misc/DynBitset.h>
#include <xmsinterp/triangulate/TrTin.h>
#include <xmsmesh/meshing/MeMeshUtils.h>

//----- Namespace declaration --------------------------------------------------
namespace py = pybind11;

//----- Python Interface -------------------------------------------------------
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void initMeshUtils(py::module &m) {

    py::module modMeshUtils = m.def_submodule("MeshUtils")
    .def("size_function_from_depth", [](py::iterable depths, double min_size,
                                         double max_size) -> py::iterable {

        xms::VecDbl vec_depths, vec_size;
        for (auto item : depths) {
          vec_depths.push_back(item.cast<double>());
        }
        xms::meSizeFunctionFromDepth(vec_depths, vec_size, min_size, max_size);

        if (py::isinstance<py::array>(depths)) {
          // NOTE: This is a copy operation
          return py::array(vec_size.size(), vec_size.data());
        } else {
          // NOTE: This is a copy operation
          auto tuple_ret = py::tuple(vec_size.size());
          for (size_t i = 0; i < vec_size.size(); ++i) {
            tuple_ret[i] = vec_size.at(i);
          }
          return tuple_ret;
        }
    })
    .def("smooth_size_function", [](boost::shared_ptr<xms::TrTin> tin, py::iterable sizes,
                                    double size_ratio, double min_size, int anchor_type,
                                    py::iterable pts_flag) -> py::iterable {
        xms::VecFlt vec_sizes, vec_smooth_sizes;
        for (auto item : sizes) {
          vec_sizes.push_back(item.cast<float>());
        }
        xms::DynBitset bitset;
        std::vector<unsigned char> bitvals;
        for (auto item : pts_flag) {
          py::bool_ flag = item.cast<py::bool_>();
          if (flag) {
            bitvals.push_back(1);
          } else {
            bitvals.push_back(0);
          }
        }
        xms::VecBooleanToDynBitset(bitvals, bitset);

        xms::meSmoothSizeFunction(tin, vec_sizes, size_ratio, min_size, anchor_type, bitset, vec_smooth_sizes);

        if (py::isinstance<py::array>(sizes)) {
          // NOTE: This is a copy operation
          return py::array(vec_smooth_sizes.size(), vec_smooth_sizes.data());
        } else {
          // NOTE: This is a copy operation
          auto tuple_ret = py::tuple(vec_smooth_sizes.size());
          for (size_t i = 0; i < vec_smooth_sizes.size(); ++i) {
            tuple_ret[i] = vec_smooth_sizes.at(i);
          }
          return tuple_ret;
        }
    })
    .def("smooth_elev_by_slope", [](boost::shared_ptr<xms::TrTin> tin, py::iterable sizes,
                                    double max_size, int anchor_type,
                                    py::iterable pts_flag) -> py::iterable {
        xms::VecFlt vec_sizes, vec_smooth_sizes;
        for (auto item : sizes) {
          vec_sizes.push_back(item.cast<float>());
        }
        xms::DynBitset bitset;
        std::vector<unsigned char> bitvals;
        for (auto item : pts_flag) {
          py::bool_ flag = item.cast<py::bool_>();
          if (flag) {
            bitvals.push_back(1);
          } else {
            bitvals.push_back(0);
          }
        }
        xms::VecBooleanToDynBitset(bitvals, bitset);

        xms::meSmoothElevBySlope(tin, vec_sizes, max_size, anchor_type, bitset, vec_smooth_sizes);

        if (py::isinstance<py::array>(sizes)) {
          // NOTE: This is a copy operation
          return py::array(vec_smooth_sizes.size(), vec_smooth_sizes.data());
        } else {
          // NOTE: This is a copy operation
          auto tuple_ret = py::tuple(vec_smooth_sizes.size());
          for (size_t i = 0; i < vec_smooth_sizes.size(); ++i) {
            tuple_ret[i] = vec_smooth_sizes.at(i);
          }
          return tuple_ret;
        }
    })
    ;
}