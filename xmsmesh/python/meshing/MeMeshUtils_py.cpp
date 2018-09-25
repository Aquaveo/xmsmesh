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

void initMeMeshUtils(py::module &m) {

    py::module modMeshUtils = m.def_submodule("MeMeshUtils");


  // ---------------------------------------------------------------------------
  // function: size_function_from_depth
  // ---------------------------------------------------------------------------
  const char* size_function_from_depth_doc = R"pydoc(
      Creates a size at each point based on the depth at the point and the min 
      and max sizes the equation is  min_depth + ( (depth - min_depth) / 
      (max_depth - min_depth) ) * (max_size - min_size)

      Args:
          depths (iterable): The measured depths at point locations
          min_size (float): The minimum element edge size
          max_size (float): The maximum element edge size
  )pydoc";
    modMeshUtils.def("size_function_from_depth", [](py::iterable depths, double min_size,
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
    },size_function_from_depth_doc,py::arg("depths"),py::arg("min_size"),
    py::arg("max_size"));
  // ---------------------------------------------------------------------------
  // function: smooth_size_function
  // ---------------------------------------------------------------------------
  const char* smooth_size_function_doc = R"pydoc(
      Smooths a size function. Ensures that the size function transitions over a
      sufficient distance so that the area change of adjacent elements meets the
      size ratio passed in.

      Args:
          tin (boost::shared_ptr<xms::TrTin>): Points and triangles defining the
            connectivity of the size function.
          sizes (iterable): Array of the current sizes
          size_ratio (float): Allowable size difference between adjacent 
            elements
          min_size (float): Minimum specified element size
          anchor_type (int): The minimum element edge size
          pts_flag (iterable): Flag to indicate if the value at the point should
            be adjusted (a value of true will skip the point). Leave the bitset 
            empty to process all points.

      Returns:
        iterable: Array of smoothed sizes
  )pydoc";
    modMeshUtils.def("smooth_size_function", [](boost::shared_ptr<xms::TrTin> tin, py::iterable sizes,
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
    },smooth_size_function_doc, py::arg("tin"),py::arg("sizes"),
      py::arg("size_ratio"),py::arg("min_size"), py::arg("anchor_type"),
      py::arg("pts_flag")  
    );
  // ---------------------------------------------------------------------------
  // function: smooth_elev_by_slope
  // ---------------------------------------------------------------------------
  const char* smooth_elev_by_slope_doc = R"pydoc(
      Smooths a elevations based on max specified slope (a_maxSlope) preserving 
      either the min or max based on a_anchorType

      Args:
          tin (boost::shared_ptr<xms::TrTin>): Points and triangles defining the
            connectivity of the size function.
          sizes (iterable): Array of the current sizes
          max_slope (float): Maximum allowable slope
          anchor_type (int): The minimum element edge size
          pts_flag (iterable): Flag to indicate if the value at the point should
            be adjusted (a value of true will skip the point). Leave the bitset 
            empty to process all points.

      Returns:
        iterable: Array of smoothed elevations
  )pydoc";
    modMeshUtils.def("smooth_elev_by_slope", [](boost::shared_ptr<xms::TrTin> tin, py::iterable sizes,
                                    double max_slope, int anchor_type,
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

        xms::meSmoothElevBySlope(tin, vec_sizes, max_slope, anchor_type, bitset, vec_smooth_sizes);

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
    },smooth_elev_by_slope_doc, py::arg("tin"),py::arg("sizes"),
    py::arg("max_slope"),py::arg("anchor_type"),py::arg("pts_flag"));
}