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
#include <xmscore/misc/DynBitset.h>
#include <xmscore/misc/XmError.h>
#include <xmscore/python/misc/PyUtils.h>
#include <xmsinterp/triangulate/TrTin.h>
#include <xmsmesh/meshing/MeMeshUtils.h>
#include <xmsmesh/meshing/MeMultiPolyMesher.h>
#include <xmsmesh/meshing/MeMultiPolyMesherIo.h>
#include <xmsmesh/meshing/MeMultiPolyTo2dm.h>


//----- Namespace declaration --------------------------------------------------
namespace py = pybind11;

//----- Python Interface -------------------------------------------------------
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void initMeMeshUtils(py::module &m) {

    py::module modMeshUtils = m.def_submodule("mesh_utils");


  // ---------------------------------------------------------------------------
  // function: size_function_from_depth
  // ---------------------------------------------------------------------------
  const char* size_function_from_depth_doc = R"pydoc(
      Creates a size at each point based on the depth at the point and the min 
      and max sizes the equation is  min_depth + ( (depth - min_depth) / 
      (max_depth - min_depth) ) * (max_size - min_size). This is often useful for
      coastal numerical model simulations.

      Args:
          depths (iterable): The measured depths at point locations
          min_size (float): The minimum element edge size
          max_size (float): The maximum element edge size

      Returns:
        iterable: Array of sizes based on depth
  )pydoc";
    modMeshUtils.def("size_function_from_depth", [](py::iterable depths, double min_size,
                                         double max_size) -> py::iterable {

        xms::VecDbl vec_depths, vec_size;
        vec_depths = *xms::VecDblFromPyIter(depths);
        xms::meSizeFunctionFromDepth(vec_depths, vec_size, min_size, max_size);
        return xms::PyIterFromVecDbl(vec_size);
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
          tin (:class:`Tin <xmsinterp.triangulate.Tin>`): Points and triangles defining the connectivity of the size function.
          sizes (iterable): Array of the current sizes
          size_ratio (float): Allowable size difference between adjacent elements
          min_size (float): Minimum specified element size
          anchor_type (int): The minimum element edge size
          pts_flag (iterable): Flag to indicate if the value at the point should be adjusted (a value of true will skip the point). Leave the bitset empty to process all points.

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
      Smooths a elevations based on max specified slope (max_slope) preserving
      either the min or max based on anchor_type

      Args:
          tin (:class:`Tin <xmsinterp.triangulate.Tin>`): Points and triangles defining the connectivity of the elevations.
          elevations (iterable): Array of the current elevations
          max_slope (float): Maximum allowable slope
          anchor_type (int): Indicates weather you are anchoring to the top or bottom of the slope.
          pts_flag (iterable): Flag to indicate if the value at the point should be adjusted (a value of true will skip the point). Leave the bitset empty to process all points.

      Returns:
        iterable: Array of smoothed elevations
  )pydoc";
    modMeshUtils.def("smooth_elev_by_slope", [](boost::shared_ptr<xms::TrTin> tin, py::iterable elevations,
                                    double max_slope, int anchor_type,
                                    py::iterable pts_flag) -> py::iterable {
        xms::VecFlt vec_elevations, vec_smooth_elevations;
        for (auto item : elevations) {
          vec_elevations.push_back(item.cast<float>());
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

        xms::meSmoothElevBySlope(tin, vec_elevations, max_slope, anchor_type, bitset, vec_smooth_elevations);

        if (py::isinstance<py::array>(elevations)) {
          // NOTE: This is a copy operation
          return py::array(vec_smooth_elevations.size(), vec_smooth_elevations.data());
        } else {
          // NOTE: This is a copy operation
          auto tuple_ret = py::tuple(vec_smooth_elevations.size());
          for (size_t i = 0; i < vec_smooth_elevations.size(); ++i) {
            tuple_ret[i] = vec_smooth_elevations.at(i);
          }
          return tuple_ret;
        }
    },smooth_elev_by_slope_doc, py::arg("tin"),py::arg("elevations"),
    py::arg("max_slope"),py::arg("anchor_type"),py::arg("pts_flag"));

  // ---------------------------------------------------------------------------
  // function: generate_mesh
  // ---------------------------------------------------------------------------
  const char* generate_mesh_doc = R"pydoc(
      Creates a mesh from the input polygons.

      Args:
          mesh_io (:class:`MultiPolyMesherIo <xmsmesh.meshing.MultiPolyMesherIo>`): Input polygons and options for generating a mesh.

      Returns:
        tuple: true if the mesh was generated successfully false otherwise, and a string of messages.
  )pydoc";
    modMeshUtils.def("generate_mesh",
     [](xms::MeMultiPolyMesherIo &mesh_io) -> py::iterable
     {
       BSHP<xms::MeMultiPolyMesher> multiPolyMesher = xms::MeMultiPolyMesher::New();
       bool rval = multiPolyMesher->MeshIt(mesh_io);
       std::string errors = xms::XmLog::Instance().GetAndClearStackStr();
       return py::make_tuple(rval, errors);
     },generate_mesh_doc, py::arg("mesh_io"));
  // ---------------------------------------------------------------------------
  // function: generate_2dm
  // ---------------------------------------------------------------------------
    const char* generate_2dm_doc = R"pydoc(
        Creates a mesh from the input polygons and writes it to a 2dm file.

        Args:
            mesh_io (:class:`MultiPolyMesherIo <xmsmesh.meshing.MultiPolyMesherIo>`): Input polygons and options for generating a mesh.
            file_name (str): The file name of the output 2dm file.
            precision (int, optional): The decimal point precision of the resulting mesh.

        Returns:
            tuple: true if the mesh was generated successfully false otherwise, and a string of messages.
    )pydoc";
    modMeshUtils.def("generate_2dm",
     [](xms::MeMultiPolyMesherIo &mesh_io,
        std::string file_name, int precision) -> py::tuple {
        BSHP<xms::MeMultiPolyTo2dm> mesher = xms::MeMultiPolyTo2dm::New();
        if (file_name.empty()) {
          throw py::value_error("file_name not specifed. Aborting mesh procedure.");
        }
        bool result = mesher->Generate2dm(mesh_io, file_name, precision);
        std::string errors = xms::XmLog::Instance().GetAndClearStackStr();
        return py::make_tuple(result, errors);
        },generate_2dm_doc,py::arg("mesh_io"),py::arg("file_name"),py::arg("precision")=15);

}