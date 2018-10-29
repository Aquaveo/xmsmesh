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
#include <xmscore/python/misc/PyUtils.h>
#include <xmsgrid/ugrid/XmUGrid.h>
#include <xmsmesh/meshing/detail/MeQuadBlossom.h>

//----- Namespace declaration --------------------------------------------------
namespace py = pybind11;

//----- Python Interface -------------------------------------------------------
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void initMeQuadBlossom(py::module &m) {
  py::class_<xms::MeQuadBlossom, boost::shared_ptr<xms::MeQuadBlossom>> 
    quad_blossom(m, "MeQuadBlossom");

  // ---------------------------------------------------------------------------
  // function: New quad_blossom initializer
  // ---------------------------------------------------------------------------
  const char* init_doc = R"pydoc(
    Make a new quad_blossom instance

    Args:
      ugrid (xms::XmUGrid) The UGrid mesh of triangles to convert to quads.
  )pydoc";
  quad_blossom.def(py::init<>([](boost::shared_ptr<xms::XmUGrid> ugrid)
  {
    return xms::MeQuadBlossom::New(ugrid);
  }));
  // ---------------------------------------------------------------------------
  // function: pre_make_quads
  // ---------------------------------------------------------------------------
  const char* pre_make_quads_doc = R"pydoc(
    Compute edge information and returns the number of boundary edges.
  )pydoc";
  quad_blossom.def("pre_make_quads", [](xms::MeQuadBlossom &self) -> int {
    return self.PreMakeQuads();
  },pre_make_quads_doc);
  // ---------------------------------------------------------------------------
  // function: make_quads
  // ---------------------------------------------------------------------------
  const char* make_quads_doc = R"pydoc(
    Turn faces from triangles into quads by using MeWeightMatcher to identify edges to
    drop and then by calling EliminateEdges to remove them.
    On return, GetPoints() and GetFaces() returns the new points and faces.

    Args:
      split_boundary_points (bool): If true, split boundary points as needed to create quads
        from "pseudo" edges between unmatched boundary triangles separated by at least one
        other triangle.
      use_angle (float): If true use angle calculation to compute the cost for a pair of
        adjacent triangles; otherwise use a diagonal distance calculation.

    Returns:
      (xms::XmUGrid) The new ugrid with quads replacing triangles and any new points.
  )pydoc";
  quad_blossom.def("make_quads", &xms::MeQuadBlossom::MakeQuads,
    make_quads_doc,py::arg("split_boundary_points"),py::arg("use_angle"));
  // ---------------------------------------------------------------------------
  // function: estimate_run_time_in_minutes
  // ---------------------------------------------------------------------------
  const char* estimated_run_time_in_minutes_doc = R"pydoc(
      Get the estimated time to run the Quad Blossom algorithm in minutes.

    Args:
      num_points (int): The number of points in the mesh.

    Returns:
      (double): The estimated minutes to generate the quad mesh.

  )pydoc";
  quad_blossom.def_static("estimated_run_time_in_minutes", &xms::MeQuadBlossom::EstimatedRunTimeInMinutes,
                          estimated_run_time_in_minutes_doc,py::arg("num_points"));
  // ---------------------------------------------------------------------------
  // function: split_to_quads
  // ---------------------------------------------------------------------------
  const char* split_to_quads_doc = R"pydoc(
      Split every cell on the input ugrid at its centroid and edge midpoints to 
      create a new ugrid with only quad cells.

    Args:
      ugrid (XmUGrid): The ugrid to be split into quads.

    Returns:
      (XmUGrid): The new ugrid that contains half size quad cells only.

  )pydoc";
  quad_blossom.def_static("split_to_quads", &xms::MeQuadBlossom::SplitToQuads,
                          split_to_quads_doc, py::arg("ugrid"));
}