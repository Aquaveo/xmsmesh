//------------------------------------------------------------------------------
/// \file
/// \brief
/// \copyright (C) Copyright Aquaveo 2018. Distributed under FreeBSD License
/// (See accompanying file LICENSE or https://aqaveo.com/bsd/license.txt)
//------------------------------------------------------------------------------

//----- Included files ---------------------------------------------------------
#include <pybind11/pybind11.h>
#include <iostream>
#include <sstream>
#include <pybind11/numpy.h>
#include <boost/shared_ptr.hpp>

#include <xmscore/misc/boost_defines.h>
#include <xmscore/python/misc/PyUtils.h>

#include <xmsinterp/interpolate/InterpLinear.h>

#include <xmsmesh/python/meshing/meshing_py.h>
#include <xmsmesh/meshing/MeMultiPolyMesherIo.h>

//----- Namespace declaration --------------------------------------------------
namespace py = pybind11;

//----- Python Interface -------------------------------------------------------
PYBIND11_DECLARE_HOLDER_TYPE(T, BSHP<T>);

void initMeMultiPolyMesherIo(py::module &m) {

    const char* multi_poly_mesher_doc_init = R"pydoc(
        Creates a mesh from one or more PolyInputs, and other settings that are defined
        in the PolyMesherIo that is passed into this class.

        Args:
            poly_inputs (collections.Iterable[xmsmesh.meshing.PolyInput]): A list of PolyInputs.
            refine_points (collections.Iterable[xmsmesh.meshing.RefinePoint]): A list of RefinePoints.
            check_topology (bool): Check polygon input topology for errors.
            return_cell_polygons (bool): Return the polygon index of each cell.
    )pydoc";

    py::class_<xms::MeMultiPolyMesherIo, BSHP<xms::MeMultiPolyMesherIo>> polyMesherIo(m, "MultiPolyMesherIo");

    polyMesherIo.def(py::init<>([](py::iterable poly_inputs, py::iterable refine_points,
            bool check_topology, bool return_cell_polygons) {
      BSHP<xms::MeMultiPolyMesherIo> rval(new xms::MeMultiPolyMesherIo());
      if (!poly_inputs.is_none())
      {
        std::vector<xms::MePolyInput> &vecPolys = rval->m_polys;
        vecPolys.clear();
        vecPolys.reserve(py::len(poly_inputs));
        for (auto item : poly_inputs) {
          vecPolys.push_back(item.cast<xms::MePolyInput>());
        }
      }
      if (!refine_points.is_none())
      {
        std::vector<xms::MeRefinePoint> &vecRefinePoints = rval->m_refPts;
        vecRefinePoints.clear();
        vecRefinePoints.reserve(py::len(refine_points));
        for (auto item : refine_points) {
          vecRefinePoints.push_back(item.cast<xms::MeRefinePoint>());
        }
      }
      rval->m_checkTopology = check_topology;
      rval->m_returnCellPolygons = return_cell_polygons;
      return rval;
    }), multi_poly_mesher_doc_init, py::arg("poly_inputs"), py::arg("refine_points") = py::make_tuple(),
        py::arg("check_topology") = false, py::arg("return_cell_polygons") = true);
    // ---------------------------------------------------------------------------
    // function: check_topology
    // ---------------------------------------------------------------------------
    const char* check_topology_doc = R"pydoc(
        If True, checks PolyInput topology for errors.
    )pydoc";
    polyMesherIo.def_readwrite("check_topology", &xms::MeMultiPolyMesherIo::m_checkTopology,check_topology_doc);
    // ---------------------------------------------------------------------------
    // function: return_cell_polygons
    // ---------------------------------------------------------------------------
    const char* return_cell_polygons_doc = R"pydoc(
        If True, the cell_polygons list will be be filled when meshing occurs.
    )pydoc";
    polyMesherIo.def_readwrite("return_cell_polygons", 
        &xms::MeMultiPolyMesherIo::m_returnCellPolygons, 
        return_cell_polygons_doc);
    // ---------------------------------------------------------------------------
    // function: points
    // ---------------------------------------------------------------------------
    const char* points_doc = R"pydoc(
        A list of (x, y, z) coordinates of the resulting mesh. (Populated by meshing functions)
    )pydoc";
    polyMesherIo.def_property("points",
        [](xms::MeMultiPolyMesherIo &self) -> py::iterable {
            return xms::PyIterFromVecPt3d(self.m_points);
        },
        [](xms::MeMultiPolyMesherIo &self, py::iterable outside_polygon) {
            self.m_points = *xms::VecPt3dFromPyIter(outside_polygon);
        },points_doc);
    // ---------------------------------------------------------------------------
    // function: cells
    // ---------------------------------------------------------------------------
    const char* cells_doc = R"pydoc(
        A cell stream representing the mesh. (Populated by meshing functions)
    )pydoc";
    polyMesherIo.def_property("cells",
            [](xms::MeMultiPolyMesherIo &self) -> py::iterable {
                return xms::PyIterFromVecInt(self.m_cells);
            },
            [](xms::MeMultiPolyMesherIo &self, py::iterable cells) {
                 self.m_cells = *xms::VecIntFromPyIter(cells);
            },
            cells_doc);
    // ---------------------------------------------------------------------------
    // function: cell_polygons
    // ---------------------------------------------------------------------------
    const char* cell_polygons_doc = R"pydoc(
        The index of the PolyInput in cell_polygons that each cell was generated from.
    )pydoc";
    polyMesherIo.def_property("cell_polygons",
            [](xms::MeMultiPolyMesherIo &self) -> py::iterable {
                return xms::PyIterFromVecInt(self.m_cellPolygons);
            },
            [](xms::MeMultiPolyMesherIo &self, py::iterable cell_polygons) {
                self.m_cellPolygons = *xms::VecIntFromPyIter(cell_polygons);
            },cell_polygons_doc);
    // ---------------------------------------------------------------------------
    // function: poly_inputs
    // ---------------------------------------------------------------------------
    const char* poly_inputs_doc = R"pydoc(
        A list of PolyInput objects to be used by the meshing functions.
    )pydoc";
    polyMesherIo.def_property("poly_inputs",
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
            },poly_inputs_doc);
    // ---------------------------------------------------------------------------
    // function: refine_points
    // ---------------------------------------------------------------------------
    const char* refine_points_doc = R"pydoc(
        A list of RefinePoint objects used for mesh refinement.
    )pydoc";
    polyMesherIo.def_property("refine_points",
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
            },refine_points_doc);
    // -------------------------------------------------------------------------
    // function: __repr__
    // -------------------------------------------------------------------------
    polyMesherIo.def("__repr__", [](xms::MeMultiPolyMesherIo &self) {
        std::stringstream ss;
        ss << "PolyMesherIo(s):\n";
        for (int i = 0; i < self.m_polys.size(); ++i)
        {
          ss << "PolyInput #" << i + 1 << ":\n";
          ss << PyReprStringFromMePolyInput(self.m_polys[i]) << "\n";
        }
        ss << "RefinePoint(s):\n";
        for (int i = 0; i < self.m_polys.size(); ++i)
        {
          ss << "RefinePoint #" << i + 1 << ":\n";
          ss << PyReprStringFromMeRefinePoint(self.m_refPts[i]) << "\n";
        }
        std::string offOn[2] = {"False", "True"};
        ss << "Check Topology: " << offOn[(int)self.m_checkTopology] << "\n";
        ss << "Return Cell Polygons: " << offOn[(int)self.m_returnCellPolygons] << "\n";
        return ss.str();
    });
}
    // -------------------------------------------------------------------------
    // function: __init__
    // -------------------------------------------------------------------------
void initMePolyInput(py::module &m) {

    const char* PolyInput_init_doc = R"pydoc(
        PolyInput initializer

        Args:
            outside_polygon (iterable): Point locations of outer polygon. Clockwise
            inside_polygons (iterable optional): Point locations of inner polygons (holes). Counter clockwise. 1st pt != last. Defaults to empty tuple.
            bias (Float optional): Factor for transitioning between areas of high refinement to less refinement. Defaults to 0.3.
            const_size_bias (Float optional): Transition factor for constant size function.
            const_size_function (Float optional): Constat to be used for size function.
            bound_pts_to_remove (iterable optional): Outer boundary locations to remove after the paving process.
            relaxation_method (str optional): The relaxation method to be used.
            remove_internal_four_triangle_pts (bool optional): Remove internal points that are only connected to 4 cells.
            seed_points (iterable optional): A list of seed points.
            size_function (:class:`InterpBase <xmsinterp.interpolate.InterpBase>` optional): Size function for scalar paving. Default to None
            patch_polygon_corners (iterable optional): Corner nodes for creating meshes using the patch algorithm. 3 per outer poly (not 4 - outer poly index point [0] is assumed to be a corner). Defaults to empty tuple.
            elev_function (:class:`InterpBase <xmsinterp.interpolate.InterpBase>` optional): Elevation function for interpolating z coordinate of mesh points. Deafults to None.

    )pydoc";

    py::class_<xms::MePolyInput, BSHP<xms::MePolyInput>> polyInput(m, "PolyInput");

    polyInput.def(py::init<>([](py::iterable outside_polygon, py::iterable inside_polygons, double bias,
                            py::object const_size_bias, py::object const_size_function, py::object bound_pts_to_remove,
                            py::object relaxation_method, py::object remove_internal_four_triangles_pts, py::object seed_points,
                            py::object size_function, py::iterable patch_polygon_corners, py::object elev_function) {
            xms::VecPt3d vec_outside_polygon = *xms::VecPt3dFromPyIter(outside_polygon);
            xms::VecPt3d2d vec_inside_polygons = *xms::VecPt3d2dFromPyIter(inside_polygons);
            xms::VecInt vec_poly_corners = *xms::VecIntFromPyIter(patch_polygon_corners);
            BSHP<xms::InterpBase> c_size_function;
            BSHP<xms::InterpBase> c_elev_function;
            if (!size_function.is_none())
            {
              c_size_function = size_function.cast<BSHP<xms::InterpBase>>();
            }
            if (!elev_function.is_none())
            {
              c_elev_function = elev_function.cast<BSHP<xms::InterpBase>>();
            }
            BSHP<xms::MePolyInput> rval(new xms::MePolyInput(vec_outside_polygon, vec_inside_polygons, bias, c_size_function,
                                                         vec_poly_corners, c_elev_function));
            if(!const_size_bias.is_none())
            {
                rval->m_constSizeBias = const_size_bias.cast<float>();
            }
            if(!const_size_function.is_none())
            {
                rval->m_constSizeFunction = const_size_function.cast<float>();
            }
            if(!bound_pts_to_remove.is_none())
            {
                // TODO: This might need to be list rather than none as default.
                rval->m_boundPtsToRemove = *xms::VecPt3dFromPyIter(bound_pts_to_remove);
            }
            if(!relaxation_method.is_none())
            {
                rval->m_relaxationMethod = relaxation_method.cast<std::string>();
            }
            if(!remove_internal_four_triangles_pts.is_none())
            {
                rval->m_removeInternalFourTrianglePts = remove_internal_four_triangles_pts.cast<bool>();
            }
            if(!seed_points.is_none())
            {
                // TODO: This might need to be list rather than none as default.
                rval->m_seedPoints = *xms::VecPt3dFromPyIter(seed_points);
            }
            return rval;
        }), PolyInput_init_doc, py::arg("outside_polygon"), py::arg("inside_polygons") = py::make_tuple(), py::arg("bias") = 1.0,
            py::arg("const_size_bias") = py::none(), py::arg("const_size_fuction") = py::none(),
            py::arg("bound_pts_to_remove") = py::none(), py::arg("relaxation_method") = py::none(),
            py::arg("remove_interal_four_triangle_pts") = py::none(), py::arg("seed_points") = py::none(),
            py::arg("size_function") = py::none(), py::arg("patch_polygon_corners") = py::make_tuple(),
            py::arg("elev_function") = py::none());
    // ---------------------------------------------------------------------------
    // function: outside_polygon
    // ---------------------------------------------------------------------------
    const char* outside_polygon_doc = R"pydoc(
        List of points defining the outside polygon.

        Warning:
            These points must be in clockwise order, and the first point must not equal the last point.
    )pydoc";
    polyInput.def_property("outside_polygon",
            [](xms::MePolyInput &self) -> py::iterable {
                return xms::PyIterFromVecPt3d(self.m_outPoly);
            },
            [](xms::MePolyInput &self, py::iterable outside_polygon) {
                 self.m_outPoly = * xms::VecPt3dFromPyIter(outside_polygon);
            },outside_polygon_doc);
    // ---------------------------------------------------------------------------
    // function: inside_polygons
    // ---------------------------------------------------------------------------
    const char* inside_polygons_doc = R"pydoc(
        A list of polygons representing holes in the PolyInput

        The polygons should be in clockwise order and the first point must not equal the last point.
    )pydoc";
    polyInput.def_property("inside_polygons",
            [](xms::MePolyInput &self) -> py::iterable {
                return xms::PyIterFromVecPt3d2d(self.m_insidePolys);
            },
            [](xms::MePolyInput &self, py::iterable inside_polygons) {
                self.m_insidePolys = *xms::VecPt3d2dFromPyIter(inside_polygons);
            },inside_polygons_doc);
    // -------------------------------------------------------------------------
    // function: patch_polygon_corners
    // -------------------------------------------------------------------------
    const char* patch_polygon_corners_doc = R"pydoc(
        Corner nodes for creating meshes using the patch algorithm.

        There can be 3 patch_polygon_corners per outer_poly not 4. The outer_poly point at index 0 is assumed to be a corner
    )pydoc";
    polyInput.def_property("patch_polygon_corners",
            [](xms::MePolyInput &self) -> py::iterable {
                return xms::PyIterFromVecInt(self.m_polyCorners);
            },
            [](xms::MePolyInput &self, py::iterable patch_polygon_corners) {
                 self.m_polyCorners = *xms::VecIntFromPyIter(patch_polygon_corners);
            },patch_polygon_corners_doc);
    // -------------------------------------------------------------------------
    // function: bound_pts_to_remove
    // -------------------------------------------------------------------------
    const char* bound_pts_to_remove_doc = R"pydoc(
        Outer boundary locations to remove after the paving process.
    )pydoc";
    polyInput.def_property("bound_pts_to_remove",
            [](xms::MePolyInput &self) -> py::iterable {
                return xms::PyIterFromVecPt3d(self.m_boundPtsToRemove);
            },
            [](xms::MePolyInput &self, py::iterable outside_polygon) {
                 self.m_boundPtsToRemove = *xms::VecPt3dFromPyIter(outside_polygon);
            },bound_pts_to_remove_doc);
    // -------------------------------------------------------------------------
    // function: bias
    // -------------------------------------------------------------------------
    const char* bias_doc = R"pydoc(
        Factor for transitioning between areas of high refinement to less refinement.
    )pydoc";
    polyInput.def_readwrite("bias", &xms::MePolyInput::m_bias, bias_doc);
    // -------------------------------------------------------------------------
    // function: size_function
    // -------------------------------------------------------------------------
    const char* size_function_doc = R"pydoc(
        Size function for scalar paving.
    )pydoc";
    polyInput.def_readwrite("size_function", &xms::MePolyInput::m_sizeFunction,
            size_function_doc);
    // -------------------------------------------------------------------------
    // function: elev_function
    // -------------------------------------------------------------------------
    const char* elev_function_doc = R"pydoc(
        Elevation function for interpolating z coordinate of mesh points.
    )pydoc";
    polyInput.def_readwrite("elev_function", &xms::MePolyInput::m_elevFunction,
        elev_function_doc);
    // -------------------------------------------------------------------------
    // function: const_size_function
    // -------------------------------------------------------------------------
    const char* const_size_function_doc = R"pydoc(
        Constant to be used for size function.
    )pydoc";
    polyInput.def_readwrite("const_size_function", &xms::MePolyInput::m_constSizeFunction,
        const_size_function_doc);
    // -------------------------------------------------------------------------
    // function: const_size_bias
    // -------------------------------------------------------------------------
    const char* const_size_bias_doc = R"pydoc(
        Transition factor for constant size function.
    )pydoc";
    polyInput.def_readwrite("const_size_bias", &xms::MePolyInput::m_constSizeBias,
        const_size_bias_doc);
    // -------------------------------------------------------------------------
    // function: remove_internal_four_triangle_pts
    // -------------------------------------------------------------------------
    const char* remove_internal_four_triangle_pts_doc = R"pydoc(
        Remove internal points that are only connected to 4 cells.
    )pydoc";
    polyInput.def_readwrite("remove_internal_four_triangle_pts", &xms::MePolyInput::m_removeInternalFourTrianglePts,
            remove_internal_four_triangle_pts_doc);
    // -------------------------------------------------------------------------
    // function: poly_id
    // -------------------------------------------------------------------------
//    const char* poly_id_doc = R"pydoc(
//        Optional. Set when needed. Can be useful for classes who need an ID.
//    )pydoc";
    polyInput.def_readwrite("poly_id", &xms::MePolyInput::m_polyId);
    // -------------------------------------------------------------------------
    // function: seed_points
    // -------------------------------------------------------------------------
    const char* seed_points_doc = R"pydoc(
        A list of seed points. If the user has some methodology for
        creating points inside the polygon then those points can be specified 
        here. If these points are specified then the paving is not performed.
        These points will not be used if the meshing option is patch.
    )pydoc";
    polyInput.def_property("seed_points",
            [](xms::MePolyInput &self) -> py::iterable {
                return xms::PyIterFromVecPt3d(self.m_seedPoints);
            },
            [](xms::MePolyInput &self, py::iterable seed_points) {
                 self.m_seedPoints = *xms::VecPt3dFromPyIter(seed_points);
            }, seed_points_doc);
    // ---------------------------------------------------------------------------
    // property: relaxation_method
    // ---------------------------------------------------------------------------
    const char* relaxation_method_doc = R"pydoc(
        Relaxation method. The default relaxation method is an area
        relax. Set the value to "spring_relaxation". See MeRelaxer.cpp for
        details on spring relaxation.
    )pydoc";
    polyInput.def_readwrite("relaxation_method",
      &xms::MePolyInput::m_relaxationMethod,relaxation_method_doc);
    // -------------------------------------------------------------------------
    // function: __str__
    // -------------------------------------------------------------------------
    const char* __str__doc = R"pydoc(
        outputs contents as string

        Returns:
            str: Contents as a string
    )pydoc";
    polyInput.def("__str__", [](xms::MePolyInput &self) {
             std::string szf = self.m_sizeFunction == nullptr ? "none" : self.m_sizeFunction->ToString();
             std::string elevf = self.m_elevFunction == nullptr ? "none" : self.m_elevFunction->ToString();
             std::stringstream ss;
             ss <<  "bias: " << self.m_bias << std::endl <<
                    "size_func: " << szf << std::endl <<
                    "elev_func: " << elevf << std::endl <<
                    "const_size_function: " << self.m_constSizeFunction << std::endl <<
                    "const_size_bias: " << self.m_constSizeBias << std::endl <<
                    "outside_polygon size: " << self.m_outPoly.size() << std::endl <<
                    "inside_polygons size: " << self.m_insidePolys.size() << std::endl;
             return ss.str();
        },__str__doc);
    // -------------------------------------------------------------------------
    // function: __repr__
    // -------------------------------------------------------------------------
    polyInput.def("__repr__", [](xms::MePolyInput &self) {
      return PyReprStringFromMePolyInput(self);
    });
}

void initMeRefinePoint(py::module &m) {

    const char* refine_point_doc = R"pydoc(
        A location and options used for refining a mesh.

        Args:
            point (iterable): An (x, y, z) location for the point.
            size (float): A size for the refinement. A negative value represents a hard point.
            create_mesh_point (bool optional): Force a mesh point at the refine point location.
    )pydoc";

    py::class_<xms::MeRefinePoint, BSHP<xms::MeRefinePoint>> refinePoint(m, "RefinePoint");

    refinePoint.def(py::init<>([](py::tuple point, double size, bool create_mesh_point) {
            xms::Pt3d _point = xms::Pt3dFromPyIter(point);
            return new xms::MeRefinePoint(_point, size, create_mesh_point);
        }),refine_point_doc, py::arg("point"), py::arg("size"), py::arg("create_mesh_point") = true);
    // -------------------------------------------------------------------------
    // function: point
    // -------------------------------------------------------------------------
    const char* point_doc = R"pydoc(
        Location of refine refine point

    )pydoc";
    refinePoint.def_property("point",
          [](xms::MeRefinePoint &self) -> py::tuple {
            return py::make_tuple(self.m_pt.x, self.m_pt.y, self.m_pt.z);
          },
          [](xms::MeRefinePoint &self, py::tuple pt) {
            xms::Pt3d point = xms::Pt3dFromPyIter(pt);
            self.m_pt = point;
          },point_doc);
    // -------------------------------------------------------------------------
    // function: size
    // -------------------------------------------------------------------------
    const char* size_doc = R"pydoc(
        Element size at the refine point.
        point.
    )pydoc";
    refinePoint.def_readwrite("size", &xms::MeRefinePoint::m_size,size_doc);
    // -------------------------------------------------------------------------
    // function: create_mesh_point
    // -------------------------------------------------------------------------
    const char* create_mesh_point_doc = R"pydoc(
        Create a mesh point at this location with attached cells at the specified size.
    )pydoc";
    refinePoint.def_readwrite("create_mesh_point", &xms::MeRefinePoint::m_createMeshPoint,
      create_mesh_point_doc);
    // -------------------------------------------------------------------------
    // function: __str__
    // -------------------------------------------------------------------------
    const char* __str___doc = R"pydoc(
        Returns contents as string

        Returns:
            str: contents as a string
    )pydoc";
    refinePoint.def("__str__", [](xms::MeRefinePoint &self) {
            std::stringstream ss;
            ss << "point: (" << self.m_pt << ")" <<  std::endl <<
                   "size: " << self.m_size <<  std::endl <<
                   "create_mesh_point: " << self.m_createMeshPoint << std::endl;
            return ss.str();
        }, __str___doc);
    // -------------------------------------------------------------------------
    // function: __repr__
    // -------------------------------------------------------------------------
    refinePoint.def("__repr__", [](xms::MeRefinePoint &self) {
      return PyReprStringFromMeRefinePoint(self);
    });
}
