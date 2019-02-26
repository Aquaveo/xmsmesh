"""Test InterpLinear_py.cpp."""
import numpy as np
import os
import unittest
import filecmp

from xmsinterp.triangulate import Tin

from xmsmesh.meshing import mesh_utils
from xmsmesh.meshing import MultiPolyMesherIo
from xmsmesh.meshing import PolyInput


class TestMeshUtils(unittest.TestCase):
    """Test MeshUtils functions."""
    @staticmethod
    def array_to_vec_pt3d(a_array):
        return [(a_array[i], a_array[i+1], 0) for i in range(0, len(a_array), 2)]

    def setUp(self):
        pass

    def test_size_function_from_depth(self):
        depths = (0, 5, 10, 20, 25, 5, 0)
        min_elem = 2
        max_elem = 102
        sizes = mesh_utils.size_function_from_depth(depths, min_elem, max_elem)
        base_elem_sizes = (2, 22, 42, 82, 102, 22, 2)
        self.assertTupleEqual(base_elem_sizes, sizes)

    def test_size_function_from_depth_numpy(self):
        depths = np.array([0, 5, 10, 20, 25, 5, 0])
        min_elem = 2
        max_elem = 102
        sizes = mesh_utils.size_function_from_depth(depths, min_elem, max_elem)
        base_elem_sizes = np.array([2, 22, 42, 82, 102, 22, 2])
        np.testing.assert_array_equal(base_elem_sizes, sizes)

    def test_smooth_size_func_01(self):
        pts = ((0, 0, 0), (10, 0, 0), (20, 0, 0), (30, 0, 0), (0, 10, 0), (10, 10, 0),
               (20, 10, 0), (30, 10, 0), (0, 20, 0), (10, 20, 0), (20, 20, 0), (30, 20, 0))

        sizes = [100 for _ in range(0, 12)]
        sizes[4] = 1

        tris = ()
        adj_tris = ()

        tin = Tin(pts, tris)
        tin.set_triangles_adjacent_to_points(adj_tris)
        tin.triangulate()

        size_ratio = 0.5
        min_size = 1.0
        anchor_type = 0
        pt_flags = ()
        smooth_sizes = mesh_utils.smooth_size_function(tin, sizes, size_ratio, min_size,
                                                      anchor_type, pt_flags)
        base = (4.46, 5.90,  9.36, 12.83, 1.0,   4.46,
                7.93, 11.39, 4.46, 7.93,  11.39, 14.86)
        np.testing.assert_almost_equal(base, smooth_sizes, 2)

    def test_smooth_size_func_02(self):
        pts = ((0, 0, 0), (10, 0, 0), (20, 0, 0), (30, 0, 0), (0, 10, 0), (10, 10, 0),
               (20, 10, 0), (30, 10, 0), (0, 20, 0), (10, 20, 0), (20, 20, 0), (30, 20, 0))

        sizes = [1 for _ in range(0, 12)]
        sizes[4] = 100

        tris = ()
        adj_tris = ()

        tin = Tin(pts, tris)
        tin.set_triangles_adjacent_to_points(adj_tris)
        tin.triangulate()

        size_ratio = 0.5
        min_size = 1.0
        anchor_type = 1
        pt_flags = ()
        smooth_sizes = mesh_utils.smooth_size_function(tin, sizes, size_ratio, min_size,
                                                      anchor_type, pt_flags)
        base = (96.53, 95.10, 91.63, 88.17, 100.0, 96.53,
                93.07, 89.60, 96.53, 93.07, 89.60, 86.14)
        np.testing.assert_almost_equal(base, smooth_sizes, 2)

    def test_smooth_elev_by_slope_01(self):
        pts = ((0, 0, 0), (10, 0, 0), (20, 0, 0), (30, 0, 0), (0, 10, 0), (10, 10, 0),
               (20, 10, 0), (30, 10, 0), (0, 20, 0), (10, 20, 0), (20, 20, 0), (30, 20, 0))

        elevations = [100 for _ in range(0, 12)]
        elevations[4] = 1

        tris = ()
        adj_tris = ()

        tin = Tin(pts, tris)
        tin.set_triangles_adjacent_to_points(adj_tris)
        tin.triangulate()

        min_size = 0.5
        anchor_type = 0
        pt_flags = ()
        smooth_sizes = mesh_utils.smooth_elev_by_slope(tin, elevations, min_size,
                                                      anchor_type, pt_flags)
        base = (6.00,  8.07,  13.07, 18.07, 1.0,   6.00,
                11.00, 16.00, 6.00,  11.00, 16.00, 21.00)
        np.testing.assert_almost_equal(base, smooth_sizes, 2)

    def test_smooth_elev_by_slope_02(self):
        pts = ((0, 0, 0), (10, 0, 0), (20, 0, 0), (30, 0, 0), (0, 10, 0), (10, 10, 0),
               (20, 10, 0), (30, 10, 0), (0, 20, 0), (10, 20, 0), (20, 20, 0), (30, 20, 0))

        elevations = [1 for _ in range(0, 12)]
        elevations[4] = 100

        tris = ()
        adj_tris = ()

        tin = Tin(pts, tris)
        tin.set_triangles_adjacent_to_points(adj_tris)
        tin.triangulate()

        min_size = 0.5
        anchor_type = 1
        pt_flags = ()
        smooth_sizes = mesh_utils.smooth_elev_by_slope(tin, elevations, min_size,
                                                      anchor_type, pt_flags)
        base = (95.00, 92.92, 87.92, 82.92, 100.0, 95.00,
                90.00, 85.00, 95.00, 90.00, 85.00, 80.00)
        np.testing.assert_almost_equal(base, smooth_sizes, 2)

    def test_check_for_intersections(self):
        io = MultiPolyMesherIo(())
        poly_input = PolyInput(outside_polygon=((0, 0, 0), (100, 0, 0), (100, 10, 0), (0, -10, 0)))
        io.poly_inputs = (poly_input,)

        expected = \
            "Error: Input polygon segments intersect. The segment defined by points 0 and 1 of outer " \
            "polygon 0 intersects with the segment defined by points 2 and 3 of outer polygon 0.\n"

        (success, errors) = mesh_utils.check_mesh_input_topology(io)
        self.assertEqual(False, success)
        self.assertEqual(expected, errors)

    def test_check_for_intersections_1(self):
        io = MultiPolyMesherIo(())
        io.check_topology = True

        poly_input = PolyInput(outside_polygon=((0, 0, 0), (100, 0, 0), (100, 10, 0), (0, -10, 0)))
        io.poly_inputs = (poly_input,)

        expected = \
            "---Error: Input polygon segments intersect. The segment defined by points 0 and 1 of outer " \
            "polygon 0 intersects with the segment defined by points 2 and 3 of outer polygon 0.\n" \
            "\n\n"

        (success, errors) = mesh_utils.generate_mesh(io)
        self.assertEqual(False, success)
        self.assertEqual(expected, errors)

    def test_check_for_intersections_2(self):
        io = MultiPolyMesherIo(())
        io.check_topology = True

        outside_poly = ((0, 0, 0), (100, 0, 0), (100, 100, 0), (0, 100, 0))
        inside_polys = (((10, 50, 0), (90, 50, 0), (90, 90, 0), (10, 10, 0)),)
        poly_input = PolyInput(outside_poly, inside_polys)
        io.poly_inputs = (poly_input,)

        expected = \
            "---Error: Input polygon segments intersect. The segment defined by points 0 and 1 of inner " \
            "polygon 0 of outer polygon 0 intersects with the segment defined by points 2 and 3 of inner " \
            "polygon 0 of outer polygon 0.\n" \
            "\n\n"

        (success, errors) = mesh_utils.generate_mesh(io)
        self.assertEqual(False, success)
        self.assertEqual(expected, errors)

    def test_check_for_intersections_3(self):
        io = MultiPolyMesherIo(())
        io.check_topology = True

        outside_poly = ((0, 0, 0), (100, 0, 0), (100, 100, 0), (0, 100, 0))
        inside_polys = (((90, 10, 0), (110, 10, 0), (110, 20, 0), (90, 20, 0)),)
        poly_inputs = PolyInput(outside_poly, inside_polys)
        io.poly_inputs = (poly_inputs,)

        expected = "---Error: Input polygon segments intersect. The segment defined by points 1 and 2 of outer " \
                   "polygon 0 intersects with the segment defined by points 0 and 1 of inner polygon 0 of outer " \
                   "polygon 0.\n" \
                   "Error: Input polygon segments intersect. The segment defined by points 1 and 2 of outer " \
                   "polygon 0 intersects with the segment defined by points 2 and 3 of inner polygon 0 of outer " \
                   "polygon 0.\n" \
                   "\n\n"

        (success, errors) = mesh_utils.generate_mesh(io)
        self.assertEqual(False, success)
        self.assertEqual(expected, errors)

    def test_check_for_intersections_4(self):
        io = MultiPolyMesherIo(())
        io.check_topology = True

        outside_poly1 = ((0, 0, 0), (0, 100, 0), (100, 100, 0), (100, 0, 0))
        outside_poly2 = ((10, 10, 0), (10, 110, 0), (110, 110, 0), (110, 10, 0))
        poly_input1 = PolyInput(outside_poly1)
        poly_input2 = PolyInput(outside_poly2)
        io.poly_inputs = (poly_input1, poly_input2)

        expected = \
            "---Error: Input polygon segments intersect. The segment defined by points 1 and 2 of outer " \
            "polygon 0 intersects with the segment defined by points 0 and 1 of outer polygon 1.\n" \
            "Error: Input polygon segments intersect. The segment defined by points 2 and 3 of outer " \
            "polygon 0 intersects with the segment defined by points 3 and 0 of outer polygon 1.\n" \
            "\n\n"

        (success, errors) = mesh_utils.generate_mesh(io)
        self.assertEqual(False, success)
        self.assertEqual(expected, errors)

    def test_check_for_intersections_5(self):
        io = MultiPolyMesherIo(())
        io.check_topology = True

        outside_poly = ((0, 0, 0), (100, 0, 0), (100, 100, 0), (0, 100, 0))
        inside_polys = (((10, 10, 0), (60, 10, 0), (60, 60, 0), (10, 60, 0)),
                             ((40, 40, 0), (90, 40, 0), (90, 90, 0), (40, 90, 0)))
        poly_input = PolyInput(outside_poly, inside_polys)
        io.poly_inputs = (poly_input,)

        expected = \
            "---Error: Input polygon segments intersect. The segment defined by points 1 and 2 of inner " \
            "polygon 0 of outer polygon 0 intersects with the segment defined by points 0 and 1 of inner " \
            "polygon 1 of outer polygon 0.\n" \
            "Error: Input polygon segments intersect. The segment defined by points 2 and 3 of inner " \
            "polygon 0 of outer polygon 0 intersects with the segment defined by points 3 and 0 of inner " \
            "polygon 1 of outer polygon 0.\n" \
            "\n\n"

        (success, errors) = mesh_utils.generate_mesh(io)
        self.assertEqual(False, success)
        self.assertEqual(expected, errors)

    def test_simple_polygon(self):
        outside_poly = [
            (0, 10, 0), (0, 20, 0), (0, 30, 0), (0, 40, 0), (0, 50, 0), (0, 60, 0), (0, 70, 0), (0, 80, 0),
            (0, 90, 0), (0, 100, 0), (10, 100, 0), (20, 100, 0), (30, 100, 0), (40, 100, 0), (50, 100, 0), (60, 100, 0),
            (70, 100, 0), (80, 100, 0), (90, 100, 0), (100, 100, 0), (100, 90, 0), (100, 80, 0), (100, 70, 0),
            (100, 60, 0), (100, 50, 0), (100, 40, 0), (100, 30, 0), (100, 20, 0), (100, 10, 0), (100, 0, 0),
            (90, 0, 0), (80, 0, 0), (70, 0, 0), (60, 0, 0), (50, 0, 0), (40, 0, 0), (30, 0, 0), (20, 0, 0), (10, 0, 0),
            (0, 0, 0)
        ]
        inside_polys = [
            [(40, 40, 0), (50, 40, 0), (60, 40, 0), (60, 50, 0),
             (60, 60, 0), (50, 60, 0), (40, 60, 0), (40, 50, 0)]
        ]
        input_poly = PolyInput(outside_poly, inside_polys)
        input = MultiPolyMesherIo(())
        input.poly_inputs = [input_poly]
        status, error = mesh_utils.generate_mesh(input)
        self.assertEqual(139, len(input.points))
        self.assertEqual(1150, len(input.cells))
        self.assertTrue(status)
        self.assertEqual(error, '')

    def test_simple_polygon_reverse(self):
        outside_poly = [
            (0, 10, 0), (0, 20, 0), (0, 30, 0), (0, 40, 0), (0, 50, 0), (0, 60, 0), (0, 70, 0), (0, 80, 0),
            (0, 90, 0), (0, 100, 0), (10, 100, 0), (20, 100, 0), (30, 100, 0), (40, 100, 0), (50, 100, 0), (60, 100, 0),
            (70, 100, 0), (80, 100, 0), (90, 100, 0), (100, 100, 0), (100, 90, 0), (100, 80, 0), (100, 70, 0),
            (100, 60, 0), (100, 50, 0), (100, 40, 0), (100, 30, 0), (100, 20, 0), (100, 10, 0), (100, 0, 0),
            (90, 0, 0), (80, 0, 0), (70, 0, 0), (60, 0, 0), (50, 0, 0), (40, 0, 0), (30, 0, 0), (20, 0, 0), (10, 0, 0),
            (0, 0, 0)
        ]
        inside_polys = [
            [(40, 40, 0), (50, 40, 0), (60, 40, 0), (60, 50, 0),
             (60, 60, 0), (50, 60, 0), (40, 60, 0), (40, 50, 0)]
        ]
        outside_poly.reverse()
        inside_polys[0].reverse()
        input_poly = PolyInput(outside_poly, inside_polys)
        input = MultiPolyMesherIo(())
        input.poly_inputs = [input_poly]
        status, error = mesh_utils.generate_mesh(input)
        self.assertEqual(139, len(input.points))
        self.assertEqual(1150, len(input.cells))
        self.assertTrue(status)
        self.assertEqual(error, '')

    def test_generate_2dm(self):
        io = MultiPolyMesherIo(())
        _ = mesh_utils.generate_2dm(io, "fname.2dm")
        self.assertTrue(os.path.isfile("fname.2dm"))
        self.assertTrue(filecmp.cmp("../test_files/python/fname.2dm", "fname.2dm"), "Files not equal")

    def test_case_4(self):
        # build test case 4 polys
        out_a = (0, 60, 10, 60, 20, 60, 30, 60, 30, 50, 30, 40, 30, 30, 30, 20, 30, 10,
                 30, 0, 20, 0, 10, 0, 0, 0, 0, 10, 0, 20, 0, 30, 0, 40, 0, 50)
        in_a1 = (10, 50, 10, 40, 20, 40, 20, 50)
        in_a2 = (10, 20, 10, 10, 20, 10, 20, 20)
        outside_poly = self.array_to_vec_pt3d(out_a)
        inside_polys = (self.array_to_vec_pt3d(in_a1), self.array_to_vec_pt3d(in_a2))
        bias = 1.0
        poly_input_a = PolyInput(outside_polygon=outside_poly, inside_polygons=inside_polys, bias=bias)

        out_b = (30, 60, 40, 60, 50, 60, 60, 60, 70, 60, 70, 50, 70, 40, 70, 30,
                 70, 20, 70, 10, 70, 0, 60, 0, 50, 0, 40, 0, 40, 10, 30, 10,
                 30, 20, 40, 20, 40, 30, 30, 30, 30, 40, 40, 40, 40, 50, 30, 50)
        in_b1 = (50, 50, 50, 40, 60, 40, 60, 50)
        in_b2 = (50, 20, 50, 10, 60, 10, 60, 20)

        outside_poly_b = self.array_to_vec_pt3d(out_b)
        inside_polys_b = (self.array_to_vec_pt3d(in_b1), self.array_to_vec_pt3d(in_b2))
        bias_b = 1.0
        poly_input_b = PolyInput(outside_poly_b, inside_polys_b, bias=bias_b)

        io = MultiPolyMesherIo(())
        io.poly_inputs = (poly_input_a, poly_input_b)

        # Value Error when file not specified
        with self.assertRaises(ValueError) as context:
            mesh_utils.generate_2dm(io, '')
        print("'{}'".format(context.exception))
        self.assertTrue('file_name not specifed. Aborting mesh procedure.' in str(context.exception))

        # mesh the polys
        (success, result) = mesh_utils.generate_2dm(io, "out_file.2dm", 3)
        self.assertTrue(success)
        self.assertTrue(os.path.isfile("out_file.2dm"))
        self.assertTrue(filecmp.cmp("../test_files/python/out_file.2dm", "out_file.2dm"), "Files not equal")

    def test_repeated_first_and_last(self):
        # build test case 4 polys
        out_a = (0, 60, 10, 60, 20, 60, 30, 60, 30, 50, 30, 40, 30, 30, 30, 20, 30, 10,
                 30, 0, 20, 0, 10, 0, 0, 0, 0, 10, 0, 20, 0, 30, 0, 40, 0, 50, 0, 60)
        in_a1 = (10, 50, 10, 40, 20, 40, 20, 50, 10, 50)
        in_a2 = (10, 20, 10, 10, 20, 10, 20, 20, 10, 20)
        outside_poly = self.array_to_vec_pt3d(out_a)
        inside_polys = (self.array_to_vec_pt3d(in_a1), self.array_to_vec_pt3d(in_a2))
        bias = 1.0
        poly_input_a = PolyInput(outside_polygon=outside_poly, inside_polygons=inside_polys, bias=bias)

        out_b = (30, 60, 40, 60, 50, 60, 60, 60, 70, 60, 70, 50, 70, 40, 70, 30,
                 70, 20, 70, 10, 70, 0, 60, 0, 50, 0, 40, 0, 40, 10, 30, 10,
                 30, 20, 40, 20, 40, 30, 30, 30, 30, 40, 40, 40, 40, 50, 30, 50)
        in_b1 = (50, 50, 50, 40, 60, 40, 60, 50)
        in_b2 = (50, 20, 50, 10, 60, 10, 60, 20)

        outside_poly_b = self.array_to_vec_pt3d(out_b)
        inside_polys_b = (self.array_to_vec_pt3d(in_b1), self.array_to_vec_pt3d(in_b2))
        bias_b = 1.0
        poly_input_b = PolyInput(outside_poly_b, inside_polys_b, bias=bias_b)

        io = MultiPolyMesherIo(())
        io.poly_inputs = (poly_input_a, poly_input_b)

        # mesh the polys
        (success, result) = mesh_utils.generate_2dm(io, "out_file_02.2dm", 8)
        self.assertTrue(success)
        self.assertTrue(os.path.isfile("out_file_02.2dm"))
        self.assertTrue(filecmp.cmp("../test_files/python/out_file.2dm", "out_file_02.2dm"), "Files not equal")


    def test_redistribute_polyline(self):
        polygon_corners = [(0, 0, 0), (0, 100, 0), (100, 100, 0),
                           (100, 0, 0), (0, 0, 0)]

        polygon_boundary = mesh_utils.redistribute_poly_line(polygon_corners, 25)
        base_poly_boundary = [(0, 0, 0), (0, 25, 0), (0, 50, 0), (0, 75, 0), (0, 100, 0),
                              (25, 100, 0), (50, 100, 0), (75, 100, 0), (100, 100, 0),
                              (100, 75, 0), (100, 50, 0), (100, 25, 0), (100, 0, 0),
                              (75, 0, 0), (50, 0, 0), (25, 0, 0), (0, 0, 0)]
        np.testing.assert_array_equal(base_poly_boundary, polygon_boundary)

