"""Test InterpLinear_py.cpp."""
import numpy as np
import os
import unittest

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

        sizes = [100 for _ in range(0, 12)]
        sizes[4] = 1

        tris = ()
        adj_tris = ()

        tin = Tin(pts, tris)
        tin.set_triangles_adjacent_to_points(adj_tris)
        tin.triangulate()

        min_size = 0.5
        anchor_type = 0
        pt_flags = ()
        smooth_sizes = mesh_utils.smooth_elev_by_slope(tin, sizes, min_size,
                                                      anchor_type, pt_flags)
        base = (6.00,  8.07,  13.07, 18.07, 1.0,   6.00,
                11.00, 16.00, 6.00,  11.00, 16.00, 21.00)
        np.testing.assert_almost_equal(base, smooth_sizes, 2)

    def test_smooth_elev_by_slope_02(self):
        pts = ((0, 0, 0), (10, 0, 0), (20, 0, 0), (30, 0, 0), (0, 10, 0), (10, 10, 0),
               (20, 10, 0), (30, 10, 0), (0, 20, 0), (10, 20, 0), (20, 20, 0), (30, 20, 0))

        sizes = [1 for _ in range(0, 12)]
        sizes[4] = 100

        tris = ()
        adj_tris = ()

        tin = Tin(pts, tris)
        tin.set_triangles_adjacent_to_points(adj_tris)
        tin.triangulate()

        min_size = 0.5
        anchor_type = 1
        pt_flags = ()
        smooth_sizes = mesh_utils.smooth_elev_by_slope(tin, sizes, min_size,
                                                      anchor_type, pt_flags)
        base = (95.00, 92.92, 87.92, 82.92, 100.0, 95.00,
                90.00, 85.00, 95.00, 90.00, 85.00, 80.00)
        np.testing.assert_almost_equal(base, smooth_sizes, 2)

    def test_check_for_intersections_1(self):
        io = MultiPolyMesherIo(())
        io.check_topology = True

        poly_input = PolyInput(out_poly=((0, 0, 0), (100, 0, 0), (100, 10, 0), (0, -10, 0)))
        io.poly_inputs = (poly_input,)

        expected = \
            "---Error: Input polygon segments intersect. The segment defined by points 0 and 1 of outer " \
            "polygon 0 intersects with the segment defined by points 2 and 3 of outer polygon 0.\n" \
            "\n\n"

        (success, errors) = mesh_utils.generate_triangle_mesh(io)
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

        (success, errors) = mesh_utils.generate_triangle_mesh(io)
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

        (success, errors) = mesh_utils.generate_triangle_mesh(io)
        self.assertEqual(False, success)
        self.assertEqual(expected, errors)

    def test_check_for_intersections_4(self):
        io = MultiPolyMesherIo(())
        io.check_topology = True

        outside_poly1 = ((0, 0, 0), (100, 0, 0), (100, 100, 0), (0, 100, 0))
        outside_poly2 = ((10, 10, 0), (110, 10, 0), (110, 110, 0), (10, 110, 0))
        poly_input1 = PolyInput(outside_poly1)
        poly_input2 = PolyInput(outside_poly2)
        io.poly_inputs = (poly_input1, poly_input2)

        expected = \
            "---Error: Input polygon segments intersect. The segment defined by points 1 and 2 of outer " \
            "polygon 0 intersects with the segment defined by points 0 and 1 of outer polygon 1.\n" \
            "Error: Input polygon segments intersect. The segment defined by points 2 and 3 of outer " \
            "polygon 0 intersects with the segment defined by points 3 and 0 of outer polygon 1.\n" \
            "\n\n"

        (success, errors) = mesh_utils.generate_triangle_mesh(io)
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

        (success, errors) = mesh_utils.generate_triangle_mesh(io)
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
        status, error = mesh_utils.generate_triangle_mesh(input)
        self.assertTrue(status)
        self.assertEqual(error, '')

    def test_generate_2dm(self):
        io = MultiPolyMesherIo(())
        _ = mesh_utils.generate_2dm(io, "fname.2dm")
        self.assertTrue(os.path.isfile("fname.2dm"))

    def test_case_4(self):
        # build test case 4 polys
        out_a = (0, 60, 10, 60, 20, 60, 30, 60, 30, 50, 30, 40, 30, 30, 30, 20, 30, 10,
                 30, 0, 20, 0, 10, 0, 0, 0, 0, 10, 0, 20, 0, 30, 0, 40, 0, 50)
        in_a1 = (10, 50, 10, 40, 20, 40, 20, 50)
        in_a2 = (10, 20, 10, 10, 20, 10, 20, 20)
        outside_poly = self.array_to_vec_pt3d(out_a)
        inside_polys = (self.array_to_vec_pt3d(in_a1), self.array_to_vec_pt3d(in_a2))
        bias = 1.0
        poly_input_a = PolyInput(out_poly=outside_poly, inside_polys=inside_polys, bias=bias)

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
        (success, result) = mesh_utils.generate_2dm(io, "", 8)

        expected = "MESH2D\n" \
                   "E3T     1     1    18     2     1\n" \
                   "E3T     2     2    18    24     1\n" \
                   "E3T     3     2    24     3     1\n" \
                   "E3T     4     3    23     4     1\n" \
                   "E3T     5     3    24    23     1\n" \
                   "E3T     6     4    23    27     1\n" \
                   "E3T     7     4    27     5     1\n" \
                   "E3T     8     5    20     6     1\n" \
                   "E3T     9     5    27    20     1\n" \
                   "E3T    10     6    19     7     1\n" \
                   "E3T    11     6    20    19     1\n" \
                   "E3T    12     7    19     8     1\n" \
                   "E3T    13     8    19    22     1\n" \
                   "E3T    14     8    22     9     1\n" \
                   "E3T    15     9    11    10     1\n" \
                   "E3T    16     9    22    11     1\n" \
                   "E3T    17    10    11    33     1\n" \
                   "E3T    18    10    33    34     1\n" \
                   "E3T    19    11    22    12     1\n" \
                   "E3T    20    12    13    31     1\n" \
                   "E3T    21    12    21    28     1\n" \
                   "E3T    22    12    22    21     1\n" \
                   "E3T    23    12    28    13     1\n" \
                   "E3T    24    12    31    32     1\n" \
                   "E3T    25    13    26    14     1\n" \
                   "E3T    26    13    28    26     1\n" \
                   "E3T    27    14    15    47     1\n" \
                   "E3T    28    14    26    15     1\n" \
                   "E3T    29    14    47    30     1\n" \
                   "E3T    30    15    25    16     1\n" \
                   "E3T    31    15    26    25     1\n" \
                   "E3T    32    16    25    17     1\n" \
                   "E3T    33    17    24    18     1\n" \
                   "E3T    34    17    25    24     1\n" \
                   "E3T    35    20    27    28     1\n" \
                   "E3T    36    20    28    21     1\n" \
                   "E3T    37    23    26    29     1\n" \
                   "E3T    38    23    29    27     1\n" \
                   "E3T    39    26    28    29     1\n" \
                   "E3T    40    27    29    28     1\n" \
                   "E3T    41    30    47    53     1\n" \
                   "E3T    42    30    52    31     1\n" \
                   "E3T    43    30    53    52     1\n" \
                   "E3T    44    31    52    56     1\n" \
                   "E3T    45    31    56    32     1\n" \
                   "E3T    46    32    49    33     1\n" \
                   "E3T    47    32    56    49     1\n" \
                   "E3T    48    33    48    34     1\n" \
                   "E3T    49    33    49    48     1\n" \
                   "E3T    50    34    48    35     1\n" \
                   "E3T    51    35    48    51     1\n" \
                   "E3T    52    35    51    36     1\n" \
                   "E3T    53    36    38    37     1\n" \
                   "E3T    54    36    51    38     1\n" \
                   "E3T    55    38    51    39     1\n" \
                   "E3T    56    39    50    57     1\n" \
                   "E3T    57    39    51    50     1\n" \
                   "E3T    58    39    57    40     1\n" \
                   "E3T    59    40    55    41     1\n" \
                   "E3T    60    40    57    55     1\n" \
                   "E3T    61    41    55    42     1\n" \
                   "E3T    62    42    54    43     1\n" \
                   "E3T    63    42    55    54     1\n" \
                   "E3T    64    43    54    44     1\n" \
                   "E3T    65    44    53    45     1\n" \
                   "E3T    66    44    54    53     1\n" \
                   "E3T    67    45    47    46     1\n" \
                   "E3T    68    45    53    47     1\n" \
                   "E3T    69    49    56    57     1\n" \
                   "E3T    70    49    57    50     1\n" \
                   "E3T    71    52    55    58     1\n" \
                   "E3T    72    52    58    56     1\n" \
                   "E3T    73    55    57    58     1\n" \
                   "E3T    74    56    58    57     1\n" \
                   "ND     1      0.0      0.0      0.0\n" \
                   "ND     2      0.0     10.0      0.0\n" \
                   "ND     3      0.0     20.0      0.0\n" \
                   "ND     4      0.0     30.0      0.0\n" \
                   "ND     5      0.0     40.0      0.0\n" \
                   "ND     6      0.0     50.0      0.0\n" \
                   "ND     7      0.0     60.0      0.0\n" \
                   "ND     8     10.0     60.0      0.0\n" \
                   "ND     9     20.0     60.0      0.0\n" \
                   "ND    10     30.0     60.0      0.0\n" \
                   "ND    11     30.0     50.0      0.0\n" \
                   "ND    12     30.0     40.0      0.0\n" \
                   "ND    13     30.0     30.0      0.0\n" \
                   "ND    14     30.0     20.0      0.0\n" \
                   "ND    15     30.0     10.0      0.0\n" \
                   "ND    16     30.0      0.0      0.0\n" \
                   "ND    17     20.0      0.0      0.0\n" \
                   "ND    18     10.0      0.0      0.0\n" \
                   "ND    19     10.0     50.0      0.0\n" \
                   "ND    20     10.0     40.0      0.0\n" \
                   "ND    21     20.0     40.0      0.0\n" \
                   "ND    22     20.0     50.0      0.0\n" \
                   "ND    23     10.0     20.0      0.0\n" \
                   "ND    24     10.0     10.0      0.0\n" \
                   "ND    25     20.0     10.0      0.0\n" \
                   "ND    26     20.0     20.0      0.0\n" \
                   "ND    27  8.51082 31.75796      0.0\n" \
                   "ND    28  19.9189 32.42692      0.0\n" \
                   "ND    29 14.90632 26.65924      0.0\n" \
                   "ND    30     40.0     20.0      0.0\n" \
                   "ND    31     40.0     30.0      0.0\n" \
                   "ND    32     40.0     40.0      0.0\n" \
                   "ND    33     40.0     50.0      0.0\n" \
                   "ND    34     40.0     60.0      0.0\n" \
                   "ND    35     50.0     60.0      0.0\n" \
                   "ND    36     60.0     60.0      0.0\n" \
                   "ND    37     70.0     60.0      0.0\n" \
                   "ND    38     70.0     50.0      0.0\n" \
                   "ND    39     70.0     40.0      0.0\n" \
                   "ND    40     70.0     30.0      0.0\n" \
                   "ND    41     70.0     20.0      0.0\n" \
                   "ND    42     70.0     10.0      0.0\n" \
                   "ND    43     70.0      0.0      0.0\n" \
                   "ND    44     60.0      0.0      0.0\n" \
                   "ND    45     50.0      0.0      0.0\n" \
                   "ND    46     40.0      0.0      0.0\n" \
                   "ND    47     40.0     10.0      0.0\n" \
                   "ND    48     50.0     50.0      0.0\n" \
                   "ND    49     50.0     40.0      0.0\n" \
                   "ND    50     60.0     40.0      0.0\n" \
                   "ND    51     60.0     50.0      0.0\n" \
                   "ND    52     50.0     20.0      0.0\n" \
                   "ND    53     50.0     10.0      0.0\n" \
                   "ND    54     60.0     10.0      0.0\n" \
                   "ND    55     60.0     20.0      0.0\n" \
                   "ND    56 48.37114 31.77082      0.0\n" \
                   "ND    57 59.82848 32.41625      0.0\n" \
                   "ND    58 54.80369 26.66427      0.0\n"
        self.assertTrue(success)
        self.assertEqual(expected, result)

