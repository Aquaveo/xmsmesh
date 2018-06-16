"""Test InterpLinear_py.cpp."""
import unittest
import numpy as np
from xmsmesh_py.meshing import MeshUtils
from xmsinterp_py.triangulate import TrTin, TrTriangulatorPoints


class TestMeshUtils(unittest.TestCase):
    """Test MeshUtils functions."""

    def setUp(self):
        pass

    def test_size_function_from_depth(self):
        depths = (0, 5, 10, 20, 25, 5, 0)
        min_elem = 2
        max_elem = 102
        sizes = MeshUtils.size_function_from_depth(depths, min_elem, max_elem)
        base_elem_sizes = (2, 22, 42, 82, 102, 22, 2)
        self.assertTupleEqual(base_elem_sizes, sizes)

    def test_size_function_from_depth_numpy(self):
        depths = np.array([0, 5, 10, 20, 25, 5, 0])
        min_elem = 2
        max_elem = 102
        sizes = MeshUtils.size_function_from_depth(depths, min_elem, max_elem)
        base_elem_sizes = np.array([2, 22, 42, 82, 102, 22, 2])
        np.testing.assert_array_equal(base_elem_sizes, sizes)

    def test_smooth_size_func_01(self):
        pts = ((0, 0, 0), (10, 0, 0), (20, 0, 0), (30, 0, 0), (0, 10, 0), (10, 10, 0),
               (20, 10, 0), (30, 10, 0), (0, 20, 0), (10, 20, 0), (20, 20, 0), (30, 20, 0))

        sizes = [100 for _ in range(0, 12)]
        sizes[4] = 1

        tris = ()
        adj_tris = ()

        pts, tris, adj_tris = TrTriangulatorPoints.triangulate(pts, tris, adj_tris)
        tin = TrTin()
        tin.set_geometry(pts, tris, adj_tris)

        size_ratio = 0.5
        min_size = 1.0
        anchor_type = 0
        pt_flags = ()
        smooth_sizes = MeshUtils.smooth_size_function(tin, sizes, size_ratio, min_size,
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

        pts, tris, adj_tris = TrTriangulatorPoints.triangulate(pts, tris, adj_tris)
        tin = TrTin()
        tin.set_geometry(pts, tris, adj_tris)

        size_ratio = 0.5
        min_size = 1.0
        anchor_type = 1
        pt_flags = ()
        smooth_sizes = MeshUtils.smooth_size_function(tin, sizes, size_ratio, min_size,
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

        pts, tris, adj_tris = TrTriangulatorPoints.triangulate(pts, tris, adj_tris)
        tin = TrTin()
        tin.set_geometry(pts, tris, adj_tris)

        min_size = 0.5
        anchor_type = 0
        pt_flags = ()
        smooth_sizes = MeshUtils.smooth_elev_by_slope(tin, sizes, min_size,
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

        pts, tris, adj_tris = TrTriangulatorPoints.triangulate(pts, tris, adj_tris)
        tin = TrTin()
        tin.set_geometry(pts, tris, adj_tris)

        min_size = 0.5
        anchor_type = 1
        pt_flags = ()
        smooth_sizes = MeshUtils.smooth_elev_by_slope(tin, sizes, min_size,
                                                      anchor_type, pt_flags)
        base = (95.00, 92.92, 87.92, 82.92, 100.0, 95.00,
                90.00, 85.00, 95.00, 90.00, 85.00, 80.00)
        np.testing.assert_almost_equal(base, smooth_sizes, 2)
