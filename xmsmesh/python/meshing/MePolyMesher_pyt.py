"""Test MePolyMesher_py.cpp."""
import unittest
import numpy as np
from xmscore_py.misc import Observer
from xmsmesh_py.meshing import MePolyMesher
from xmsmesh_py.meshing import MeMultiPolyMesherIo
from xmsmesh_py.meshing import MePolyInput


class TestMePolyMesher(unittest.TestCase):
    """Test MePolyMesher functions."""

    def setUp(self):
        pass

    @staticmethod
    def array_to_vec_pt3d(a_array):
        return [(a_array[i], a_array[i+1], 0) for i in range(0, len(a_array), 2)]

    @staticmethod
    def merge(points, indices):
        return [(points[i]) for i in indices]

    @staticmethod
    def assertArraysEqual(base, out):
        np.testing.assert_array_equal(np.array(base), out)

    def test_creating_MePolyMesher(self):
        mesher = MePolyMesher()
        self.assertIsInstance(mesher, MePolyMesher)

    def test_mesh_it(self):
        outside_poly = ((5, 0, 0), (10, 5, 0), (10, 15, 0),
                        (5, 20, 0), (0, 15, 0), (0, 5, 0))
        inside_polys = (((10, 5, 0), (10, 15, 0), (0, 15, 0), (0, 5, 0)),)
        poly_input = MePolyInput()
        poly_input.outside_poly = outside_poly
        poly_input.inside_polys = inside_polys
        mesher = MePolyMesher()
        mesher_io = MeMultiPolyMesherIo()
        mesher_io.poly_inputs = (poly_input,)
        status, points, tris, cells = mesher.mesh_it(mesher_io, 0)
        self.assertEqual(True, status)
        np.testing.assert_array_equal([[0.,  5.,  0.], [5.,  0.,  0.], [10.,  5.,  0.],
                                       [10., 15.,  0.], [5., 20.,  0.], [0., 15.,  0.]], points)
        np.testing.assert_array_equal((2, 0, 1, 5, 3, 4), tris)
        np.testing.assert_array_equal((), cells)

    def test_get_processed_refine_pts(self):
        # TODO: How do we test this?
        pass
