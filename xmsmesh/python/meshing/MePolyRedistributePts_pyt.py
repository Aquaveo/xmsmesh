"""Test PolyRedistributePts_py.cpp."""
import unittest
import numpy as np
from xmscore.misc import Observer
from xmsmesh.meshing import PolyRedistributePts
from xmsmesh.meshing import MultiPolyMesherIo
from xmsmesh.meshing import PolyInput
from xmsmesh.meshing import RefinePoint
from xmsinterp_py.interpolate import InterpLinear
from xmsinterp_py.interpolate import InterpIdw


class TestPolyRedistributePts(unittest.TestCase):
    """Test PolyRedistributePts functions."""

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

    def test_creating_PolyRedistributePts(self):
        mesher = PolyRedistributePts()
        self.assertIsInstance(mesher, PolyRedistributePts)

    def test_set_size_func_01(self):
        r = PolyRedistributePts()
        sf = InterpLinear()
        r.set_size_func(sf)
        # TODO: No way to test if there size function was set correctly

    def test_set_size_func_02(self):
        r = PolyRedistributePts()
        sf = InterpIdw()
        r.set_size_func(sf)
        # TODO: No way to test if there size function was set correctly

    def test_set_size_fun_from_poly(self):
        out_poly = ((0, 0, 0), (0, 10, 0), (10, 10, 0), (10, 0, 0))
        in_polys = ()
        r = PolyRedistributePts()
        size_bias = 1.0
        r.set_size_func_from_poly(out_poly, in_polys, size_bias)
        # TODO: No way to test if there size function was set correctly

    def test_constant_size_func(self):
        r = PolyRedistributePts()
        r.set_constant_size_func(0.75)
        poly_line = ((0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0))
        new_poly_line = r.redistribute(poly_line)
        base_poly_line = [
            [0., 0., 0.], [0.75, 0., 0.], [1.5, 0., 0.],
            [2.25, 0., 0.], [3., 0., 0.]
        ]
        np.testing.assert_array_equal(base_poly_line, new_poly_line)

    # TODO: This Crashes when using size bias
    # def test_constant_size_bias(self):
    #     r = PolyRedistributePts()
    #     r.set_constant_size_bias(1.5)
    #     poly_line = ((0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0))
    #     new_poly_line = r.redistrubte(poly_line)
    #     base_poly_line = []  # Haven't gotten a base yet
    #     np.testing.assert_array_equal(base_poly_line, new_poly_line)
   
    def test_set_use_curvature_redistribution(self):
        r = PolyRedistributePts()
        r.set_use_curvature_redistribution(3.0, 4.0, 0.001, False)
        poly_line = ((0, 0, 0), (5, 5, 0), (10, 10, 0), (15, 5, 0),
                     (20, 10, 0), (21, 10, 0), (25, 0, 0))
        new_poly_line = r.redistribute(poly_line)
        base_poly_line = ((0, 0, 0), (9.961, 9.961, 0), (10.457, 9.543, 0),
                          (14.892, 5.108, 0), (15.388, 5.388, 0),  (19.685, 9.685, 0),
                          (19.987, 9.987, 0), (20.906, 10.000, 0), (21.122, 9.695, 0),
                          (25.000, 0.000, 0))
        np.testing.assert_array_almost_equal(base_poly_line, new_poly_line, decimal=2)

    def test_redistribute(self):
        r = PolyRedistributePts()
        r.set_constant_size_func(7)
        # r.set_constant_size_bias(3.5)
        poly_line = [(x, 0, 0) for x in range(0, 100, 2)]
        new_poly_line = r.redistribute(poly_line)
        base_poly_line = np.array([
            [0., 0., 0.], [7., 0., 0.], [14., 0., 0.], [21., 0., 0.], [28., 0., 0.],
            [35., 0., 0.], [42., 0., 0.], [49., 0., 0.], [56., 0., 0.], [63., 0., 0.],
            [70., 0., 0.], [77., 0., 0.], [84., 0., 0.], [91., 0., 0.], [98., 0., 0.]
        ])
        np.testing.assert_array_almost_equal(np.array(base_poly_line), new_poly_line, decimal=7)

