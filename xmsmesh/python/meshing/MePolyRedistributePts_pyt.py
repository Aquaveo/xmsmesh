"""Test MePolyRedistributePts_py.cpp."""
import unittest
import numpy as np
from xmscore_py.misc import Observer
from xmsmesh_py.meshing import MePolyRedistributePts
from xmsmesh_py.meshing import MeMultiPolyMesherIo
from xmsmesh_py.meshing import MePolyInput
from xmsmesh_py.meshing import MeRefinePoint
from xmsinterp_py.interpolate import InterpLinear
from xmsinterp_py.interpolate import InterpIdw


class TestMePolyRedistributePts(unittest.TestCase):
    """Test MePolyRedistributePts functions."""

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

    def test_creating_MePolyRedistributePts(self):
        mesher = MePolyRedistributePts()
        self.assertIsInstance(mesher, MePolyRedistributePts)

    def test_set_size_func_01(self):
        r = MePolyRedistributePts()
        sf = InterpLinear()
        r.set_size_func(sf)
        # TODO: No way to test if there size function was set correctly

    def test_set_size_func_02(self):
        r = MePolyRedistributePts()
        sf = InterpIdw()
        r.set_size_func(sf)
        # TODO: No way to test if there size function was set correctly

    def test_set_size_fun_from_poly(self):
        out_poly = ((0, 0, 0), (0, 10, 0), (10, 10, 0), (10, 0, 0))
        in_polys = ()
        r = MePolyRedistributePts()
        size_bias = 1.0
        r.set_size_func_from_poly(out_poly, in_polys, size_bias)
        # TODO: No way to test if there size function was set correctly

    def test_constant_size_func(self):
        r = MePolyRedistributePts()
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
    #     r = MePolyRedistributePts()
    #     r.set_constant_size_bias(1.5)
    #     poly_line = ((0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0))
    #     new_poly_line = r.redistrubte(poly_line)
    #     base_poly_line = []  # Haven't gotten a base yet
    #     np.testing.assert_array_equal(base_poly_line, new_poly_line)

    def test_redistribute(self):
        r = MePolyRedistributePts()
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

