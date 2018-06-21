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

    def test_interp_edge_lengths(self):
        out_poly = ((0, 0, 0), (0, 10, 0), (10, 10, 0), (10, 0, 0))
        in_polys = ()
        r = MePolyRedistributePts()
        size_bias = 1.0
        r.set_size_func_from_poly(out_poly, in_polys, size_bias)

        pts = ((1, 1, 0), (1, 9, 0), (9, 9, 0), (9, 1, 0))
        # interp_edge_lengths is not exposed, only on the C++ impl
        # lengths = (,)
        # r.interp_edge_lengths(pts, lengths)

        # base_lengths = (10.0, 10.0, 10.0)
        # self.assertArraysEqual(base_lengths, lengths)
