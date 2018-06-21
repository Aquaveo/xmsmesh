"""Test MePolyMesher_py.cpp."""
import unittest
import numpy as np
from xmscore_py.misc import Observer
from xmsmesh_py.meshing import MePolyMesher
from xmsmesh_py.meshing import MeMultiPolyMesherIo
from xmsmesh_py.meshing import MePolyInput
from xmsmesh_py.meshing import MeRefinePoint
from xmsinterp_py.interpolate import InterpLinear
from xmsinterp_py.interpolate import InterpIdw


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

    def test_set_observer(self):
        mesher = MePolyMesher()
        obs = Observer()
        mesher.set_observer(obs)
        # TODO: Expand testing here
