"""Test MeMultiPolyMesher_py.cpp."""
import unittest
import numpy as np
from xmscore_py.misc import Observer
from xmsmesh_py.meshing import MeMultiPolyMesher
from xmsmesh_py.meshing import MeMultiPolyMesherIo
from xmsmesh_py.meshing import MePolyInput
from xmsmesh_py.meshing import MeRefinePoint
from xmsinterp_py.interpolate import InterpLinear
from xmsinterp_py.interpolate import InterpIdw


class TestMeMultiPolyMesher(unittest.TestCase):
    """Test MeMultiPolyMesher functions."""

    def setUp(self):
        pass

    def test_creating_MeMultiPolyMesher(self):
        poly_mesher = MeMultiPolyMesher()
        self.assertIsInstance(poly_mesher, MeMultiPolyMesher)

    def test_set_observer(self):
        mesher = MeMultiPolyMesher()
        obs = Observer()
        obs.pytest_id = "ObserverTestMeMultiPolyMesher"
        mesher.set_observer(obs)
        # TODO we don't have a way to test this yet.

    def test_mesh_it(self):
        pass

    def test_check_for_intersections_1(self):
        io = MeMultiPolyMesherIo()
        io.check_topology = True

        poly_input = MePolyInput()
        poly_input.outside_poly = ((0, 0, 0), (100, 0, 0), (100, 10, 0), (0, -10, 0))
        io.poly_inputs = (poly_input,)
        mesher = MeMultiPolyMesher()

        expected = \
            "---Error: Input polygon segments intersect. The segment defined by points 0 and 1 of outer " \
            "polygon 0 intersects with the segment defined by points 2 and 3 of outer polygon 0.\n" \
            "\n\n"

        (success, errors) = mesher.mesh_it(io)
        self.assertEqual(False, success)
        self.assertEqual(expected, errors)

    def test_check_for_intersections_2(self):
        io = MeMultiPolyMesherIo()
        io.check_topology = True

        poly_input = MePolyInput()
        poly_input.outside_poly = ((0, 0, 0), (100, 0, 0), (100, 100, 0), (0, 100, 0))
        poly_input.inside_polys = (((10, 50, 0), (90, 50, 0), (90, 90, 0), (10, 10, 0)),)
        io.poly_inputs = (poly_input,)

        mesher = MeMultiPolyMesher()

        expected = \
            "---Error: Input polygon segments intersect. The segment defined by points 0 and 1 of inner " \
            "polygon 0 of outer polygon 0 intersects with the segment defined by points 2 and 3 of inner " \
            "polygon 0 of outer polygon 0.\n" \
            "\n\n"

        (success, errors) = mesher.mesh_it(io)
        self.assertEqual(False, success)
        self.assertEqual(expected, errors)

    def test_check_for_intersections_3(self):
        io = MeMultiPolyMesherIo()
        io.check_topology = True

        poly_inputs = MePolyInput()
        poly_inputs.outside_poly = ((0, 0, 0), (100, 0, 0), (100, 100, 0), (0, 100, 0))
        poly_inputs.inside_polys = (((90, 10, 0), (110, 10, 0), (110, 20, 0), (90, 20, 0)),)
        io.poly_inputs = (poly_inputs,)

        expected = "---Error: Input polygon segments intersect. The segment defined by points 1 and 2 of outer " \
                   "polygon 0 intersects with the segment defined by points 0 and 1 of inner polygon 0 of outer " \
                   "polygon 0.\n" \
                   "Error: Input polygon segments intersect. The segment defined by points 1 and 2 of outer " \
                   "polygon 0 intersects with the segment defined by points 2 and 3 of inner polygon 0 of outer " \
                   "polygon 0.\n" \
                   "\n\n"

        mesher = MeMultiPolyMesher()
        (success, errors) = mesher.mesh_it(io)
        self.assertEqual(False, success)
        self.assertEqual(expected, errors)

    def test_check_for_intersections_4(self):
        io = MeMultiPolyMesherIo()
        io.check_topology = True

        poly_input1 = MePolyInput()
        poly_input1.outside_poly = ((0, 0, 0), (100, 0, 0), (100, 100, 0), (0, 100, 0))
        poly_input2 = MePolyInput()
        poly_input2.outside_poly = ((10, 10, 0), (110, 10, 0), (110, 110, 0), (10, 110, 0))
        io.poly_inputs = (poly_input1, poly_input2)

        mesher = MeMultiPolyMesher()
        expected = \
            "---Error: Input polygon segments intersect. The segment defined by points 1 and 2 of outer " \
            "polygon 0 intersects with the segment defined by points 0 and 1 of outer polygon 1.\n" \
            "Error: Input polygon segments intersect. The segment defined by points 2 and 3 of outer " \
            "polygon 0 intersects with the segment defined by points 3 and 0 of outer polygon 1.\n" \
            "\n\n"

        (success, errors) = mesher.mesh_it(io)
        self.assertEqual(False, success)
        self.assertEqual(expected, errors)


    def test_check_for_intersections_5(self):
        io = MeMultiPolyMesherIo()
        io.check_topology = True

        poly_input = MePolyInput()
        poly_input.outside_poly = ((0, 0, 0), (100, 0, 0), (100, 100, 0), (0, 100, 0))
        poly_input.inside_polys = (((10, 10, 0), (60, 10, 0), (60, 60, 0), (10, 60, 0)),
                             ((40, 40, 0), (90, 40, 0), (90, 90, 0), (40, 90, 0)))
        io.poly_inputs = (poly_input,)

        mesher = MeMultiPolyMesher()

        expected = \
            "---Error: Input polygon segments intersect. The segment defined by points 1 and 2 of inner " \
            "polygon 0 of outer polygon 0 intersects with the segment defined by points 0 and 1 of inner " \
            "polygon 1 of outer polygon 0.\n" \
            "Error: Input polygon segments intersect. The segment defined by points 2 and 3 of inner " \
            "polygon 0 of outer polygon 0 intersects with the segment defined by points 3 and 0 of inner " \
            "polygon 1 of outer polygon 0.\n" \
            "\n\n"

        (success, errors) = mesher.mesh_it(io)
        self.assertEqual(False, success)
        self.assertEqual(expected, errors)

