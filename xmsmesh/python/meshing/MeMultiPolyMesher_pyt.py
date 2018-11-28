"""Test MultiPolyMesher_py.cpp."""
import unittest
import numpy as np
from xmscore_py.misc import Observer
from xmsmesh.meshing import MultiPolyMesher
from xmsmesh.meshing import MultiPolyMesherIo
from xmsmesh.meshing import PolyInput
from xmsmesh.meshing import RefinePoint
from xmsinterp_py.interpolate import InterpLinear
from xmsinterp_py.interpolate import InterpIdw


class TestMultiPolyMesher(unittest.TestCase):
    """Test MultiPolyMesher functions."""

    def setUp(self):
        pass

    def test_creating_MultiPolyMesher(self):
        poly_mesher = MultiPolyMesher()
        self.assertIsInstance(poly_mesher, MultiPolyMesher)

    def test_mesh_it(self):
        pass

    def test_check_for_intersections_1(self):
        io = MultiPolyMesherIo()
        io.check_topology = True

        poly_input = PolyInput()
        poly_input.outside_poly = ((0, 0, 0), (100, 0, 0), (100, 10, 0), (0, -10, 0))
        io.poly_inputs = (poly_input,)
        mesher = MultiPolyMesher()

        expected = \
            "---Error: Input polygon segments intersect. The segment defined by points 0 and 1 of outer " \
            "polygon 0 intersects with the segment defined by points 2 and 3 of outer polygon 0.\n" \
            "\n\n"

        (success, errors) = mesher.mesh_it(io)
        self.assertEqual(False, success)
        self.assertEqual(expected, errors)

    def test_check_for_intersections_2(self):
        io = MultiPolyMesherIo()
        io.check_topology = True

        poly_input = PolyInput()
        poly_input.outside_poly = ((0, 0, 0), (100, 0, 0), (100, 100, 0), (0, 100, 0))
        poly_input.inside_polys = (((10, 50, 0), (90, 50, 0), (90, 90, 0), (10, 10, 0)),)
        io.poly_inputs = (poly_input,)

        mesher = MultiPolyMesher()

        expected = \
            "---Error: Input polygon segments intersect. The segment defined by points 0 and 1 of inner " \
            "polygon 0 of outer polygon 0 intersects with the segment defined by points 2 and 3 of inner " \
            "polygon 0 of outer polygon 0.\n" \
            "\n\n"

        (success, errors) = mesher.mesh_it(io)
        self.assertEqual(False, success)
        self.assertEqual(expected, errors)

    def test_check_for_intersections_3(self):
        io = MultiPolyMesherIo()
        io.check_topology = True

        poly_inputs = PolyInput()
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

        mesher = MultiPolyMesher()
        (success, errors) = mesher.mesh_it(io)
        self.assertEqual(False, success)
        self.assertEqual(expected, errors)

    def test_check_for_intersections_4(self):
        io = MultiPolyMesherIo()
        io.check_topology = True

        poly_input1 = PolyInput()
        poly_input1.outside_poly = ((0, 0, 0), (100, 0, 0), (100, 100, 0), (0, 100, 0))
        poly_input2 = PolyInput()
        poly_input2.outside_poly = ((10, 10, 0), (110, 10, 0), (110, 110, 0), (10, 110, 0))
        io.poly_inputs = (poly_input1, poly_input2)

        mesher = MultiPolyMesher()
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
        io = MultiPolyMesherIo()
        io.check_topology = True

        poly_input = PolyInput()
        poly_input.outside_poly = ((0, 0, 0), (100, 0, 0), (100, 100, 0), (0, 100, 0))
        poly_input.inside_polys = (((10, 10, 0), (60, 10, 0), (60, 60, 0), (10, 60, 0)),
                             ((40, 40, 0), (90, 40, 0), (90, 90, 0), (40, 90, 0)))
        io.poly_inputs = (poly_input,)

        mesher = MultiPolyMesher()

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

    def test_simple_polygon(self):
        input_poly = PolyInput()
        input_poly.outside_poly = [
            (0, 10, 0), (0, 20, 0), (0, 30, 0), (0, 40, 0), (0, 50, 0), (0, 60, 0), (0, 70, 0), (0, 80, 0),
            (0, 90, 0), (0, 100, 0), (10, 100, 0), (20, 100, 0), (30, 100, 0), (40, 100, 0), (50, 100, 0), (60, 100, 0),
            (70, 100, 0), (80, 100, 0), (90, 100, 0), (100, 100, 0), (100, 90, 0), (100, 80, 0), (100, 70, 0),
            (100, 60, 0), (100, 50, 0), (100, 40, 0), (100, 30, 0), (100, 20, 0), (100, 10, 0), (100, 0, 0),
            (90, 0, 0), (80, 0, 0), (70, 0, 0), (60, 0, 0), (50, 0, 0), (40, 0, 0), (30, 0, 0), (20, 0, 0), (10, 0, 0),
            (0, 0, 0)
        ]
        input_poly.inside_polys = [
            [(40, 40, 0), (50, 40, 0), (60, 40, 0), (60, 50, 0),
             (60, 60, 0), (50, 60, 0), (40, 60, 0), (40, 50, 0)]
        ]
        input = MultiPolyMesherIo()
        input.poly_inputs = [input_poly]
        mesher = MultiPolyMesher()
        status, error = mesher.mesh_it(input)
        self.assertTrue(status)
        self.assertEqual(error, '')
