"""Test MeMultiPolyMesherIo_py.cpp."""
import unittest
import numpy as np
from xmsmesh_py.meshing import MeMultiPolyMesherIo
from xmsmesh_py.meshing import MePolyInput
from xmsmesh_py.meshing import MeRefinePoint
from xmsinterp_py.interpolate import InterpBase
from xmsinterp_py.interpolate import InterpLinear
from xmsinterp_py.interpolate import InterpIdw


class TestMeMultiPolyMesherIo(unittest.TestCase):
    """Test MeMultiPolyMesherIo functions."""

    def setUp(self):
        pass

    def assertInsidePolysEqual(self, base, out):
        self.assertEqual(len(base), len(out), "Base InsidePolys and Out InsidePolys lengths do not match.")
        for i in range(len(base)):
            np.testing.assert_array_equal(base[i], out[i],
                                          "inside_poly[{}] is not the same, -- {} != {}".format(i, base[i], out[i]))

    @staticmethod
    def assertArraysEqual(base, out):
        np.testing.assert_array_equal(np.array(base), out)

    def assertTupleStringsEqual(self, base, out):
        self.assertEqual(len(base), len(out), "Base PolyInputs and Out PolyInputs lengths do not match.")
        for i in range(0, len(base)):
            self.assertEqual(str(base[i]), str(out[i]))


    def test_creating_MeMultiPolyMesherIo(self):
        io = MeMultiPolyMesherIo()
        self.assertIsInstance(io, MeMultiPolyMesherIo)
        self.assertEqual(False, io.check_topology)
        self.assertEqual(True, io.return_cell_polygons)
        self.assertEqual(0, len(io.points))
        self.assertEqual(0, len(io.cells))
        self.assertEqual(0, len(io.cell_polygons))
        self.assertEqual(0, len(io.poly_inputs))
        self.assertEqual(0, len(io.refine_points))

    def test_properties_MeMultiPolyMesherIO(self):
        io = MeMultiPolyMesherIo()
        self.assertIsInstance(io, MeMultiPolyMesherIo)

        io.check_topology = True
        self.assertEqual(True, io.check_topology)

        io.return_cell_polygons = False
        self.assertEqual(False, io.return_cell_polygons)

        points = ((1, 1, 2), (1, 2, 3), (2, 3, 4), (3, 4, 5))
        io.points = points
        self.assertArraysEqual(points, io.points)

        cells = (1, 5, 9, 13)
        io.cells = cells
        self.assertArraysEqual(cells, io.cells)

        cell_polygons = (2, 4)
        io.cell_polygons = cell_polygons
        self.assertArraysEqual(cell_polygons, io.cell_polygons)

        pi1 = MePolyInput()
        pi1.bias = 2.718
        pi2 = MePolyInput()
        pi2.bias = 0.618
        io.poly_inputs = (pi1, pi2)
        self.assertTupleStringsEqual((pi1, pi2), io.poly_inputs)

        rp1 = MeRefinePoint((5, 0, -3), 3.1, False)
        rp2 = MeRefinePoint((-2, -2, 1), -0.4, True)
        io.refine_points = (rp1, rp2)
        self.assertTupleStringsEqual((rp1, rp2), io.refine_points)

class TestMePolyInput(unittest.TestCase):
    """Test MePolyInput functions."""

    def setUp(self):
        pass

    def assertInsidePolysEqual(self, base, out):
        self.assertEqual(len(base), len(out), "Base InsidePolys and Out InsidePolys lengths do not match.")
        for i in range(len(base)):
            np.testing.assert_array_equal(base[i], out[i],
                                          "inside_poly[{}] is not the same, -- {} != {}".format(i, base[i], out[i]))

    @staticmethod
    def assertArraysEqual(base, out):
        np.testing.assert_array_equal(np.array(base), out)

    def test_creating_default_MePolyInput(self):
        pi = MePolyInput()
        self.assertIsInstance(pi, MePolyInput)
        self.assertEqual(0, len(pi.outside_poly))
        self.assertEqual(0, len(pi.inside_polys))
        self.assertEqual(1.0, pi.bias)
        self.assertEqual(None, pi.size_function)
        self.assertEqual(None, pi.elev_function)
        self.assertEqual(-1, pi.const_size_bias)
        self.assertEqual(-1, pi.const_size_function)
        self.assertEqual(False, pi.remove_internal_four_triangle_pts)

    def test_creating_MePolyInput(self):
        outside_poly = ((1, 2, 0), (5, 2, 0), (5, 9, 0), (1, 9, 0))
        inside_polys = (((3, 3, 0), (2.5, 4, 0), (2, 3, 0)),
                        ((4, 8, 0), (3, 7, 0), (2, 8, 0))
                        )
        poly_corners = (0, 1, 2, 3)
        bias = 3.14159
        size_func = InterpLinear()
        elev_func = InterpIdw()

        pi = MePolyInput(outside_poly, inside_polys, bias, size_func, poly_corners, elev_func)
        self.assertIsInstance(pi, MePolyInput)
        self.assertArraysEqual(outside_poly, pi.outside_poly)
        self.assertInsidePolysEqual(inside_polys, pi.inside_polys)
        self.assertEqual(bias, pi.bias)
        self.assertEqual(size_func.to_string(), pi.size_function.to_string())
        self.assertEqual(elev_func.to_string(), pi.elev_function.to_string())
        self.assertEqual(-1, pi.const_size_bias)
        self.assertEqual(-1, pi.const_size_function)
        self.assertEqual(False, pi.remove_internal_four_triangle_pts)

    def test_properties_MePolyInput(self):
        pi = MePolyInput()
        outside_poly = ((1, 2, 0), (5, 2, 0), (5, 9, 0), (1, 9, 0))
        inside_polys = (((3, 3, 0), (2.5, 4, 0), (2, 3, 0)),
                        ((4, 8, 0), (3, 7, 0), (2, 8, 0))
                        )
        poly_corners = (0, 1, 2, 3)
        size_func = InterpLinear()
        elev_func = InterpIdw()

        self.assertEqual(0, len(pi.outside_poly))
        pi.outside_poly = outside_poly
        self.assertArraysEqual(outside_poly, pi.outside_poly)

        self.assertEqual(0, len(pi.inside_polys))
        pi.inside_polys = inside_polys
        self.assertInsidePolysEqual(inside_polys, pi.inside_polys)

        self.assertEqual(0, len(pi.poly_corners))
        pi.poly_corners = poly_corners
        self.assertArraysEqual(poly_corners, pi.poly_corners)

        self.assertEqual(1.0, pi.bias)
        pi.bias = 0.3
        self.assertEqual(0.3, pi.bias)

        self.assertEqual(None, pi.size_function)
        pi.size_function = size_func
        self.assertEqual(size_func.to_string(), pi.size_function.to_string())

        self.assertEqual(None, pi.elev_function)
        pi.elev_function = elev_func
        self.assertEqual(elev_func.to_string(), pi.elev_function.to_string())

        self.assertEqual(-1, pi.const_size_bias)
        pi.const_size_bias = 4.0
        self.assertEqual(4.0, pi.const_size_bias)

        self.assertEqual(-1, pi.const_size_function)
        pi.const_size_function = 1.2
        self.assertEqual(1.2, pi.const_size_function)

        self.assertEqual(False, pi.remove_internal_four_triangle_pts)
        pi.remove_internal_four_triangle_pts = True
        self.assertEqual(True, pi.remove_internal_four_triangle_pts)

class TestMeRefinePoint(unittest.TestCase):
    """Test MeRefinePoint functions."""

    def setUp(self):
        pass

    def test_creating_MeRefinePoint(self):
        rp = MeRefinePoint((1, 2, 3), -2.0, True)
        self.assertIsInstance(rp, MeRefinePoint)
        self.assertEqual((1, 2, 3), rp.point)
        self.assertEqual(-2.0, rp.size)
        self.assertEqual(True, rp.create_mesh_point)

    def test_properties_MeRefinePoint(self):
        rp = MeRefinePoint((4, 5, 3), 2.0, False)

        self.assertEqual(2.0, rp.size)
        rp.size = 4.5
        self.assertEqual(4.5, rp.size)

        self.assertEqual(False, rp.create_mesh_point)
        rp.create_mesh_point = True
        self.assertEqual(True, rp.create_mesh_point)

        self.assertEqual((4, 5, 3), rp.point)
        rp.point = (-2.0, 4.15, -900.001)
        self.assertEqual((-2.0, 4.15, -900.001), rp.point)

