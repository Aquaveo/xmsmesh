"""Test MeQuadBlossom_py.cpp."""
import unittest
import numpy as np
from xmsgrid_py.ugrid import XmUGrid
from xmsmesh.meshing import MeQuadBlossom

def list_int2d_to_cell_stream(cells):
    new_cells = []
    for cell in cells:
        if len(cell) == 3:
            new_cells += [5, 3]
            new_cells += cell
        elif len(cell) == 4:
            new_cells += [9, 4]
            new_cells += cell
    return new_cells

def make_XmUGrid(points, faces):
    return XmUGrid(points, list_int2d_to_cell_stream(faces))

class TestMeQuadBlossom(unittest.TestCase):
    """Test MeQuadBlossom functions."""

    @staticmethod
    def array_to_vec_pt3d(a_array):
        return [(a_array[i], a_array[i+1], 0) for i in range(0, len(a_array), 2)]

    @staticmethod
    def merge(points, indices):
        return [(points[i]) for i in indices]

    @staticmethod
    def assertArraysEqual(base, out):
        np.testing.assert_array_equal(np.array(base), out)

    def test_simple_triangle(self):
        points = [
            [0.0, 0.0, 0.0],  [10.0, 0.0, 0.0],  [20.0, 0.0, 0.0], [30.0, 0.0, 0.0],
            [0.0, 10.0, 0.0], [10.0, 10.0, 0.0], [20.0, 10.0, 0.0],
            [0.0, 20.0, 0.0], [10.0, 20.0, 0.0],
            [0.0, 30.0, 0.0]]
        triangles = [
            5, 3, 0, 1, 4,
            5, 3, 1, 5, 4,
            5, 3, 1, 2, 5,
            5, 3, 2, 6, 5,
            5, 3, 2, 3, 6,
            5, 3, 4, 5, 7,
            5, 3, 5, 8, 7,
            5, 3, 5, 6, 8,
            5, 3, 7, 8, 9]
        ugrid = XmUGrid(points, triangles)
        blossom = MeQuadBlossom(ugrid)
        split_boundary_points = False
        use_angle = False
        new_ugrid = blossom.make_quads(split_boundary_points, use_angle)
        cells = new_ugrid.get_cellstream()
        new_locations = new_ugrid.get_locations()
        expected_locations = points
        expected_cells = [
            9, 4, 1, 5, 4, 0,   # 0
            9, 4, 2, 6, 5, 1,   # 1
            5, 3, 2, 3, 6,      # 2
            9, 4, 5, 8, 7, 4,   # 3
            5, 3, 5, 6, 8,      # 4
            5, 3, 7, 8, 9]      # 5
        self.assertArraysEqual(expected_cells, cells)
        np.testing.assert_allclose(expected_locations, new_locations, atol=1.0e-4)

        blossom = MeQuadBlossom(ugrid)
        split_boundary_points = True
        use_angle = False
        new_ugrid = blossom.make_quads(split_boundary_points, use_angle)
        cells = new_ugrid.get_cellstream()
        new_locations = new_ugrid.get_locations()
        expected_locations = points
        expected_locations.append([50.0/3.0, 20.0/3.0, 0.0])
        expected_cells = (
            9, 4, 1, 5, 4, 0,
            9, 4, 2, 10, 5, 1,
            9, 4, 2, 3, 6, 10,
            9, 4, 5, 8, 7, 4,
            9, 4, 6, 8, 5, 10,
            5, 3, 7, 8, 9)
        self.assertArraysEqual(expected_cells, cells)
        np.testing.assert_allclose(expected_locations, new_locations, atol=1.0e-4)

    def test_simple_quad(self):
        points = [
            [0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [20.0, 0.0, 0.0], [30.0, 0.0, 0.0],
            [0.0, 10.0, 0.0], [10.0, 10.0, 0.0], [20.0, 10.0, 0.0], [30.0, 10.0, 0.0],
            [0.0, 20.0, 0.0], [10.0, 20.0, 0.0], [20.0, 20.0, 0.0], [30.0, 20.0, 0.0],
            [0.0, 30.0, 0.0], [10.0, 30.0, 0.0], [20.0, 30.0, 0.0], [30.0, 30.0, 0.0]]
        triangles = [
            [0, 1, 4], [1, 5, 4], [1, 2, 5], [2, 6, 5], [2, 3, 6], [3, 7, 6],
            [4, 9, 8], [4, 5, 9], [5, 10, 9], [5, 6, 10], [6, 7, 11], [6, 11, 10],
            [8, 9, 12], [9, 13, 12], [9, 10, 13], [10, 14, 13], [10, 11, 14],
            [11, 15, 14]]

        ugrid = make_XmUGrid(points, triangles)
        blossom = MeQuadBlossom(ugrid)
        split_boundary_points = False
        use_angle = False
        new_ugrid = blossom.make_quads(split_boundary_points, use_angle)
        cells = new_ugrid.get_cellstream()
        actual_locations = new_ugrid.get_locations()
        expected_cells = [
          9, 4, 1, 5, 4, 0,     # 0
          9, 4, 2, 6, 5, 1,     # 1
          9, 4, 3, 7, 6, 2,     # 2
          9, 4, 4, 5, 9, 8,     # 3
          9, 4, 5, 6, 10, 9,    # 4
          9, 4, 6, 7, 11, 10,   # 5
          9, 4, 9, 13, 12, 8,   # 6
          9, 4, 10, 14, 13, 9,  # 7
          9, 4, 11, 15, 14, 10 # 8
        ]
        expected_locations = points
        self.assertArraysEqual(expected_cells, cells)
        np.testing.assert_allclose(expected_locations, actual_locations, atol=1.0e-5)


        blossom = MeQuadBlossom(ugrid)
        spilt_boundary_points = True
        use_angle = False
        new_ugrid = blossom.make_quads(spilt_boundary_points, use_angle)
        cells = new_ugrid.get_cellstream()
        actual_locations = new_ugrid.get_locations()
    
        expected_cells = [
          9, 4, 1, 5, 4, 0,     # 0
          9, 4, 2, 6, 5, 1,     # 1
          9, 4, 3, 7, 6, 2,     # 2
          9, 4, 4, 5, 9, 8,     # 3
          9, 4, 5, 6, 10, 9,    # 4
          9, 4, 6, 7, 11, 10,   # 5
          9, 4, 9, 13, 12, 8,   # 6
          9, 4, 10, 14, 13, 9,  # 7
          9, 4, 11, 15, 14, 10, # 8
        ]
        expected_locations = points
        self.assertArraysEqual(expected_cells, cells)
        np.testing.assert_allclose(expected_locations, actual_locations, atol=1.0e-5)

    def test_complex_quad(self):
        points = [
            [-10.0, 0.0, 0.0], [0.0, 0.0, 0.0], [10.0, 0.0,0.0],
            [-15.0, 10.0, 0.0], [-5.0, 10.0, 0.0], [5.0, 10.0, 0.0], [15.0, 10.0, 0.0],
            [-20.0, 20.0, 0.0], [-10.0, 20.0, 0.0], [0.0, 20.0, 0.0], [10.0, 20.0, 0.0], [20.0, 20.0, 0.0],
            [-25.0, 30.0, 0.0], [-15.0, 30.0, 0.0], [-5.0, 30.0, 0.0], [5.0, 30.0, 0.0], [15.0, 30.0, 0.0], [25.0, 30.0, 0.0],
            [-30.0, 40.0, 0.0], [-20.0, 40.0, 0.0], [-10.0, 40.0, 0.0], [0.0, 40.0, 0.0], [10.0, 40.0, 0.0], [20.0, 40.0, 0.0], [30, 40, 0.0]]
        triangles = [
            [0, 4, 3], [0, 1, 4], [1, 5, 4], [1, 2, 5], [2, 6, 5],
            [3, 8, 7], [3, 4, 8], [4, 9, 8], [4, 5, 9], [5, 10, 9], [5, 11, 10], [5, 6, 11],
            [7, 13, 12], [7, 8, 13], [8, 14, 13], [8, 9, 14], [9, 15, 14], [9, 10, 15], [10, 16, 15], [10, 11, 16], [11, 17, 16],
            [12, 19, 18], [12, 13, 19], [13, 20, 19], [13, 14, 20], [14, 21, 20], [14, 15, 21], [15, 22, 21], [15, 23, 22], [15, 16, 23], [16, 17, 23], [17, 24, 23]]

        ugrid = make_XmUGrid(points, triangles)
        blossom = MeQuadBlossom(ugrid)
        num_boundary_edges = blossom.pre_make_quads()
        self.assertEqual(16, num_boundary_edges)
        spilt_boundary_points = False
        use_angle = False
        new_ugrid = blossom.make_quads(spilt_boundary_points, use_angle)
        cells = new_ugrid.get_cellstream()
        expected_cells = [
            9, 4, 0, 1, 4, 3,     #  0
            9, 4, 1, 2, 5, 4,     #  1
            9, 4, 3, 4, 8, 7,     #  2
            9, 4, 4, 5, 9, 8,     #  3
            9, 4, 5, 11, 10, 9,   #  4
            9, 4, 5, 2, 6, 11,    #  5
            9, 4, 7, 8, 13, 12,   #  6
            9, 4, 8, 9, 14, 13,   #  7
            9, 4, 9, 10, 15, 14,  #  8
            9, 4, 11, 17, 16, 10, #  9
            9, 4, 12, 13, 19, 18, # 10
            9, 4, 13, 14, 20, 19, # 11
            9, 4, 14, 15, 21, 20, # 12
            9, 4, 15, 23, 22, 21, # 13
            9, 4, 15, 10, 16, 23, # 14
            9, 4, 17, 24, 23, 16  # 15
        ]
        actual_locations = new_ugrid.get_locations()
        expected_locations = points
        self.assertArraysEqual(expected_cells, cells)
        np.testing.assert_allclose(expected_locations, actual_locations, atol=1.0e-5)

        blossom = MeQuadBlossom(ugrid)
        num_boundary_edges = blossom.pre_make_quads()
        self.assertEqual(16, num_boundary_edges)
        spilt_boundary_points = False
        use_angle = True
        new_ugrid = blossom.make_quads(spilt_boundary_points, use_angle)
        cells = new_ugrid.get_cellstream()
        expected_cells = [
            5, 3, 0, 4, 3,        #  0
            9, 4, 1, 5, 4, 0,     #  1
            9, 4, 2, 6, 5, 1,     #  2
            5, 3, 3, 8, 7,        #  3
            9, 4, 4, 9, 8, 3,     #  4
            9, 4, 5, 10, 9, 4,    #  5
            9, 4, 5, 6, 11, 10,   #  6
            5, 3, 7, 13, 12,      #  7
            9, 4, 8, 14, 13, 7,   #  8
            9, 4, 9, 15, 14, 8,   #  9
            9, 4, 10, 16, 15, 9,  # 10
            9, 4, 11, 17, 16, 10, # 11
            5, 3, 12, 19, 18,     # 12
            9, 4, 13, 20, 19, 12, # 13
            9, 4, 14, 21, 20, 13, # 14
            9, 4, 15, 22, 21, 14, # 15
            9, 4, 15, 16, 23, 22, # 16
            9, 4, 17, 24, 23, 16  # 17
        ]
        actual_locations = new_ugrid.get_locations()
        expected_locations = points
        self.assertArraysEqual(expected_cells, cells)
        np.testing.assert_allclose(expected_locations, actual_locations, atol=1.0e-5)

        blossom = MeQuadBlossom(ugrid)
        num_boundary_edges = blossom.pre_make_quads()
        self.assertEqual(16, num_boundary_edges)
        spilt_boundary_points = True
        use_angle = True
        new_ugrid = blossom.make_quads(spilt_boundary_points, use_angle)
        cells = new_ugrid.get_cellstream()
        expected_cells = (
            9, 4, 3, 0, 4, 25,
            9, 4, 1, 5, 4, 0,
            9, 4, 2, 6, 5, 1,
            9, 4, 8, 7, 3, 25,
            9, 4, 4, 9, 8, 25,
            9, 4, 5, 10, 9, 4,
            9, 4, 5, 6, 11, 10,
            9, 4, 12, 7, 13, 26,
            9, 4, 8, 14, 13, 7,
            9, 4, 9, 15, 14, 8,
            9, 4, 10, 16, 15, 9,
            9, 4, 11, 17, 16, 10,
            9, 4, 19, 18, 12, 26,
            9, 4, 13, 20, 19, 26,
            9, 4, 14, 21, 20, 13,
            9, 4, 15, 22, 21, 14,
            9, 4, 15, 16, 23, 22,
            9, 4, 17, 24, 23, 16)
        actual_locations = new_ugrid.get_locations()
        expected_locations = points
        expected_locations.append([-10, 13.33333, 0.0])
        expected_locations.append([-20.0, 33.33333, 0.0])
        self.assertArraysEqual(expected_cells, cells)
        np.testing.assert_allclose(expected_locations, actual_locations, atol=1.0e-5)

    def test_pre_make_quads(self):
        # single triangle
        #  2-----1
        #  |    /
        #  |   /
        #  |  /
        #  | /
        #  0)
        triangles = [[0, 1, 2]]
        points = [[0.0, 0.0, 0.0], [10.0, 10.0, 0.0], [0.0, 10.0, 0.0]]

        ugrid = make_XmUGrid(points, triangles)
        blossom = MeQuadBlossom(ugrid)
        num_boundary_edges = blossom.pre_make_quads()
        self.assertEqual(3, num_boundary_edges)
        split_boundary_points = True
        use_angle = False
        new_ugrid = blossom.make_quads(split_boundary_points, use_angle)
        expected_cells = [5, 3, 0, 1, 2]
        expected_locations = points
        actual_locations = new_ugrid.get_locations()
        self.assertArraysEqual(expected_cells, new_ugrid.get_cellstream())
        np.testing.assert_allclose(expected_locations, actual_locations, atol=1.0e-5)

        # two adjacent
        #  2-----3
        #  |    /|
        #  |   / |
        #  |  /  |
        #  | /   |
        #  0-----1
        triangles = [[2, 0, 3], [3, 0, 1]]
        points = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [10.0, 10.0, 0.0]]

        ugrid = make_XmUGrid(points, triangles)
        blossom = MeQuadBlossom(ugrid)
        num_boundary_edges = blossom.pre_make_quads()
        self.assertEqual(4, num_boundary_edges)

        new_ugrid = blossom.make_quads(split_boundary_points, use_angle)
        expected_cells = [9, 4, 0, 1, 3, 2]
        expected_locations = points
        actual_locations = new_ugrid.get_locations()
        self.assertArraysEqual(expected_cells, new_ugrid.get_cellstream())
        np.testing.assert_allclose(expected_locations, actual_locations, atol=1.0e-5)

        # three adjacent
        triangles = [[2, 9, 6], [9, 4, 6], [9, 3, 4]]
        points = [[0.0, 0.0, 0.0] for i in range(0, 10)]
        points[2] = [0.0, 0.0, 0.0]
        points[3] = [10.0, 20.0, 0.0]
        points[4] = [0.0, 20.0, 0.0]
        points[6] = [0.0, 10.0, 0.0]
        points[9] = [10.0, 10.0, 0.0]
    
        ugrid = make_XmUGrid(points, triangles)
        blossom = MeQuadBlossom(ugrid)
        num_boundary_edges = blossom.pre_make_quads()
        self.assertEqual(5, num_boundary_edges)
        split_boundary_points = True
        use_angle = False
        new_ugrid = blossom.make_quads(split_boundary_points, use_angle)
        expected_cells = [5, 3, 2, 9, 6, 9, 4, 4, 6, 9, 3]
        expected_locations = points
        actual_locations = new_ugrid.get_locations()
        self.assertArraysEqual(expected_cells, new_ugrid.get_cellstream())
        np.testing.assert_allclose(expected_locations, actual_locations, atol=1.0e-5)

    def test_estimated_run_time(self):
      self.assertAlmostEqual(0.0085, MeQuadBlossom.estimated_run_time_in_minutes(500), places=3)
      self.assertAlmostEqual(0.0680, MeQuadBlossom.estimated_run_time_in_minutes(1000), places=3)
      self.assertAlmostEqual(8.5000, MeQuadBlossom.estimated_run_time_in_minutes(5000), places=3)
      self.assertAlmostEqual(68.0, MeQuadBlossom.estimated_run_time_in_minutes(10000), places=3)
      self.assertAlmostEqual(68000.0, MeQuadBlossom.estimated_run_time_in_minutes(100000), places=3)

    def test_split_quads(self):
        #  8 ---- 9 ----10 ---- 11 
        #  |   \  |      |      |
        #  4 ---- 5 ---- 6 ---- 7
        #  |   /  |  \   |      |
        #  0 ---- 1 ---- 2 ---- 3
        points = [
            [ 0.0, 0.0, 0.0 ], [ 10.0, 0.0, 0.0 ], [ 20.0, 0.0, 0.0 ], [ 30.0, 0.0, 0.0 ],
            [0.0, 10.0, 0.0],  [10.0, 10.0, 0.0], [20.0, 10.0, 0.0], [30.0, 10.0, 0.0],
            [0.0, 20.0, 0.0], [10.0, 20.0, 0.0], [20.0, 20.0, 0.0], [30.0, 20.0, 0.0]]

        cells = [
            # row 1
            5, 3, 0, 5, 4,
            5, 3, 0, 1, 5,
            5, 3, 1, 2, 5,
            5, 3, 2, 6, 5,
            9, 4, 2, 3, 7, 6,
            # row 2
            5, 3, 4, 5, 8,
            7, 3, 5, 9, 8,
            9, 4, 5, 6, 10, 9,
            7, 4, 6, 7, 11, 10
        ]

        ugrid = XmUGrid(points, cells)
        # blossom = MeQuadBlossom(ugrid)
        new_ugrid = MeQuadBlossom.split_to_quads(ugrid)
        actual_locations = new_ugrid.get_locations()
        actual_cells = new_ugrid.get_cellstream()
       
        expected_locations = [
            [ 0.0, 0.0, 0.0 ], [ 10.0, 0.0, 0.0 ], [ 20.0, 0.0, 0.0 ], [ 30.0, 0.0, 0.0 ],
            [0.0, 10.0, 0.0],  [10.0, 10.0, 0.0], [20.0, 10.0, 0.0], [30.0, 10.0, 0.0],
            [0.0, 20.0, 0.0], [10.0, 20.0, 0.0], [20.0, 20.0, 0.0], [30.0, 20.0, 0.0],
            # The split edges
            [5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [5.0, 5.0, 0.0], [15.0, 0.0, 0.0], [10.0, 5.0, 0.0], # 12-16
            [25.0, 0.0, 0.0], [15.0, 5.0, 0.0], [20.0, 5.0, 0.0], [30.0, 5.0, 0.0], [5.0, 10.0, 0.0], # 17-21
            [0.0, 15.0, 0.0], [15.0, 10.0, 0.0], [5.0, 15.0, 0.0], [10.0, 15.0, 0.0], [25.0, 10.0, 0.0], # 22-26
            [20.0, 15.0, 0.0], [30.0, 15.0, 0.0], [5.0, 20.0, 0.0], [15.0, 20.0, 0.0], [25.0, 20.0, 0.0], # 27-31
             # The centroids
            [3.33333, 6.66667, 0.0], [6.66667, 3.33333, 0.0], # 32-33
            [13.33333, 3.33333, 0.0], [16.66667, 6.66667, 0.0], # 34-35
            [25.0, 5.0, 0.0], [3.33333, 13.33333, 0.0], [6.66667, 16.66667, 0.0], # 36-38
            [15.0, 15.0, 0.0], [25.0, 15.0, 0.0]  #39-40
        ]
        expected_cells = [
            9, 4, 0, 14, 32, 13,
            9, 4, 5, 21, 32, 14,
            9, 4, 4, 13, 32, 21,
            9, 4, 0, 12, 33, 14,
            9, 4, 1, 16, 33, 12,
            9, 4, 5, 14, 33, 16,
            9, 4, 1, 15, 34, 16,
            9, 4, 2, 18, 34, 15,
            9, 4, 5, 16, 34, 18,
            9, 4, 2, 19, 35, 18,
            9, 4, 6, 23, 35, 19,
            9, 4, 5, 18, 35, 23,
            9, 4, 2, 17, 36, 19,
            9, 4, 3, 20, 36, 17,
            9, 4, 7, 26, 36, 20,
            9, 4, 6, 19, 36, 26,
            9, 4, 4, 21, 37, 22,
            9, 4, 5, 24, 37, 21,
            9, 4, 8, 22, 37, 24,
            9, 4, 5, 25, 38, 24,
            9, 4, 9, 29, 38, 25,
            9, 4, 8, 24, 38, 29,
            9, 4, 5, 23, 39, 25,
            9, 4, 6, 27, 39, 23,
            9, 4, 10, 30, 39, 27,
            9, 4, 9, 25, 39, 30,
            9, 4, 6, 26, 40, 27,
            9, 4, 7, 28, 40, 26,
            9, 4, 11, 31, 40, 28,
            9, 4, 10, 27, 40, 31
        ]
  
        self.assertArraysEqual(expected_cells, actual_cells)
        np.testing.assert_allclose(expected_locations, actual_locations, atol=1.0e-5)



  
