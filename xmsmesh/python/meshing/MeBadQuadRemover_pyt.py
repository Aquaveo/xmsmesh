"""Test BadQuadRemover_py.cpp."""
import unittest
import numpy as np
from xmsmesh.meshing import BadQuadRemover
from xmsgrid_py.ugrid import XmUGrid

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

class TestBadQuadRemover(unittest.TestCase):
    """Test BadQuadRemover functions."""

    def setUp(self):
        pass

    @staticmethod
    def array_to_vec_pt3d(a_array):
        return [(a_array[i], a_array[i+1], 0) for i in range(0, len(a_array), 2)]

    @staticmethod
    def assertArraysEqual(base, out):
        np.testing.assert_array_equal(np.array(base), out)

    def test_collapse(self):
        """Testing quad removal on a more complex quad UGrid. The original mesh
        has two rows of quads with some 3 diamond shaped quads.  One if on the left
        boundary.  One is roughly in the middle.  One is near the right side.  There
        are also two pairs of quads where each pair is separated by 2 sequential
        shared edges.  One has an aspect ratio to collapse.  The other does not.
        Then toward the right is a quad that has 3 points with 3 edges attached. It
        should not collapse.  One the bottom right is a quad with 4 3 edge points,
        and it also should not collapse."""

        points = []
        for y in range(20, -10, -10):
            for x in range(0, 100, 10):
                points.append([x,y])

        points[10][0] -= 3
        points[11][0] -= 3
        points[12][0] -= 3
        points[13][0] -= 3

        points[14][0] += 3
        points[15][0] += 1.5

        points[17][0] -= 1.5
        points[18][0] -= 3
        points[19][0] += 5

        points.append([3, 10])    #30
        points.append([33, 10])   #31
        points.append([82, 10])    #32
        points.append([90, 12])  #33
        points.append([90, 8])  #34
        points.append([15, 15])  #35
        points.append([15, 5])  #36

        points.append([100, 10])  #37
        points.append([70, -10])  #38
        points.append([80, -10])  #39
        points.append([73, -3])  #40
        points.append([76, -4])  #41
        points.append([77, -7])  #42
        points.append([74, -6])  #43

        faces = [[10, 20, 30, 0],  [30, 11, 1, 0],   [35, 12, 2, 1],   [11, 12, 35, 1],
                        [12, 13, 3, 2],   [31, 14, 4, 3],   [14, 15, 5, 4],   [15, 16, 6, 5],
                        [16, 17, 7, 6],   [17, 18, 8, 7],   [32, 33, 9, 8],   [20, 21, 11, 30],
                        [36, 22, 12, 11], [21, 22, 36, 11], [22, 23, 13, 12], [13, 23, 31, 3],
                        [23, 24, 14, 31], [24, 25, 15, 14], [25, 26, 16, 15], [26, 27, 17, 16],
                        [27, 28, 18, 17], [18, 28, 32, 8],  [28, 29, 34, 32], [34, 19, 33, 32],
                        [33, 19, 37, 9],  [34, 29, 37, 19], [40, 41, 28, 27], [41, 42, 39, 28],
                        [42, 43, 38, 39], [38, 43, 40, 27], [43, 42, 41, 40]]

        remover = BadQuadRemover(make_XmUGrid(points, faces))
        collapsed_ugrid = remover.remove_bad_quads(0.7)
        expected_locations = [
            # top row (0-9)
            [0, 20, 0],
            [10, 20, 0],
            [20, 20, 0],
            [30, 20, 0],
            [40, 20, 0],
            [50, 20, 0],
            [60, 20, 0],
            [70, 20, 0],
            [80, 20, 0],
            [90, 20, 0],
            # row 2 (10-18)
            [-3, 10, 0],
            [7, 10, 0],
            [17, 10, 0],
            [30, 10, 0],
            [43, 10, 0],
            [51.5, 10, 0],
            [60, 10, 0],
            [68.5, 10, 0],
            [95, 10, 0],
            # bottom row (19-28)
            [0, 0, 0],
            [10, 0, 0],
            [20, 0, 0],
            [30, 0, 0],
            [40, 0, 0],
            [50, 0, 0],
            [60, 0, 0],
            [70, 0, 0],
            [80, 0, 0],
            [90, 0, 0],
            # extra points in the middle row
            [79.8571, 10, 0],
            [90, 12, 0],
            [90, 8, 0],
            [15, 15, 0],
            [100, 10, 0],
            # extra points at the bottom
            [70, -10, 0],
            [80, -10, 0],
            [73, -3, 0],
            [76, -4, 0],
            [77, -7, 0],
            [74, -6, 0]]
        actual_locations = collapsed_ugrid.get_locations()
        import numpy as np
        np.testing.assert_allclose(expected_locations, actual_locations, atol=1.0e-4)

        expected_cells = [
        9, 4, 10, 11, 1,  0,  9, 4, 32, 12, 2,  1,  9, 4, 11, 12, 32, 1,  9, 4, 12, 13, 3,  2,
        9, 4, 13, 14, 4,  3,  9, 4, 14, 15, 5,  4,  9, 4, 15, 16, 6,  5,  9, 4, 16, 17, 7,  6,
        9, 4, 17, 29, 8,  7,  9, 4, 29, 30, 9,  8,  9, 4, 19, 20, 11, 10, 9, 4, 20, 21, 12, 11,
        9, 4, 21, 22, 13, 12, 9, 4, 22, 23, 14, 13, 9, 4, 23, 24, 15, 14, 9, 4, 24, 25, 16, 15,
        9, 4, 25, 26, 17, 16, 9, 4, 26, 27, 29, 17, 9, 4, 27, 28, 31, 29, 9, 4, 31, 18, 30, 29,
        9, 4, 30, 18, 33, 9,  9, 4, 31, 28, 33, 18, 9, 4, 36, 37, 27, 26, 9, 4, 37, 38, 35, 27,
        9, 4, 38, 39, 34, 35, 9, 4, 34, 39, 36, 26, 9, 4, 39, 38, 37, 36]
        actual_cells = collapsed_ugrid.get_cellstream()
        self.assertArraysEqual(expected_cells, actual_cells)

    def testCollapseQuadTri(self):
        """Test simple mesh with one quad and one triangle that share two
        adjacent edges. Expected to result in one triangle."""

        faces = [[0, 2, 1, 3], [0, 1, 2]]

        points = [
            [-10, 0, 0],
            [10, 0, 0], # row 1
            [0, 10, 0], # row 2
            [0, 20, 0]  # row 3
            ]

        remover = BadQuadRemover(make_XmUGrid(points, faces))
        collapsed_ugrid = remover.remove_bad_quads(0.7)
        expected_locations = [[-10, 0, 0], [10, 0, 0], [0, 20, 0]]

        actual_locations = collapsed_ugrid.get_locations()
        import numpy as np
        np.testing.assert_allclose(expected_locations, actual_locations, atol=1.0e-4)

        expected_cells = [5, 3, 0, 1, 2]

        actual_cells = collapsed_ugrid.get_cellstream()
        self.assertArraysEqual(expected_cells, actual_cells)
