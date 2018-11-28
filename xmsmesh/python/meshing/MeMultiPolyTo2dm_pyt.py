"""Test MultiPolyTo2dm_py.cpp."""
import unittest
import numpy as np
import os
from xmscore_py.misc import Observer
from xmsmesh.meshing import MultiPolyTo2dm
from xmsmesh.meshing import MultiPolyMesherIo
from xmsmesh.meshing import PolyInput
from xmsmesh.meshing import RefinePoint
from xmsinterp_py.interpolate import InterpLinear
from xmsinterp_py.interpolate import InterpIdw


class TestMultiPolyTo2dm(unittest.TestCase):
    """Test MultiPolyTo2dm functions."""

    def setUp(self):
        self.maxDiff = None

    @staticmethod
    def array_to_vec_pt3d(a_array):
        return [(a_array[i], a_array[i+1], 0) for i in range(0, len(a_array), 2)]

    def test_creating_MultiPolyTo2dm(self):
        to2dm = MultiPolyTo2dm()
        self.assertIsInstance(to2dm, MultiPolyTo2dm)

    def test_generate_2dm(self):
        to2dm = MultiPolyTo2dm()
        io = MultiPolyMesherIo()
        rv = to2dm.generate_2dm(io, "fname.2dm")
        self.assertTrue(os.path.isfile("fname.2dm"))

    def test_case_4(self):
        # build test case 4 polys
        out_a = (0,  60, 10, 60, 20, 60, 30, 60, 30, 50, 30, 40, 30, 30, 30, 20, 30, 10,
                 30, 0,  20, 0,  10, 0,  0,  0,  0,  10, 0,  20, 0,  30, 0,  40, 0,  50)
        in_a1 = (10, 50, 10, 40, 20, 40, 20, 50)
        in_a2 = (10, 20, 10, 10, 20, 10, 20, 20)
        poly_input_a = PolyInput()
        poly_input_a.outside_poly = self.array_to_vec_pt3d(out_a)
        poly_input_a.inside_polys = (self.array_to_vec_pt3d(in_a1), self.array_to_vec_pt3d(in_a2))
        poly_input_a.bias = 1.0

        out_b = (30, 60, 40, 60, 50, 60, 60, 60, 70, 60, 70, 50, 70, 40, 70, 30,
                 70, 20, 70, 10, 70, 0,  60, 0,  50, 0,  40, 0,  40, 10, 30, 10,
                 30, 20, 40, 20, 40, 30, 30, 30, 30, 40, 40, 40, 40, 50, 30, 50)
        in_b1 = (50, 50, 50, 40, 60, 40, 60, 50)
        in_b2 = (50, 20, 50, 10, 60, 10, 60, 20)
        poly_input_b = PolyInput()
        poly_input_b.outside_poly = self.array_to_vec_pt3d(out_b)
        poly_input_b.inside_polys = (self.array_to_vec_pt3d(in_b1), self.array_to_vec_pt3d(in_b2))
        poly_input_b.bias = 1.0

        io = MultiPolyMesherIo()
        io.poly_inputs = (poly_input_a, poly_input_b)

        # mesh the polys
        to_2dm = MultiPolyTo2dm()
        (success, result) = to_2dm.generate_2dm(io, "", 8)
        
        expected = "MESH2D\n" \
                   "E3T     1     1    18     2     1\n" \
                   "E3T     2     2    18    24     1\n" \
                   "E3T     3     2    24     3     1\n" \
                   "E3T     4     3    23     4     1\n" \
                   "E3T     5     3    24    23     1\n" \
                   "E3T     6     4    23    27     1\n" \
                   "E3T     7     4    27     5     1\n" \
                   "E3T     8     5    20     6     1\n" \
                   "E3T     9     5    27    20     1\n" \
                   "E3T    10     6    19     7     1\n" \
                   "E3T    11     6    20    19     1\n" \
                   "E3T    12     7    19     8     1\n" \
                   "E3T    13     8    19    22     1\n" \
                   "E3T    14     8    22     9     1\n" \
                   "E3T    15     9    11    10     1\n" \
                   "E3T    16     9    22    11     1\n" \
                   "E3T    17    10    11    33     1\n" \
                   "E3T    18    10    33    34     1\n" \
                   "E3T    19    11    22    12     1\n" \
                   "E3T    20    12    13    31     1\n" \
                   "E3T    21    12    21    28     1\n" \
                   "E3T    22    12    22    21     1\n" \
                   "E3T    23    12    28    13     1\n" \
                   "E3T    24    12    31    32     1\n" \
                   "E3T    25    13    26    14     1\n" \
                   "E3T    26    13    28    26     1\n" \
                   "E3T    27    14    15    47     1\n" \
                   "E3T    28    14    26    15     1\n" \
                   "E3T    29    14    47    30     1\n" \
                   "E3T    30    15    25    16     1\n" \
                   "E3T    31    15    26    25     1\n" \
                   "E3T    32    16    25    17     1\n" \
                   "E3T    33    17    24    18     1\n" \
                   "E3T    34    17    25    24     1\n" \
                   "E3T    35    20    27    28     1\n" \
                   "E3T    36    20    28    21     1\n" \
                   "E3T    37    23    26    29     1\n" \
                   "E3T    38    23    29    27     1\n" \
                   "E3T    39    26    28    29     1\n" \
                   "E3T    40    27    29    28     1\n" \
                   "E3T    41    30    47    53     1\n" \
                   "E3T    42    30    52    31     1\n" \
                   "E3T    43    30    53    52     1\n" \
                   "E3T    44    31    52    56     1\n" \
                   "E3T    45    31    56    32     1\n" \
                   "E3T    46    32    49    33     1\n" \
                   "E3T    47    32    56    49     1\n" \
                   "E3T    48    33    48    34     1\n" \
                   "E3T    49    33    49    48     1\n" \
                   "E3T    50    34    48    35     1\n" \
                   "E3T    51    35    48    51     1\n" \
                   "E3T    52    35    51    36     1\n" \
                   "E3T    53    36    38    37     1\n" \
                   "E3T    54    36    51    38     1\n" \
                   "E3T    55    38    51    39     1\n" \
                   "E3T    56    39    50    57     1\n" \
                   "E3T    57    39    51    50     1\n" \
                   "E3T    58    39    57    40     1\n" \
                   "E3T    59    40    55    41     1\n" \
                   "E3T    60    40    57    55     1\n" \
                   "E3T    61    41    55    42     1\n" \
                   "E3T    62    42    54    43     1\n" \
                   "E3T    63    42    55    54     1\n" \
                   "E3T    64    43    54    44     1\n" \
                   "E3T    65    44    53    45     1\n" \
                   "E3T    66    44    54    53     1\n" \
                   "E3T    67    45    47    46     1\n" \
                   "E3T    68    45    53    47     1\n" \
                   "E3T    69    49    56    57     1\n" \
                   "E3T    70    49    57    50     1\n" \
                   "E3T    71    52    55    58     1\n" \
                   "E3T    72    52    58    56     1\n" \
                   "E3T    73    55    57    58     1\n" \
                   "E3T    74    56    58    57     1\n" \
                   "ND     1      0.0      0.0      0.0\n" \
                   "ND     2      0.0     10.0      0.0\n" \
                   "ND     3      0.0     20.0      0.0\n" \
                   "ND     4      0.0     30.0      0.0\n" \
                   "ND     5      0.0     40.0      0.0\n" \
                   "ND     6      0.0     50.0      0.0\n" \
                   "ND     7      0.0     60.0      0.0\n" \
                   "ND     8     10.0     60.0      0.0\n" \
                   "ND     9     20.0     60.0      0.0\n" \
                   "ND    10     30.0     60.0      0.0\n" \
                   "ND    11     30.0     50.0      0.0\n" \
                   "ND    12     30.0     40.0      0.0\n" \
                   "ND    13     30.0     30.0      0.0\n" \
                   "ND    14     30.0     20.0      0.0\n" \
                   "ND    15     30.0     10.0      0.0\n" \
                   "ND    16     30.0      0.0      0.0\n" \
                   "ND    17     20.0      0.0      0.0\n" \
                   "ND    18     10.0      0.0      0.0\n" \
                   "ND    19     10.0     50.0      0.0\n" \
                   "ND    20     10.0     40.0      0.0\n" \
                   "ND    21     20.0     40.0      0.0\n" \
                   "ND    22     20.0     50.0      0.0\n" \
                   "ND    23     10.0     20.0      0.0\n" \
                   "ND    24     10.0     10.0      0.0\n" \
                   "ND    25     20.0     10.0      0.0\n" \
                   "ND    26     20.0     20.0      0.0\n" \
                   "ND    27  8.51082 31.75796      0.0\n" \
                   "ND    28  19.9189 32.42692      0.0\n" \
                   "ND    29 14.90632 26.65924      0.0\n" \
                   "ND    30     40.0     20.0      0.0\n" \
                   "ND    31     40.0     30.0      0.0\n" \
                   "ND    32     40.0     40.0      0.0\n" \
                   "ND    33     40.0     50.0      0.0\n" \
                   "ND    34     40.0     60.0      0.0\n" \
                   "ND    35     50.0     60.0      0.0\n" \
                   "ND    36     60.0     60.0      0.0\n" \
                   "ND    37     70.0     60.0      0.0\n" \
                   "ND    38     70.0     50.0      0.0\n" \
                   "ND    39     70.0     40.0      0.0\n" \
                   "ND    40     70.0     30.0      0.0\n" \
                   "ND    41     70.0     20.0      0.0\n" \
                   "ND    42     70.0     10.0      0.0\n" \
                   "ND    43     70.0      0.0      0.0\n" \
                   "ND    44     60.0      0.0      0.0\n" \
                   "ND    45     50.0      0.0      0.0\n" \
                   "ND    46     40.0      0.0      0.0\n" \
                   "ND    47     40.0     10.0      0.0\n" \
                   "ND    48     50.0     50.0      0.0\n" \
                   "ND    49     50.0     40.0      0.0\n" \
                   "ND    50     60.0     40.0      0.0\n" \
                   "ND    51     60.0     50.0      0.0\n" \
                   "ND    52     50.0     20.0      0.0\n" \
                   "ND    53     50.0     10.0      0.0\n" \
                   "ND    54     60.0     10.0      0.0\n" \
                   "ND    55     60.0     20.0      0.0\n" \
                   "ND    56 48.37114 31.77082      0.0\n" \
                   "ND    57 59.82848 32.41625      0.0\n" \
                   "ND    58 54.80369 26.66427      0.0\n"
        self.assertTrue(success)
        self.assertEqual(expected, result)
