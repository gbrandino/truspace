#!/usr/bin/env python3
import sys
import unittest
import numpy as np
import os

filename = 'Eigenvalues_p_3_q_4_minLev_0_maxLev_10_Per_2_c_0.09_Subspace_0_Mom_0'
filepath1 = os.getcwd() + '/../../test/' + filename
filepath2 = os.getcwd() + '/../../test/serial/' + filename


class TestTruspace(unittest.TestCase):

    def runTest(self):
        file1 = np.genfromtxt(filepath1, delimiter=' ', dtype=np.double)
        file2 = np.genfromtxt(filepath2, delimiter=' ', dtype=np.double)
        print(file1)
        print(file2)
        diff = np.abs(file1 - file2)
        equal = not np.any(diff > 1e-14)
        self.assertEqual(equal, True)


if __name__ == '__main__':
    unittest.TextTestRunner().run(TestTruspace())
