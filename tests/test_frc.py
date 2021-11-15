import frc
import numpy as np


def test_simple_no_err():
    a = np.ones((5, 5))*20
    a[3, 2] = 16
    a[4, 3] = 18
    frc.one_frc(a)
