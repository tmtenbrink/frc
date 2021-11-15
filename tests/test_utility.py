import frc
import numpy as np


def test_square():
    a = np.ones((5, 6)) * 20
    a[2, 3] = 19
    square_a: np.ndarray = frc.utility.square_image(a)
    assert square_a.shape == (6, 6)
    # Test moves to bottom of center
    assert square_a[3, 3] == 19


def test_square_right_center():
    a = np.ones((4, 3)) * 20
    a[2, 1] = 19
    square_a: np.ndarray = frc.utility.square_image(a)
    assert square_a.shape == (4, 4)
    # Test moves to right of center
    assert square_a[2, 2] == 19


