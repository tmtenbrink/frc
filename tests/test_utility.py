import pytest
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




@pytest.mark.parametrize("length, alpha, want", [
    (10, 0.125, [0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0]),
    (11, 0.5, [0.0, 0.34549150281252633, 0.9045084971874737, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9045084971874735, 0.3454915028125262, 0.0]),
    (8, 1.0, [0.0, 0.18825509907063326, 0.6112604669781572, 0.9504844339512095, 0.9504844339512095, 0.6112604669781572, 0.18825509907063326, 0.0]),
    (5, 0.0, [1.0, 1.0, 1.0, 1.0, 1.0]),
    (1, 0.3, [1.0]),
])
def test_tukey_window_is_scipys(length, alpha, want):
    # scipy.signal.windows.tukey(length, alpha), as frc used it before 0.3.
    assert np.allclose(frc.util.tukey_window(length, alpha), want, rtol=0, atol=1e-15)


def test_tukey_square_is_the_outer_product():
    w = frc.util.tukey_window(64)
    assert np.allclose(frc.util.tukey_square(64), np.outer(w, w))


def test_gaussf_is_diplibs():
    img = np.zeros((9, 9))
    img[4, 4] = 1
    img[0, 2] = 3
    # diplib.Derivative(img, derivativeOrder=0, sigmas=[1]), as frc used it before 0.3.
    want_first_rows = [
        [0.12930565709, 0.295060902398, 0.477723377072, 0.289753875056, 0.064652828545, 0.005307027342, 0.0, 0.0, 0.0],
        [0.078427845499, 0.178983135703, 0.289993284405, 0.176817567257, 0.040982931864, 0.00429183306, 0.000239409349, 1.9651916e-05, 0.0],
    ]
    got = frc.util.gaussf(img, 1)
    assert np.allclose(got[:2], want_first_rows, rtol=0, atol=1e-11)
    assert got[4, 4] == pytest.approx(0.159241125691, abs=1e-11)
