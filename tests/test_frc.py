import pytest
import frc
import numpy as np


def test_simple_no_err():
    a = np.ones((5, 5))*20
    a[3, 2] = 16
    a[4, 3] = 18
    frc.one_frc(a)


def test_half_bit_threshold():
    from frc.frc_functions import half_bit_threshold

    threshold = half_bit_threshold(80)
    values = threshold(np.array([0.1, 0.2]))

    assert values.shape == (2,)
    assert np.all((values > 0) & (values < 1))


def test_frc_res_half_bit():
    rng = np.random.default_rng(0)
    img = rng.poisson(50, (64, 64)).astype(np.float64)
    img = frc.util.apply_tukey(img)
    frc_curve = frc.one_frc(img)
    xs = np.arange(len(frc_curve)) / img.shape[0]

    try:
        frc.frc_res(xs, frc_curve, img.shape[0], threshold='half_bit')
    except frc.deps_types.NoIntersectionException:
        pass


def pattern_image(n, mean, seed):
    y, x = np.mgrid[0:n, 0:n]
    waves = 0.4 * np.cos(2 * np.pi * (3 * x + 5 * y) / n) + 0.2 * np.cos(2 * np.pi * (17 * x - 9 * y) / n)
    return np.random.default_rng(seed).poisson(mean * (1 + waves))


def test_one_frc_sees_the_pattern():
    curve = frc.one_frc(pattern_image(128, 100.0, 1), seed=3)
    assert curve.shape == (64,) and curve.dtype == np.float64
    assert curve[5] > 0.9 and curve[19] > 0.9
    assert np.median(np.abs(curve[30:])) < 0.15


def test_seeds_repeat_curves():
    img = pattern_image(64, 30.0, 2)
    for method in (1, 2):
        assert np.array_equal(frc.one_frc(img, method, seed=5), frc.one_frc(img, method, seed=5))
        assert not np.array_equal(frc.one_frc(img, method, seed=5), frc.one_frc(img, method, seed=6))


def test_float_images_are_truncated_to_counts():
    img = pattern_image(32, 20.0, 3)
    assert np.array_equal(frc.one_frc(img + 0.5, seed=1), frc.one_frc(img, seed=1))


def test_two_frc_of_an_image_with_itself_is_one():
    img = frc.util.apply_tukey(pattern_image(32, 50.0, 4).astype(np.float64))
    assert np.allclose(frc.two_frc(img, img), 1.0)
    with pytest.raises(ValueError):
        frc.two_frc(img, img[:16, :16])


def test_half_bit_threshold_counts_the_rings_of_twice_the_size():
    from frc.frc_functions import half_bit_threshold

    # The threshold at ring 0 has one frequency: (s + f + 1) / (s + 1 + f), with s = sqrt(2)/2 - 1/2.
    assert half_bit_threshold(16)(0.0) == pytest.approx(1.0)
    # DIPlib's RadialSum counted 8 frequencies in ring 1.
    s = 0.5 * np.sqrt(2) - 0.5
    f = np.sqrt(s) * 2
    want = (s * np.sqrt(8) + (f + 1)) / ((s + 1) * np.sqrt(8) + f)
    assert half_bit_threshold(16)(1 / 16) == pytest.approx(want)
