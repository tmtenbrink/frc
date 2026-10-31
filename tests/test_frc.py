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
