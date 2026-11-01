"""Functions to compute Fourier Ring Correlation of images.
"""

# Copyright (C) 2021                Department of Imaging Physics
# All rights reserved               Faculty of Applied Sciences
#                                   TU Delft
# Tip ten Brink

import rustfrc
from frc.deps_types import Img, np, Optional, Func1D, Union, Tuple
import frc._internal_utility as _util
import frc.utility as frc_util

__all__ = ['one_frc', 'two_frc', 'frc_res']


def two_frc(img1: Img, img2: Img) -> np.ndarray:
    """
    Standard FRC using two input images. Input images must be square, have
    equal dimensions and should have uncorrelated noise.

    The FRC at ring r, the frequencies f with r <= |f| < r + 1, is
    Re(sum F1 conj F2) / sqrt(sum |F1|^2 sum |F2|^2), with F1 and F2 the
    Fourier transforms of the images. See:
    Van Heel, M. (1987). Similarity measures between images.
    Ultramicroscopy, 21(1), 95–100. doi:10.1016/0304-3991(87)90010-6

    NOTE: The input images should be windowed to prevent FFT artifacts, i.e.
    using a 1/8 Tukey window.

    :param img1: Img representing the first image (half data set)
    :param img2: Img representing the second image (half data set)
    :return: A 1D array with the FRC value for each Fourier ring, with the
    index representing the pixel distance from the origin.
    """
    img1 = np.asarray(img1)
    img2 = np.asarray(img2)
    if not _util.are_sizes_equal([img1, img2]):
        raise ValueError(
            "Images are not of exact equal dimension, or axes are not equal! Trim or pad the images first.")

    return np.asarray(rustfrc.two_frc(img1, img2))


def one_frc(img: Img, method: int = 1, seed: Optional[int] = None) -> np.ndarray:
    """
        FRC using a single input image. The image is split into two by
        sampling binomial distributions for each pixel value. The input image
        must be a square 2D image of counts; floats are cast to int32, which
        truncates them.

        NOTE: The input image should be windowed to prevent FFT artifacts, i.e.
        using a 1/8 Tukey window.

        By default, method 1 performs one split (hence, method '1'), with the
        second image being the difference between the original image and the
        binomially sampled image.
        Method 2 performs two independent splits (hence, method '2').
        Method 1 is preferred.

        :param Img img: Img representing the image (full data set)
        :param int method: Int representing which method to use. The default is
        method 1.
        :param seed: A seed for the split, from 0 to 2**64 - 1. The same image
        and seed give the same curve. Without one, every call splits with a
        fresh seed from the operating system.
        :returns: A 1D array with the FRC value for each Fourier ring, with the
        index representing the pixel distance from the origin.
        """
    img = np.asarray(img)
    if not _util.are_axes_equal(img):
        raise ValueError("Image does not have equal axes! Trim or pad it first.")
    if img.ndim != 2:
        raise ValueError("Image is not 2D.")
    if not np.issubdtype(img.dtype, np.integer):
        img = img.astype(np.int32)

    if method == 1:
        return np.asarray(rustfrc.one_frc(img, seed))
    elif method == 2:
        seeds = (None, None) if seed is None else np.random.SeedSequence(seed).generate_state(2, np.uint64).tolist()
        half1 = rustfrc.binom_split(img, seeds[0])
        half2 = rustfrc.binom_split(img, seeds[1])
        return np.asarray(rustfrc.two_frc(half1, half2))
    else:
        raise ValueError("Choose either method 1 or 2 for binomial split.")


def half_bit_threshold(img_size: int) -> Func1D:
    """
    1/2-bit information threshold curve for Fourier Ring Correlation.
    Assumes a square image of image_size x image_size pixels.

    Van Heel, M., & Schatz, M. (2005). Fourier shell correlation threshold
    criteria. Journal of Structural Biology, 151(3), 250–262.
    doi:10.1016/j.jsb.2005.05.009

    :param img_size:
    :return: A function representing the 1/2-bit information curve for a
    specific square size, where an input of 1 for x yields the threshold
    for the ring that is img_size pixels away from the origin and zero for
    the ring at the origin.
    """
    # Number of pixels in a ring at distance 'array index' from the origin
    # up to a distance 'img_size'
    nr = np.asarray(rustfrc.ring_counts(img_size * 2))[:img_size].astype(np.float64)
    nr_rt = np.sqrt(nr)

    # Exact calculations
    snr_half_set = 0.5 * np.sqrt(2) - 0.5
    factor = np.sqrt(snr_half_set) * 2
    half_bit_thres = (snr_half_set * nr_rt + (factor + 1)) / ((snr_half_set + 1) * nr_rt + factor)

    def half_bit_thres_f(x):
        """ Maps a number x (0 to 1) to threshold value for FRC value for pixel
        value (0 to img_size) """
        indexes = np.array(x * img_size).astype(int)
        return half_bit_thres[indexes]

    return half_bit_thres_f


def frc_res(xs: np.ndarray, frc_curve: np.ndarray, img_size: int,
            threshold: Optional[Union[str, Func1D]] = None) -> Tuple[float, float, Func1D]:
    """
    Calculate the resolution from a FRC curve.

    For threshold, choose '1/7' (default) for the 1/7 threshold value or
    'half_bit' for the 1/2-bit information curve threshold.

    :param xs: Array of FRC curve x-values.
    :param frc_curve: Array of FRC y-values corresponding to xs.
    :param img_size: Length of the sides of the square input image
    :param threshold: String representing which threshold curve to use
    to calculate the resolution or actual Func1D threshold function.
    :return: A tuple of FRC resolution as 1/intersection-x (float),
    intersection-y (float) and the exact threshold Func1D used.
    """
    if threshold is None or threshold == '1/7':
        def threshold_f(x):
            return x*0 + 1./7
    elif threshold == 'half_bit':
        threshold_f = half_bit_threshold(img_size)
    elif callable(threshold):
        threshold_f = threshold
    else:
        raise ValueError("Unknown threshold!")

    # Calculate intersection using FRC utility module
    intersect_x, intersect_y = frc_util.intersection(xs, frc_curve, threshold_f)

    return 1. / intersect_x, intersect_y, threshold_f
