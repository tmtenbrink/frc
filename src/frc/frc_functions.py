"""Functions to compute Fourier Ring Correlation of images.
"""

# Copyright (C) 2021                Department of Imaging Physics
# All rights reserved               Faculty of Applied Sciences
#                                   TU Delft
# Tip ten Brink

import rustfrc as r_frc
from frc.deps_types import Img, np, dip, Callable, Optional, Func1D, Union, Tuple
import frc._internal_utility as _util
import frc.utility as frc_util

__all__ = ['one_frc', 'two_frc', 'frc_res']


def _frc(img1: dip.Image, img2: dip.Image) -> np.ndarray:
    """
    Internal function to perform Fourier Ring Correlation on a 2D dip.Image.
    Arguments img1 and img2 are assumed to be of the same dimensions and
    square. Do not use this function directly, but use two_frc or one_frc.

    See below paper for FRC formula.
    Van Heel, M. (1987). Similarity measures between images.
    Ultramicroscopy, 21(1), 95–100. doi:10.1016/0304-3991(87)90010-6

    :param img1: dip.Image representing the first image.
    :param img2: dip.Image representing the second image.
    :return: A 1D array with the FRC value for Fourier ring, with the index
    representing the pixel distance from the origin.
    """
    axis_size = img1.Size(0)

    # Fourier transformation of each image
    fourier1 = dip.FourierTransform(img1, set(""))
    fourier2 = dip.FourierTransform(img2, set(""))

    # Numerator of the FRC formula
    frc_num = np.real(dip.RadialSum(fourier1 * dip.Conjugate(fourier2), None))
    frc_denom1 = dip.RadialSum(r_frc.sqr_abs(np.array(fourier1).astype(np.complex64)), None)
    frc_denom2 = dip.RadialSum(r_frc.sqr_abs(np.array(fourier2).astype(np.complex64)), None)
    frc_denom_prod = frc_denom1 * frc_denom2
    frc_denom = np.sqrt(frc_denom_prod)
    frc = frc_num / frc_denom

    # Find last sensible FRC value (beyond axis_size / 2 there are no full
    # rings)
    circle_edge = int(axis_size / 2)

    return frc[:circle_edge]


def two_frc(img1: Img, img2: Img) -> np.ndarray:
    """
    Standard FRC using two input images, which can be either DIP images or
    NumPy arrays. Input images must be square, have equal dimensions and
    should have uncorrelated noise.

    NOTE: The input images should be windowed to prevent FFT artifacts, i.e.
    using a 1/8 Tukey window.

    :param img1: Img representing the first image (half data set)
    :param img2: Img representing the second image (half data set)
    :return: A 1D array with the FRC value for Fourier ring, with the index
    representing the pixel distance from the origin.
    """
    # Check if array to test if images must be converted
    is_arr = isinstance(img1, np.ndarray)

    if is_arr:
        dip1 = dip.Image(img1)
        dip2 = dip.Image(img2)
    else:
        dip1 = img1
        dip2 = img2
        img1 = np.array(img1)
        img2 = np.array(img2)

    if not _util.are_sizes_equal([img1, img2]):
        raise ValueError(
            "Images are not of exact equal dimension, or axes are not equal! Trim or pad the images first.")

    return _frc(dip1, dip2)


def one_frc(img: Img, method: int = 1) -> np.ndarray:
    """
        FRC using a single input image, which can be either DIP images or NumPy
        arrays. The image is split into two by sampling binomial distributions
        for each pixel value. The input image must be a square 2D image.

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
        :returns: A 1D array with the FRC value for Fourier ring, with the index
        representing the pixel distance from the origin.
        """
    # Check if array to test if images must be converted
    is_arr = isinstance(img, np.ndarray)

    # Converting is necessary for the binomial split
    if not is_arr:
        img = np.array(img)

    if not _util.are_axes_equal(img):
        raise ValueError("Image does not have equal axes! Trim or pad it first.")

    split_method: Callable[[np.ndarray], Tuple[np.ndarray, np.ndarray]]
    if method == 1:
        split_method = _one_frc1
    elif method == 2:
        split_method = _one_frc2
    else:
        raise ValueError("Choose either method 1 or 2 for binomial split.")

    # Split is performed
    img_half1, img_half2 = split_method(img)

    # _frc function requires DIP Image inputs
    dip1 = dip.Image(img_half1)
    dip2 = dip.Image(img_half2)

    return _frc(dip1, dip2)


def _one_frc2(img: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """ Performs two independent binomial splits. """
    return r_frc.binom_split(img), r_frc.binom_split(img)


def _one_frc1(img: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """ Performs a single binomial split, with the second image equal to the
    difference of the sampled image and the input image. """
    img_half1 = r_frc.binom_split(img)
    return img_half1, img - img_half1


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
    nr = np.array(dip.RadialSum(np.ones((img_size*2, img_size*2)), None))[:img_size]
    nr_rt = np.sqrt(nr)

    # Exact calculations
    snr_half_set = 0.5 * np.sqrt(2) - 0.5
    factor = np.sqrt(snr_half_set) * 2
    half_bit_thres = (snr_half_set * nr_rt + (factor + 1)) / ((snr_half_set + 1) * nr_rt + factor)

    def half_bit_thres_f(x):
        """ Maps a number x (0 to 1) to threshold value for FRC value for pixel
        value (0 to img_size) """
        indexes = np.array(x * img_size).astype(np.int)
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
