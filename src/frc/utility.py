"""General useful utilities for use cases related to the library."""

# Copyright (C) 2021                Department of Imaging Physics
# All rights reserved               Faculty of Applied Sciences
#                                   TU Delft
# Tip ten Brink

from frc.deps_types import Img, dip, np, Union, Func1D, NoIntersectionException
import frc._internal_utility as _util
from scipy.signal import windows as wins

__all__ = ['gaussf', 'square_image', 'tukey_square', 'apply_tukey', 'intersection']


def square_image(img: Img, add_padding=True) -> Img:
    """
    Transforms a 2D input image into a square image, using zero-padding by
    default and preserving the shortest axis.
    The padding places the input right/bottom of center. If trimming, top/
    left portion of the image is preserved.
    :param img: Array or DIP image that will be transformed.
    :param add_padding: If True, add zero-padding to make the image square, otherwise trim the image
    :return: Array or DIP image (based on input).
    """
    is_arr = isinstance(img, np.ndarray)
    if not is_arr:
        img = np.array(img)

    if img.ndim == 2:
        x_len = img.shape[0]
        y_len = img.shape[1]
        lengths = [x_len, y_len]
        min_axis = np.argmin(lengths)
        min_len = lengths[min_axis]
        if not add_padding:
            return _util.dipornp(img[:min_len, :min_len], is_arr, True)
        else:
            max_axis = 1 - min_axis
            max_len = lengths[max_axis]
            new_arr = np.zeros((max_len, max_len), dtype=img.dtype)
            ax_dif = (max_len - min_len)
            min_ax_start = int(np.ceil(ax_dif / 2.))
            min_ax_end = max_len - min_ax_start + (0 if ax_dif % 2 == 0 else 1)
            x_slice_start = 0 if max_axis == 0 else min_ax_start
            x_slice_end = x_len if max_axis == 0 else min_ax_end
            y_slice_start = 0 if max_axis == 1 else min_ax_start
            y_slice_end = y_len if max_axis == 1 else min_ax_end
            new_arr[x_slice_start:x_slice_end, y_slice_start:y_slice_end] = img
            return _util.dipornp(new_arr, is_arr, True)
    else:
        raise ValueError("Image not 2D.")


def gaussf(img: Img, sigmas=1) -> Img:
    """
    Gaussian filter using dip.Derivative.
    :param img: np.ndarray or dip.Image, input image
    :param sigmas: Gaussian kernel width, determines smoothing effect.
    :return: Array or DIP Image representing filtered input.
    """
    is_arr = isinstance(img, np.ndarray)

    if is_arr:
        img = dip.Image(img)
    # Gaussian filter
    img = dip.Derivative(img, derivativeOrder=0, sigmas=[sigmas])

    return _util.dipornp(img, is_arr, False)


def tukey_square(length, alpha=0.125) -> np.ndarray:
    """ Creates a square Tukey window mask of length x length of width
    alpha. """
    tukey = wins.tukey(length, alpha=alpha)
    return np.ones((length, length)) * tukey.reshape((length, 1)) * tukey.reshape((1, length))


def apply_tukey(img: np.ndarray, alpha=0.125) -> np.ndarray:
    """ Apply a Tukey window to a 2D, square input image. """
    main_size = img.shape[0]
    if img.ndim == 2 and main_size == np.ma.size(img, axis=1):
        return tukey_square(main_size, alpha) * img
    else:
        raise ValueError("Image not 2D or not square.")


def intersection(xs: np.ndarray, ys1_f: Union[np.ndarray, Func1D], ys2_f: Func1D, tol=0.08):
    """
    Calculates an intersection point between a series of (x,y) points ys1
    and a float function ys2 defined for xs.

    The algorithm is not designed to be commutative. For an accurate y-
    value the algorithm assumes uniform xs, which is also recommended.

    The algorithm is not designed to be commutative. Its general workings
    are as follows:

    It first calculates points for both functions for all x values. It then
    defines a tolerance region around the second function (y2) based on the
    maximums of both functions (8%). To prevent strange results, the
    functions should not diverge within the input range. It calculates a
    series of points that are this tolerance above and below y2.

    It then looks for indices in y1 that are within this region ("close").
    It also calculates all the crossings, i.e. the y1 indices where the
    next value will be below (if y1 started above y2) or above (if y1
    started below y2) y2. These crossings do not take the threshold into
    account.

    The crossings, the indices after each crossing as well as the close
    values are then looked at together as "points of interest". They are
    then separated into contiguous groups (i.e indexes following each
    other).

    Each group (a series of indexes) is evaluated independently. The y1 and
    y2 values of each group are evaluated and the difference computed. If
    there are both positive and negative values (i.e. there is a crossing),
    it is seen as a valid crossing group and no more groups are evaluated.

    It then calculates the average x-value of the group and returns that as
    the intersection point. It also returns the y-value (by inputting the
    x-value in the y2 function).

    :param xs: Array of x values for the y functions. If ys1_f is a series
    of points, xs must correspond to those points.
    :param ys1_f: Either a series of points or a float function.
    :param ys2_f: A 1D array function (ndarray -> ndarray).
    :param tol: Optional parameter that determines the tolerance for
    the region of interest.
    :return: x and y value of the intersection point, if it can be found.
    """

    if callable(ys1_f):
        ys1: np.ndarray = ys2_f(xs)
    else:
        ys1 = ys1_f

    ys2: np.ndarray = ys2_f(xs)
    max_ys1 = np.amax(ys1)
    max_ys2 = np.amax(ys2)
    full_max = max(max_ys1, max_ys2)
    tolerance = tol * full_max

    hi = ys2 + tolerance
    lo = ys2 - tolerance

    close_indices = np.asarray((ys1 < hi) & (ys1 > lo)).nonzero()[0]
    next_y = np.diff(ys1)
    thres_y = (ys2 - ys1)[:-1]

    # artificially increase first value to ensure proper crossings calculated
    if ys1[0] == ys2[0]:
        ys1[0] += tolerance * 0.5

    if ys1[0] > ys2[0]:
        # down
        crossings = np.asarray((next_y < thres_y) & (thres_y < 0)).nonzero()[0]
    else:
        # up
        crossings = np.asarray((next_y > thres_y) & (thres_y > 0)).nonzero()[0]

    crossings = np.concatenate((crossings, crossings + 1))

    close_indices = np.concatenate((close_indices, crossings))
    close_indices = np.unique(close_indices)
    index_groups = _util.separate_noncontiguous(close_indices)

    closest = []
    y_closest = []
    for group in index_groups:
        group_ys1 = ys1[group]
        group_ys2 = ys2[group]
        group_diff: np.ndarray = group_ys1 - group_ys2

        if np.any(group_diff >= 0) and np.any(group_diff <= 0):
            # An accurate y-value assumes uniform distribution of xs
            x_val = np.average(xs[group])
            y_closest.append(ys2_f(x_val))
            closest.append(x_val)
            break

    if closest:
        return closest[0], y_closest[0]
    else:
        raise NoIntersectionException("Could not find intersection!")
