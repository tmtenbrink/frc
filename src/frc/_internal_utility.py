"""Internal utilities for use by other library modules."""

# Copyright (C) 2021                Department of Imaging Physics
# All rights reserved               Faculty of Applied Sciences
#                                   TU Delft
# Tip ten Brink

from frc.deps_types import Img, np, dip, List


def dipornp(img: Img, was_arr: bool, is_arr: bool) -> Img:
    if was_arr:
        if is_arr:
            return img
        else:
            return np.array(img)
    else:
        if not is_arr:
            return img
        else:
            return dip.Image(img)


def are_axes_equal(a: np.ndarray) -> bool:
    shape = a.shape
    shape_arr = np.array(shape)
    first_value = np.full(shape_arr.shape, shape[0])
    axes_equal = np.all(shape_arr == first_value)
    return axes_equal


def are_sizes_equal(arrs: List[np.ndarray], check_equal_axes=True) -> bool:
    sizes = np.empty((len(arrs)), dtype=tuple)
    first_value = sizes.copy()
    for i, a in enumerate(arrs):
        sizes[i] = a.shape

    first_value.fill(sizes[0])

    equal_sizes = np.all(sizes == first_value)

    return equal_sizes and (not check_equal_axes or are_axes_equal(arrs[0]))


def separate_noncontiguous(a: np.ndarray) -> List[np.ndarray]:
    diff = np.diff(a)
    # Indices of a that start a new noncontiguous group
    diff_i_gt1 = np.asarray(diff > 1).nonzero()[0] + 1

    return np.split(a, diff_i_gt1)
