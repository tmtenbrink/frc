"""Shared dependencies and types module for use by other modules. Only
includes dependencies used in two or more related modules.
"""

# Copyright (C) 2021                Department of Imaging Physics
# All rights reserved               Faculty of Applied Sciences
#                                   TU Delft
# Tip ten Brink

import typing
import diplib
import numpy

__all__ = ['Img', 'Func1D', 'NoIntersectionException']

# Imports are redefined to prevent warnings that imports are unused
dip = diplib
np = numpy

Callable = typing.Callable
Union = typing.Union
Optional = typing.Optional
List = typing.List
Tuple = typing.Tuple

# Img type that allows interchangeable use of diplib Image and numpy array
Img = typing.TypeVar('Img', dip.Image, np.ndarray)

# One-dimensional function that maps (an array of) floats to (an array of) floats
Func1D = Callable[[Union[np.ndarray, float]], Union[np.ndarray, float]]


class NoIntersectionException(ValueError):
    """ Intersection could not be found, input did not take correct values for
    intersection calculation. """
    pass
