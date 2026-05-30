from collections.abc import Callable
from typing import Union, Optional
import numpy.typing as npt
import numpy as np
from dataclasses import dataclass
from enum import Enum
from typing import List

Vectorizable = Union[float, npt.NDArray[np.float64]]

"""
inputs are x, y in solar magnetic coordinates
(k or z information should be baked into the function directly)
"""
PotentialFunction = Callable[
    [
        Vectorizable,
        Vectorizable
    ],
    Vectorizable
]

"""
inputs are x, y in solar magnetic coordinates
(k or z information should be baked into the function directly)
"""
MagneticAmplitudeFunction = Callable[
    [
        float,
        float
    ],
    float
]

"""
inputs are x, y in solar magnetic coordinates AND a list of k values
the return value is the magnetic intensity at each k value
"""
MagneticAmplitudeFunctionWithK = Callable[
    [
        float,
        float,
        List[float]
    ],
    List[float]
]

W0ContourFunction = Callable[
    [
        Vectorizable
    ],
    Vectorizable
]

"""
First input is the array of magnetic field intensity,
second is the electric potential.
"""
HamiltonianFunction = Callable[
    [
        Vectorizable,
        Vectorizable
    ],
    Vectorizable
]


@dataclass
class UBCoord:
    B: float
    U: float
 

@dataclass
class UBTrajectory:
    """
    Stores the information about a particle trajectory in UB coordinates for a given 
    value of K. Typically calculated from the :func:'findContinuousUBTrajectory' where
    gradient is supplied as input, hence the gradient is only for completeness.

    Attributes:
        intercept (kV): the intercept of this given trajectory in the UB space
        lower_intercept: the UB coordinates of the interception point of this trajectory
                         with the lower W=0 contour, is None if the curve does not intercept this contour.
        upper_intercept: the same as lower_intercept but for the upper W=0 contour.
    """
    intercept: float = 0
    lower_intercept: Optional[UBCoord] = None
    upper_intercept: Optional[UBCoord] = None


class W0ContourType(Enum):
    """
    Defines the type of the W = 0 contour found, whether its the lower contour, the upper contour, 
    or some random one.
    """
    UPPER = 0
    LOWER = 1
    NONE = 2


class TurningPointType(Enum):
    MAXIMUM = 0
    MINIMUM = 1


@dataclass
class TurningPoint:
    type: TurningPointType
    x: float
    y: float
    B: float
    U: float = 0


@dataclass
class TurningPointWithErrors:
    pos: TurningPoint
    err: TurningPoint
