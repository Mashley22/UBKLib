from collections.abc import Callable
from typing import Union
import numpy.typing as npt
import numpy as np
from dataclasses import dataclass
from enum import Enum

Vectorizable = Union[float, npt.NDArray[np.float64]]

PotentialFunction = Callable[
    [
        Vectorizable,
        Vectorizable
    ],
    Vectorizable
]

MagneticAmplitudeFunction = Callable[
    [
        float,
        float
    ],
    float
]

W0ContourFunction = Callable[
    [
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
        gradient(kV/nT):  the value of -mu/q
        intercept(kV): the intercept of this given trajectory in the UB space
        lower_intercept: the UB coordinates of the interception point of this trajectory
                         with the lower W=0 contour, is None if the curve does not intercept this contour.
        upper_intercept: the same as lower_intercept but for the upper W=0 contour.
    """
    gradient: float = 0
    intercept: float = 0
    lower_intercept: Optional[UBCoord] = None
    upper_intercept: Optional[UBCoord] = None

class TurningPointType(Enum):
    MAXIMUM = 0
    MINIMUM = 1

@dataclass
class TurningPoint:
    type : TurningPointType
    x : float
    y : float
    B : float
    U : float = 0
