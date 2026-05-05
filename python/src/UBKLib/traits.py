from collections.abc import Callable
from typing import Union
import numpy.typing as npt
import numpy as np
from dataclasses import dataclass

Vectorizable = Union[float, npt.NDArray[np.float64]]

PotentialFunction = Callable[
    [
        Vectorizable,
        Vectorizable
    ],
    Vectorizable
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
