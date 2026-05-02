from collections.abc import Callable
from typing import Union
import numpy.typing as npt
import numpy as np

Vectorizable = Union[float, npt.NDArray[np.float64]]

PotentialFunction = Callable[
    [
        Vectorizable,
        Vectorizable
    ],
    Vectorizable
]
