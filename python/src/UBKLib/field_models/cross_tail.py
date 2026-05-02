import numpy as np
from typing import Union

CROSS_TAIL_STRENGTH = 10
SURFACE_POTENTIAL = 94.2

def cross_tail_potential(
    x: Union[float, np.ndarray], 
    y: Union[float, np.ndarray], 
) -> Union[float, np.ndarray]:
    """
    Calculates the equatorial electric potential for a uniform cross-tail field 
    combined with corotation for a field line with equitorial crossing (x, y)
    
    Parameters:
        x: Sunward-pointing coordinate in Earth Radii (Re).
        y: Duskward-pointing coordinate in Earth Radii (Re).

    Raises:
        ValueError: If the coordinates supplied are smaller than a dist of 1 Re

    Returns:
        Total electric potential in kV.
    """
    
    r = np.hypot(x, y)

    if np.any(r < 1.0):
        raise ValueError("Coordinates must supply a radial distance larger than 1 Re")
    
    corotation = - SURFACE_POTENTIAL / r
    convection = - CROSS_TAIL_STRENGTH * y
    
    return corotation + convection
