import numpy as np
from typing import Union

SURFACE_POTENTIAL = 94.2
CONVECTION_NUMERATOR = 0.045
CONVECTION_KP_LINEAR_COEFF = -0.0159
CONVECTION_KP_QUADRATIC_COEFF = 0.0093

def volland_stern_potential(
    x: Union[float, np.ndarray], 
    y: Union[float, np.ndarray], 
    kp: float = 3.0, 
) -> Union[float, np.ndarray]:
    """
    Calculates the Volland Stern electric potential along a field line with equitorial
    crossing at (x,  y)
    
    Parameters:
        x: Sunward-pointing coordinate in Earth Radii (Re).
        y: Duskward-pointing coordinate in Earth Radii (Re).
        kp: Kp index for geomagnetic activity.
        
    Raises:
        ValueError: If the coordinates supplied are not of a radial distance larger than 1 Re

    Returns:
        Total electric potential in kV.
    """
    r = np.hypot(x, y)

    if np.any(r < 1.0):
        raise ValueError("Coordinates must supply a radial distance larger than 1 Re")
    
    e0 = CONVECTION_NUMERATOR / (
        1.0 + CONVECTION_KP_LINEAR_COEFF * kp + CONVECTION_KP_QUADRATIC_COEFF * kp**2
    )**3
    
    corotation = -SURFACE_POTENTIAL / r
    convection = -e0 * r * y
    
    return corotation + convection
