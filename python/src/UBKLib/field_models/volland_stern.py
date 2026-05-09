import numpy as np
from ..types import Vectorizable

SURFACE_POTENTIAL = 94.2
CONVECTION_NUMERATOR = 0.045
CONVECTION_KP_LINEAR_COEFF = -0.0159
CONVECTION_KP_QUADRATIC_COEFF = 0.0093


def stagnation_potential(
        kp: float = 3.0,                 
        surface_potential: float = SURFACE_POTENTIAL,
        convection_numerator: float = CONVECTION_NUMERATOR
) -> float:
    """
    Calculate the stagnation potential of the Volland-Stern electric field model.
    
    The stagnation point occurs at dusk (x = 0, y < 0) where the corotation
    and convection electric fields balance, resulting in zero E×B drift.
    
    Parameters:
        kp: Kp index for geomagnetic activity.
    
    Returns:
        Stagnation potential in kV.
    """
    e0 = CONVECTION_NUMERATOR / (
        1.0 + CONVECTION_KP_LINEAR_COEFF * kp + CONVECTION_KP_QUADRATIC_COEFF * kp**2
    )**3
    
    stagnation_distance = (SURFACE_POTENTIAL / (2 * e0))**(1 / 3)
    
    return volland_stern_potential(0, stagnation_distance, kp, surface_potential, convection_numerator)


def volland_stern_potential(
        x: Vectorizable, 
        y: Vectorizable, 
        kp: float = 3.0, 
        surface_potential: float = SURFACE_POTENTIAL,
        convection_numerator: float = CONVECTION_NUMERATOR
) -> Vectorizable:
    """
    Calculates the Volland Stern electric potential along a field line with equatorial
    crossing at (x,  y) (solar magnetic)
    
    Parameters:
        x: x solar magnetic component (Re).
        y: y solar magnetic component (Re).
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
