import numpy as np
from ..types import Vectorizable

CROSS_TAIL_STRENGTH = 10
SURFACE_POTENTIAL = 94.2

def cross_tail_potential(
        x: Vectorizable,
        y: Vectorizable,
        cross_tail_strength: float = CROSS_TAIL_STRENGTH,
        surface_potential: float = SURFACE_POTENTIAL
    ) -> Vectorizable:
    """
    Calculates the equatorial electric potential for a uniform cross-tail field 
    combined with corotation for a field line with equatorial crossing (x, y) (solar magnetic)
    
    Parameters:
        x: x solar magnetic component (Re).
        y: y solar magnetic component (Re).

    Raises:
        ValueError: If the coordinates supplied are smaller than a dist of 1 Re

    Returns:
        Total electric potential in kV.
    """
    
    r = np.hypot(x, y)

    if np.any(r < 1.0):
        raise ValueError("Coordinates must supply a radial distance larger than 1 Re")
    
    corotation = - SURFACE_POTENTIAL / r
    convection = CROSS_TAIL_STRENGTH * y
    
    return corotation + convection
