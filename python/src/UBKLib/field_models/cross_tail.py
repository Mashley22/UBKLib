import numpy as np
from ..types import Vectorizable

CROSS_TAIL_STRENGTH = 10
SURFACE_POTENTIAL = 94.2


def stagnation_potential(
        cross_tail_strength: float = CROSS_TAIL_STRENGTH,
        surface_potential: float = SURFACE_POTENTIAL
) -> float:
    """
    Calculates the electric potential at the stagnation point for a uniform 
    cross-tail field combined with corotation.
    
    The stagnation point is located on the dawn flank at:
        x = 0
        y = -sqrt(surface_potential / cross_tail_strength)
    
    Parameters:
        cross_tail_strength: Cross-tail electric field component (kV/Re). 
                             Default is 10 kV/Re.
        surface_potential: Corotation potential at Earth's surface (kV). 
                           Default is 94.2 kV.

    Returns:
        Stagnation potential in kV.
    """
    
    # (dawn flank, y < 0)
    y_stag = -np.sqrt(surface_potential / cross_tail_strength)
    
    # At the stagnation point, x = 0
    r_stag = np.abs(y_stag)
    
    corotation = -surface_potential / r_stag
    convection = cross_tail_strength * y_stag
    
    return corotation + convection


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
