import numpy as np
from ..types import Vectorizable

# (nT Re^3)
SURFACE_STRENGTH = 31200

def equatorial_dipole_amplitude(
        x: Vectorizable,
        y: Vectorizable,
        surface_strength: float = SURFACE_STRENGTH,
    ) -> Vectorizable:
    """
    Calculates the value of the magnetic field strength at the (magnetic) equator,
    i.e. z = 0 in solar magnetic coordinates
    
    Parameters:
        x: x solar magnetic component (Re).
        y: y solar magnetic component (Re).

    Raises:
        ValueError: If the coordinates supplied are smaller than a dist of 1 Re

    Returns:
        The magnetic field amplitude (nT)
    """
    
    r = np.hypot(x, y)

    if np.any(r < 1.0):
        raise ValueError("Coordinates must supply a radial distance larger than 1 Re")
    
    return surface_strength / np.pow(r, 3)
