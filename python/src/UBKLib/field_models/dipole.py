import numpy as np
from ..types import Vectorizable

# (nT Re^3)
SURFACE_STRENGTH = 31200

def equitorial_dipole_amplitude(
        x: Vectorizable,
        y: Vectorizable,
        surface_strength: float = SURFACE_STRENGTH,
    ) -> Vectorizable:
    """
    Calculates the value of the magnetic field strength at the equator
    
    Parameters:
        x: Sunward-pointing coordinate in Earth Radii (Re).
        y: Duskward-pointing coordinate in Earth Radii (Re).

    Raises:
        ValueError: If the coordinates supplied are smaller than a dist of 1 Re

    Returns:
        The magnetic field amplitude (nT)
    """
    
    r = np.hypot(x, y)

    if np.any(r < 1.0):
        raise ValueError("Coordinates must supply a radial distance larger than 1 Re")
    
    return surface_strength / np.pow(r, 3)
