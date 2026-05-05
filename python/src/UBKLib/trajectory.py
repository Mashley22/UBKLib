import numpy as np
from scipy.optimize import minimize_scalar
from dataclasses import dataclass
from typing import Optional

from .types import W0ContourFunction, Vectorizable, UBTrajectory, UBCoord
from .field_models import dipole

__MINIZER_OPTIONS = {
    'xatol': 1e-16,    
    'maxiter': 10000    
}

def __continuous_lcds_ub_intercept(
        lower_contour: W0ContourFunction,
        upper_contour: W0ContourFunction,
        gradient: float,
        B_bounds: tuple[float, float]
    ) -> float:

    if gradient >= 0:
        lower_minize_func = lambda x: gradient * x - lower_contour(x)
        upper_minize_func = lambda x: gradient * x - upper_contour(x)
    else:
        lower_minize_func = lambda x: -(gradient * x - lower_contour(x))
        upper_minize_func = lambda x: -(gradient * x - upper_contour(x))
    
    lower_result = minimize_scalar(lower_minize_func,
                                   bounds=B_bounds,
                                   method='bounded',
                                   options=__MINIZER_OPTIONS)

    if lower_result.success is False:
        raise RuntimeError(f"Scipy minimization failed: {lower_result.message}")

    upper_result = minimize_scalar(upper_minize_func,
                                   bounds=B_bounds,
                                   method='bounded',
                                   options=__MINIZER_OPTIONS)

    if upper_result.success is False:
        raise RuntimeError(f"Scipy minimization failed: {upper_result.message}")
        
    if gradient >= 0:
        return min(-lower_result.fun, -upper_result.fun)
    else:
        return max(lower_result.fun, upper_result.fun)

def continuous_lcds_ub_trajectory(
        lower_contour: W0ContourFunction,
        upper_contour: W0ContourFunction,
        gradient: float,
        B_bounds: tuple[float, float] = (dipole.SURFACE_STRENGTH / (15 **3), dipole.SURFACE_STRENGTH / (1.05 ** 3))
    ) -> UBTrajectory:
    """Calculates the the LCDS trajectory for a particle with a given gradient.
    
    Args:
        lower_contour: a function that describes the lower W=0 contour as a function U(B)
        upper_contour: same as lower_contour but for the upper W=0 contour
        gradient(kV/nT): the value of the -mu/q for the given particle
        B_bounds(nT): the bounds for the value of B to search within, by default chooses field delta as the equatorial strength
                      at distances 15 to 1.05 Re
    Raises:
        RuntimeError: if the scipy minizer fails

    Returns:
        The trajectory information
    """

    trajectory = UBTrajectory()
    trajectory.gradient = gradient
    trajectory.intercept = __continuous_lcds_ub_intercept(lower_contour, upper_contour, gradient, B_bounds)


    lower_result = minimize_scalar(lambda x: (gradient * x + trajectory.intercept - lower_contour(x)) ** 2,
                                   bounds=B_bounds,
                                   method='bounded',
                                   options=__MINIZER_OPTIONS)

    if lower_result.success is False:
        raise RuntimeError(f"Scipy minimization failed: {lower_result.message}")

    upper_result = minimize_scalar(lambda x: (gradient * x + trajectory.intercept - upper_contour(x)) ** 2,
                                   bounds=B_bounds,
                                   method='bounded',
                                   options=__MINIZER_OPTIONS)

    if upper_result.success is False:
        raise RuntimeError(f"Scipy minimization failed: {upper_result.message}")

    trajectory.lower_intercept = UBCoord(B=lower_result.x, U=lower_contour(lower_result.x))
    trajectory.upper_intercept = UBCoord(B=upper_result.x, U=upper_contour(upper_result.x))

    return trajectory
