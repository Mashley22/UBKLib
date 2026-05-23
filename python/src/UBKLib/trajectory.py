from scipy.optimize import minimize_scalar
import numpy as np
from collections.abc import Callable

from .types import (
    W0ContourFunction,
    UBTrajectory,
    UBCoord,
    Vectorizable
)
from .field_models import dipole

__MINIZER_OPTIONS = {
    'xatol': 1e-16,    
    'maxiter': 10000    
}


def classical_ub_trajectory(
        B_m: Vectorizable,
        mu: float,
        charge: int,
        intercept: float = 0
) -> Vectorizable:
    """
    Calculates the classical trajectory of the mirror points of
    a particle in U-B space as a function U(B_m)

    Args:
        B_m (nT): the magnitude of the magnetic field at the mirror point
        mu (eV/nT): the value of the first adiabatic invariant
        charge (e): the charge of the particle
        intercept (kV): the intercept/constant term i.e U(0)

    Returns:
        the values U(B_m) for the given parameters
    """

    return intercept - B_m * mu / charge


def relativistic_ub_trajectory(
        B_m: Vectorizable,
        mu: float,
        charge: int,
        rest_mass: float,
        intercept: float = 0
) -> Vectorizable:
    """
    Calculates the relivistice trajectory of the mirror points of
    a particle in U-B space as a function U(B_m)

    Args:
        B_m (nT): the magnitude of the magnetic field at the mirror point
        mu (keV/nT): the value of the first adiabatic invariant
        charge (e): the charge of the particle
        rest_mass (keV): the rest mass of the particle
        intercept (kV): the intercept/constant term i.e U(0)

    Returns:
        the values U(B_m) for the given parameters
    """

    return intercept - rest_mass * np.sqrt(1 + 2 * mu * B_m / rest_mass) / charge + rest_mass / charge


def __continuous_lcds_ub_min_intercept(
        lower_contour: W0ContourFunction,
        upper_contour: W0ContourFunction,
        trajectory: Callable[[float], float],
        B_bounds: tuple[float, float]
) -> float:

    def lower_minize_func(x):
        return -trajectory(x) + lower_contour(x)

    def upper_minize_func(x):
        return -trajectory(x) + upper_contour(x)
    
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

    return max(lower_result.fun, upper_result.fun)


def __continuous_lcds_ub_max_intercept(
        lower_contour: W0ContourFunction,
        upper_contour: W0ContourFunction,
        trajectory: Callable[[float], float],
        B_bounds: tuple[float, float]
) -> float:

    def lower_minize_func(x):
        return trajectory(x) - lower_contour(x)

    def upper_minize_func(x):
        return trajectory(x) - upper_contour(x)
    
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

    return min(-lower_result.fun, -upper_result.fun)


def __continuous_ub_trajectory(
        lower_contour: W0ContourFunction,
        upper_contour: W0ContourFunction,
        trajectory: Callable[[float], float],
        intercept: float,
        B_bounds: tuple[float, float] = (dipole.SURFACE_STRENGTH / (15 ** 3), dipole.SURFACE_STRENGTH / (1.05 ** 3))
) -> UBTrajectory:
    """Calculates the the closed trajectory W=0 intercepts for a particle with a given trajectory, charge and intercept
    
    Args:
        lower_contour: a function that describes the lower W=0 contour as a function U(B)
        upper_contour: same as lower_contour but for the upper W=0 contour
        trajectory: the U(B_m) for a general particle with the same properties
        charge (e): the charge of the particle
        intercept (kV): the constant term in the trajectory
        B_bounds (nT): the bounds for the value of B to search within, by default chooses field delta as the equatorial strength
                      at distances 15 to 1.05 Re
    Raises:
        RuntimeError: if the scipy minizer fails

    Returns:
        The trajectory information
    """

    trajectoryInf = UBTrajectory()
    trajectoryInf.intercept = intercept

    lower_result = minimize_scalar(lambda x: (trajectory(x) + trajectoryInf.intercept - lower_contour(x)) ** 2,
                                   bounds=B_bounds,
                                   method='bounded',
                                   options=__MINIZER_OPTIONS)

    if lower_result.success is False:
        raise RuntimeError(f"Scipy minimization failed: {lower_result.message}")

    upper_result = minimize_scalar(lambda x: (trajectory(x) + trajectoryInf.intercept - upper_contour(x)) ** 2,
                                   bounds=B_bounds,
                                   method='bounded',
                                   options=__MINIZER_OPTIONS)

    if upper_result.success is False:
        raise RuntimeError(f"Scipy minimization failed: {upper_result.message}")

    trajectoryInf.lower_intercept = UBCoord(B=lower_result.x, U=lower_contour(lower_result.x))
    trajectoryInf.upper_intercept = UBCoord(B=upper_result.x, U=upper_contour(upper_result.x))

    return trajectoryInf


def continuous_lcds_ub_trajectory(
        lower_contour: W0ContourFunction,
        upper_contour: W0ContourFunction,
        trajectory: Callable[[float], float],
        B_bounds: tuple[float, float] = (dipole.SURFACE_STRENGTH / (15 ** 3), dipole.SURFACE_STRENGTH / (1.05 ** 3))
) -> UBTrajectory:
    """Calculates the the LCDS trajectory for a particle with a given trajectory, charge
    
    Args:
        lower_contour: a function that describes the lower W=0 contour as a function U(B)
        upper_contour: same as lower_contour but for the upper W=0 contour
        trajectory: the U(B_m) for a general particle with the same properties
        charge (e): the charge of the particle
        B_bounds (nT): the bounds for the value of B to search within, by default chooses field delta as the equatorial strength
                      at distances 15 to 1.05 Re
    Raises:
        RuntimeError: if the scipy minizer fails

    Returns:
        The trajectory information
    """
    minIntercept = __continuous_lcds_ub_min_intercept(lower_contour, upper_contour, trajectory, B_bounds)
    maxIntercept = __continuous_lcds_ub_max_intercept(lower_contour, upper_contour, trajectory, B_bounds)

    minTrajectory = __continuous_ub_trajectory(lower_contour, upper_contour, trajectory, minIntercept, B_bounds)
    maxTrajectory = __continuous_ub_trajectory(lower_contour, upper_contour, trajectory, maxIntercept, B_bounds)

    if minTrajectory.lower_intercept.B < maxTrajectory.lower_intercept.B:
        return minTrajectory

    return maxTrajectory
