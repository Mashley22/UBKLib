from collections.abc import Callable
from typing import List
import numpy as np
import numpy.typing as npt
from skimage import measure

from .types import PotentialFunction, MagneticAmplitudeFunction, TurningPoint, TurningPointType

def __generate_potential_image(
        potential : PotentialFunction,
        x_bounds: tuple[float, float],
        y_bounds: tuple[float, float],
        resolution : int) -> npt.NDArray[np.float64]:

    x_vals = np.linspace(x_bounds[0], x_bounds[1], resolution)
    y_vals = np.linspace(y_bounds[0], y_bounds[1], resolution)
    X_2d, Y_2d = np.meshgrid(x_vals, y_vals)
    
    X_1d = X_2d.flatten()
    Y_1d = Y_2d.flatten()

    R_1d = np.hypot(X_1d, Y_1d)

    valid_mask = (R_1d >= 1.0)

    U_1d = np.full_like(X_1d, np.nan, dtype=float)
    
    U_1d[valid_mask] = potential(X_1d[valid_mask], Y_1d[valid_mask])

    return U_1d.reshape(X_2d.shape)

def __lower_turning_point(
        turning_point: TurningPoint
    ) -> bool:

    return turning_point.y < 0

def __upper_turning_point(
        turning_point: TurningPoint
    ) -> bool:

    return turning_point.y > 0

def __parse_turning_points(
        turning_points: List[List[List[TurningPoint]]],
        cond: Callable[[TurningPoint], bool]
    ) -> List[TurningPoint]:

    retVal = []

    for i in range(len(turning_points)):
        for j in range(len(turning_points[i])):
            for point in turning_points[i][j]:
                if cond(point) is True:
                    retVal.append(point)

def find_contour_turning_points(
        contour: npt.NDArray[np.float64],
        magnetic_amplitude_func: MagneticAmplitudeFunction
    ) -> List[TurningPoint]:

    """Finds the turning points of the magnetic amplitude for a given contour
    and a given magnetic ampltidue function
    """
    
    turningPoints = []
    magnetic_amplitudes = [magnetic_amplitude_func(x[0], x[1]) for x in contour]

    for i in range(1, len(magnetic_amplitudes) - 1):
        prev = magnetic_amplitudes[i - 1]
        cur = magnetic_amplitudes[i]
        nxt = magnetic_amplitudes[i + 1]

        diff_prev = cur - prev
        diff_curr = nxt - cur

        if diff_prev >= 0 and diff_curr < 0:
            turningPoints.append(
                    TurningPoint(
                        type=TurningPointType.MAXIMUM,
                        x=contour[i][0],
                        y=contour[i][1],
                        B=cur)
            )

        elif diff_prev < 0 and diff_curr >= 0:
            turningPoints.append(
                    TurningPoint(
                        type=TurningPointType.MINIMUM,
                        x=contour[i][0],
                        y=contour[i][1],
                        B=cur)
            )

    return turningPoints

def find_all_turning_points(
        contours: List[List[npt.NDArray[np.float64]]],
        magnetic_amplitude_func: MagneticAmplitudeFunction
    ) -> List[List[List[TurningPoint]]]:

    """
    Finds all the turning points of all the contours, the return value structure is the same as
    the input contour structure.

    Returns:
        Returns a nested list of lists of lists of TurningPoint,
        follows the same structure as the contours, the first index is over the different potential values,
        the second over the contours for this potential, the last list contains the turning points in this contour
    """

    turningPoints = [[] for _ in range(len(contours))]

    for i in range(len(contours)):
        for contour in contours[i]:
             turningPoints[i].append(find_contour_turning_points(contour, magnetic_amplitude_func))

    return turningPoints

def get_lower_contour_turning_points(
        turning_points: List[List[List[TurningPoint]]]
    ) -> List[TurningPoint]:
    """Retrieves the lower contour turning points from the output of find_all_turning_points
    returns them as a list
    """
    
    __parse_turning_points(turning_points, __lower_turning_point)

def get_upper_contour_turning_points(
        turning_points: List[List[List[TurningPoint]]]
    ) -> List[TurningPoint]:
    """Retrieves the upper contour turning points from the output of find_all_turning_points
    returns them as a list
    """
    
    __parse_turning_points(turning_points, __upper_turning_point)

def generate_equipotentials(
        potential : PotentialFunction,
        levels : List[float],
        x_bounds: tuple[float, float] = (-15.0, 15.0),
        y_bounds: tuple[float, float] = (-15.0, 15.0),
        resolution : int = 100
    ) -> List[List[npt.NDArray[np.float64]]]:

    """Extracts isolines (equipotentials) from a potential field.

    This function discretizes the potential into a 2D image and uses the 
    Marching Squares algorithm to find contours at specific potential levels.

    Args:
        potential: The physics model used to calculate the field.
        levels: A list of specific potential values to find contours for.
        x_bounds: Domain limits for the X-axis. Defaults to (-15.0, 15.0).
        y_bounds: Domain limits for the Y-axis. Defaults to (-15.0, 15.0).
        resolution: The grid density used for contour extraction.

    Returns:
        A nested list structure: `results[level_index][contour_index]`. 
        Each contour is an (N, 2) array, where the first index is over 
        the points in each contours, the second refers to x or y coords
        
    Raises:
        AssertionError: If bounds are provided where min >= max.
    """

    assert x_bounds[0] < x_bounds[1], "X-axis min must be less than max."
    assert y_bounds[0] < y_bounds[1], "Y-axis min must be less than max."

    potential_image = __generate_potential_image(potential, x_bounds, y_bounds, resolution)
    contours = []

    for i in range(len(levels)):
        raw_contours = measure.find_contours(potential_image, levels[i])
        mapped_contours = []
        
        for c in raw_contours:
            x = np.interp(c[:, 1], [0, resolution - 1], [x_bounds[0], x_bounds[1]])
            y = np.interp(c[:, 0], [0, resolution - 1], [y_bounds[0], y_bounds[1]])
            mapped_contours.append(np.column_stack([x, y]))
        
        contours.append(mapped_contours)

    return contours
