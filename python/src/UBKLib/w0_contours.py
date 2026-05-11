from collections.abc import Callable
import numpy as np
import numpy.typing as npt
from typing import List
from scipy.interpolate import CubicSpline
import pandas as pd

from .types import (
    TurningPoint,
    TurningPointType,
    MagneticAmplitudeFunction
)


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
                if cond(point):
                    retVal.append(point)

    return retVal


def single_contour_w0_points(
        contour: npt.NDArray[np.float64],
        magnetic_amplitude_func: MagneticAmplitudeFunction
) -> List[TurningPoint]:
    """Finds the w=0 points along a given contour of constant eletric potential,
    assumes that these points are also the points of dB/dS = 0 where S parametrises
    the contour of constant U
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
                    B=cur
                )
            )

        elif diff_prev < 0 and diff_curr >= 0:
            turningPoints.append(
                TurningPoint(
                    type=TurningPointType.MINIMUM,
                    x=contour[i][0],
                    y=contour[i][1],
                    B=cur
                )
            )

    return turningPoints


def contour_w0_points(
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
            turningPoints[i].append(single_contour_w0_points(contour, magnetic_amplitude_func))

    return turningPoints


def parse_lower_contour_w0_points(
        turning_points: List[List[List[TurningPoint]]]
) -> List[TurningPoint]:
    """Retrieves the lower contour (in U-B space) turning points from the output of find_all_turning_points
    returns them as a single list of TurningPoints
    """
    return __parse_turning_points(turning_points, __lower_turning_point)


def parse_upper_contour_w0_points(
        turning_points: List[List[List[TurningPoint]]]
) -> List[TurningPoint]:
    """Retrieves the upper contour (in U-B space) turning points from the output of find_all_turning_points
    returns them as a single list of TurningPoints
    """
    
    return __parse_turning_points(turning_points, __upper_turning_point)


def generate_ub_spline(
        turning_points: List[TurningPoint]
) -> CubicSpline:
    """Computes and returns the cubic spline for the turing points as a function U(B)
    """
    sorted_list = sorted(turning_points, key=lambda x: x.B)

    df = pd.DataFrame({'x': [x.B for x in sorted_list], 'y': [x.U for x in sorted_list]})
    grouped = df.groupby('x')['y'].mean().reset_index()

    return CubicSpline(grouped['x'].values, grouped['y'].values)


def generate_realSpace_splines(
        turning_points: List[TurningPoint]
) -> tuple[CubicSpline, CubicSpline]:
    """Computes and returns the cubic splines for the x(B) and y(B) respectively
    """

    sorted_list = sorted(turning_points, key=lambda x: x.B)

    retVal = []

    df = pd.DataFrame({'x': [x.B for x in sorted_list], 'y': [x.x for x in sorted_list]})
    grouped = df.groupby('x')['y'].mean().reset_index()

    retVal.append(CubicSpline(grouped['x'].values, grouped['y'].values))

    df = pd.DataFrame({'x': [x.B for x in sorted_list], 'y': [x.y for x in sorted_list]})
    grouped = df.groupby('x')['y'].mean().reset_index()

    retVal.append(CubicSpline(grouped['x'].values, grouped['y'].values))

    return retVal
