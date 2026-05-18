from collections.abc import Callable
import numpy as np
import numpy.typing as npt
from typing import List
from scipy.interpolate import (
    CubicSpline,
    make_splrep
)
import pandas as pd

from .types import (
    TurningPoint,
    TurningPointType,
    MagneticAmplitudeFunction,
    PotentialFunction
)

from .equipotentials import generate_equipotentials


SINGULAR_VAL = 1e-12


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


def __split_into_valid_segments(
        contour: npt.NDArray[np.float64],
        magnetic_amplitudes: MagneticAmplitudeFunction
) -> tuple[npt.NDArray[np.float64], list[float]]:
    
    contour_segments = []
    amplitude_segments = []

    ctr_segment = []
    amp_segment = []

    for amp, pos in zip(magnetic_amplitudes, contour):
        if np.isfinite(amp):
            amp_segment.append(amp)
            ctr_segment.append(pos)
        else:
            if len(amp_segment) > 0:
                amplitude_segments.append(amp_segment)
                contour_segments.append(ctr_segment)
                amp_segment = []
                ctr_segment = []

    if len(amp_segment) > 0:
        amplitude_segments.append(amp_segment)
        contour_segments.append(ctr_segment)
    
    return (contour_segments, amplitude_segments)


def __single_contour_w0_points(
        contour: npt.NDArray[np.float64],
        magnetic_amplitudes: MagneticAmplitudeFunction
) -> List[TurningPoint]:

    if len(contour) < 32:
        return []

    x_idxs = np.arange(len(contour))

    if len(magnetic_amplitudes) < 4:
        return []

    mag_amp_spl = make_splrep(x_idxs, magnetic_amplitudes, s=len(x_idxs))

    d_spl = mag_amp_spl.derivative(1)
    d2_spl = mag_amp_spl.derivative(2)
    
    x_dense = x_idxs
    deriv_vals = d_spl(x_dense)
    
    sign_changes = np.where(np.diff(np.sign(deriv_vals)))[0]

    turningPoints = []
    cidxs = []

    for idx in sign_changes:
        # only expecting maxs or mins
        x1, x2 = x_dense[idx], x_dense[idx + 1]
        if abs(deriv_vals[idx + 1] - deriv_vals[idx]) > SINGULAR_VAL:
            critical_idx = x1 - deriv_vals[idx] * (x2 - x1) / (deriv_vals[idx + 1] - deriv_vals[idx])
        else:
            critical_idx = (x1 + x2) / 2

        if critical_idx < 1:
            continue

        critical_idx = np.clip(critical_idx, 0, len(contour) - 1)
        cidxs.append(critical_idx)

        crit_floor = int(np.floor(critical_idx))
        crit_ceil = int(min(crit_floor + 1, len(contour) - 1))

        frac = critical_idx - crit_floor

        if d2_spl(critical_idx) < 0:
            tpt = TurningPointType.MAXIMUM
        else:
            tpt = TurningPointType.MINIMUM
        
        turningPoints.append(
            TurningPoint(
                type=tpt,
                x=(1 - frac) * contour[int(crit_floor)][0] + frac * contour[int(crit_ceil)][0],
                y=(1 - frac) * contour[int(crit_floor)][1] + frac * contour[int(crit_ceil)][1],
                B=mag_amp_spl(critical_idx)
            )
        )

    return turningPoints


def single_contour_w0_points(
        contour: npt.NDArray[np.float64],
        magnetic_amplitude_func: MagneticAmplitudeFunction
) -> List[TurningPoint]:
    """Finds the w=0 points along a given contour of constant eletric potential,
    assumes that these points are also the points of dB/dS = 0 where S parametrises
    the contour of constant U
    """
    magnetic_amplitudes = np.array([magnetic_amplitude_func(x[0], x[1]) for x in contour])
    mask = ~np.isnan(magnetic_amplitudes)
    magnetic_amplitudes = magnetic_amplitudes[mask]

    retVal = []

    valid_contours, valid_amps = __split_into_valid_segments(contour, magnetic_amplitudes)
    
    for ctr, amps in zip(valid_contours, valid_amps):
        retVal.extend(__single_contour_w0_points(ctr, amps))

    return retVal


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
        print(i)
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


def find_w0_points_in_region(
        x_bounds: tuple[float, float],
        y_bounds: tuple[float, float],
        potential_levels: List[float],
        magnetic_amplitude_func: MagneticAmplitudeFunction,
        electric_potential_func: PotentialFunction,
        final_resolution: int,
        initial_resolution: int,
) -> tuple[List[TurningPoint], List[TurningPoint]]:
    """
    Calculates and finds any of the W = 0 points in a specific region 
    for the given parameters, first testing at a lower resolution, then 
    if atelast one W = 0 point is found, calculations are redone using a higher
    resolution. This relies on the fact that lower resolutions tend to false positive

    Args:
        x_bounds (Re): the x_bounds to search within
        y_bounds (Re): the y_bounds to search within
        potential_levels (keV): a list of the potential levels to search for
        magnetic_ampltidue_func: returns the B assosciated with the 
                                 equatorial position
        electric_potential_func: returns the U assosciated with the 
                                 equatorial position
        initial_resolution: the grid sidelength to use for inital probing.
                            set to 0 to not do the initial probing
        final_resolution: the grid sidelength to use to generate the
                          equipotential contours if a positive is triggered 

    Returns:
        a tuple of lists of the W = 0 points found at the high resolution,
        the first list is the lower contour points, the second the upper contour
        points
    """

    if initial_resolution is not None:
        initial_contours = generate_equipotentials(
            electric_potential_func,
            potential_levels,
            x_bounds,
            y_bounds,
            initial_resolution
        )
        
        intial_w0_points = contour_w0_points(initial_contours, magnetic_amplitude_func)

        if not any(item for sublist1 in intial_w0_points for sublist2 in sublist1 for item in sublist2):
            return ([], [])
        
    final_contours = generate_equipotentials(
        electric_potential_func,
        potential_levels,
        x_bounds,
        y_bounds,
        final_resolution
    )

    final_w0_points = contour_w0_points(final_contours, magnetic_amplitude_func)

    for i in range(len(potential_levels)):
        for j in range(len(final_w0_points[i])):
            for ele in final_w0_points[i][j]:
                ele.U = potential_levels[i]

    return (
        parse_lower_contour_w0_points(final_w0_points),
        parse_upper_contour_w0_points(final_w0_points)
    )
