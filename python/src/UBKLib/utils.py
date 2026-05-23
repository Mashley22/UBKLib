import numpy as np
import numpy.typing as npt


def random_xy(
        r_min: float = 1.0,
        r_max: float = 15.0,
        n: float = 1
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """
    Generate random (x, y) coordinates with radial distance between r_min and r_max.
    
    Parameters:
        r_min: minimum radial distance (default 1.0)
        r_max: maximum radial distance (default 15.0)
        n: number of points to generate
    
    Returns:
        x, y: arrays of coordinates
    """
    theta = np.random.uniform(0, 2 * np.pi, n)
    
    r = np.sqrt(np.random.uniform(r_min**2, r_max**2, n))
    
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    
    return x, y


def generate_log_space_bins(
        min_val: float,
        max_val: float,
        num: int
) -> npt.NDArray[np.float64]:
    """
    Generates evenly spaced bins in the log space between min and max.
    It returns the bin edges in the original scale. 
    
    Args:
        min_val: minimum value (must be > 0)
        max_val: maximum value
        num: number of bins (edges will have length num + 1)
    
    Returns:
        array of bin edges in original scale, length num + 1
    """
    if min_val <= 0:
        min_val = np.finfo(float).eps

    return np.logspace(np.log10(min_val), np.log10(max_val), num + 1)


def generate_inv_cuberoot_bins(
        min_val: float,
        max_val: float,
        num: int
) -> npt.NDArray[np.float64]:
    """
    Generates bins evenly spaced in 1/cuberoot(B) space,
    which corresponds to linear spacing in r for a dipole field (B ~ 1/r³).
    
    Args:
        min_val: minimum B value (can be 0)
        max_val: maximum B value
        num: number of bins (edges will have length num + 1)
    
    Returns:
        array of bin edges in original B scale, length num + 1
    """
    if min_val <= 0:
        min_val = np.finfo(float).eps
    
    inv_cbrt_min = 1.0 / np.cbrt(max_val) 
    inv_cbrt_max = 1.0 / np.cbrt(min_val) 
    
    inv_cbrt_edges = np.linspace(inv_cbrt_min, inv_cbrt_max, num + 1)
    
    B_edges = 1.0 / inv_cbrt_edges**3
    
    return np.sort(B_edges)
