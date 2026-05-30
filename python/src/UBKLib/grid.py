from typing import (
    List,
    Tuple
)
import numpy as np
import numpy.typing as npt
import multiprocessing as mp
from concurrent.futures import (
    ProcessPoolExecutor,
    ThreadPoolExecutor
)
from skimage import measure
from scipy.interpolate import RegularGridInterpolator

from .types import (
    PotentialFunction,
    MagneticAmplitudeFunctionWithK,
    TurningPoint,
    TurningPointType,
    HamiltonianFunction,
    Vectorizable
)


def _is_closed(contour):
    return np.array_equal(contour[0], contour[-1])


def _includes_centre(contour, grid_size):
    halfway = np.array([np.array([grid_size / 2, grid_size / 2])])
    return measure.points_in_poly(halfway, contour)


def _process_row(
        row_idx: int,
        x_grid: npt.NDArray[np.float64],
        y_grid: npt.NDArray[np.float64],
        valid_mask: npt.NDArray[np.bool],
        k_values: List[float],
        magnetic_amp_func: MagneticAmplitudeFunctionWithK
) -> List[npt.NDArray[np.float64]]:
    """
    Process a single row. Returns list of arrays, one per k-value.
    """
    num_k = len(k_values)
    row_length = x_grid.shape[1]
    
    row_data = [np.full(row_length, np.nan, dtype=np.float64) for _ in range(num_k)]
    
    row_valid = valid_mask[row_idx, :]
    x_row = x_grid[row_idx, :]
    y_row = y_grid[row_idx, :]
    
    valid_cols = np.where(row_valid)[0]
    for col_idx in valid_cols:
        try:
            amplitudes = magnetic_amp_func(x_row[col_idx], y_row[col_idx], k_values)
            for k_idx, amp in enumerate(amplitudes):
                row_data[k_idx][col_idx] = amp
        except Exception:
            pass  # Keep NaN for failed points
    
    return row_data


def _get_outer_contour_by_distance(contour1, contour2):
    """
    Returns the outer contour based on average distance from centroid.
    """
    avg_dist1 = np.mean(np.linalg.norm(contour1, axis=1))
    avg_dist2 = np.mean(np.linalg.norm(contour2, axis=1))
    
    return avg_dist1 > avg_dist2


class GridWithInterp:
    __discrete_grid: npt.NDArray[np.float64]
    __interp_grid: RegularGridInterpolator

    def __init__(
            self,
            discrete_grid: npt.NDArray[np.float64],
            x_vals: npt.NDArray[np.float64],
            y_vals: npt.NDArray[np.float64]
    ):
        self.__discrete_grid = discrete_grid
        self.__interp_grid = RegularGridInterpolator(
            (x_vals, y_vals),
            self.__discrete_grid,
            bounds_error=False,
            fill_value=np.nan
        )

    @property
    def discrete(self) -> npt.NDArray[np.float64]:

        return self.__discrete_grid

    def __call__(
            self,
            x: float, 
            y: float
    ) -> float:
        return self.__interp_grid((y, x))


class FullGrid:
    
    __x_bounds: Tuple[float, float]
    __y_bounds: Tuple[float, float]
    __k_values: List[float]
    __potential_func: PotentialFunction
    __magnetic_amp_func: MagneticAmplitudeFunctionWithK

    __x_grid: npt.NDArray[np.float64]
    __y_grid: npt.NDArray[np.float64]
    __potential_grid: GridWithInterp
    __magnetic_amp_grids: List[GridWithInterp]
    __cross_product_grids: List[npt.NDArray[np.float64]]
    __hamiltonian_grids: List[npt.NDArray[np.float64]]
    __valid_mask: npt.NDArray[np.bool]
    __INVALID_RADIUS: float = 1.05

    def __init__(
            self,
            k_values: List[float],
            potential_func: PotentialFunction,
            magnetic_amp_func: MagneticAmplitudeFunctionWithK,
            x_bounds: Tuple[float, float] = (-15, 15),
            y_bounds: Tuple[float, float] = (-15, 15),
            resolution: int = 1000
    ):
        self.__potential_func = potential_func
        self.__magnetic_amp_func = magnetic_amp_func
        self.__x_bounds = x_bounds
        self.__y_bounds = y_bounds
        self.__k_values = k_values.copy()

        self.__init_real_space_grid(x_bounds, y_bounds, resolution)

    def __grid_sides(
            self,
            bounds: Tuple[float, float],
            resolution: int
    ) -> npt.NDArray[np.float64]:
        return np.linspace(bounds[0], bounds[1], resolution, endpoint=True, dtype=np.float64)

    def __x_vals(
            self,
            bounds: Tuple[float, float],
            resolution: int
    ) -> npt.NDArray[np.float64]:
        return self.__grid_sides(bounds, resolution)

    def __y_vals(
            self,
            bounds: Tuple[float, float],
            resolution: int
    ) -> npt.NDArray[np.float64]:
        return self.__grid_sides(bounds, resolution)

    def __init_real_space_grid(
            self, 
            x_bounds: Tuple[float, float],
            y_bounds: Tuple[float, float],
            resolution: int,
    ) -> None:
            
        x_vals = self.__x_vals(x_bounds, resolution)
        y_vals = self.__y_vals(y_bounds, resolution)
        self.__x_grid, self.__y_grid = np.meshgrid(x_vals, y_vals)
        r = np.hypot(self.__x_grid, self.__y_grid)
        self.__valid_mask = (r >= self.__INVALID_RADIUS)

    def calc_potential_grid(self) -> None:
        potential_grid = np.full_like(self.__x_grid, np.nan, dtype=np.float64)
        potential_grid[self.__valid_mask] = self.__potential_func(
            self.__x_grid[self.__valid_mask],
            self.__y_grid[self.__valid_mask]
        )
        self.__potential_grid = GridWithInterp(
            potential_grid,
            self.__x_grid[0, :],
            self.__y_grid[:, 1]
        )

    def calc_magnetic_amp_grids(self) -> None:
        mag_amp_grids = [
            np.full_like(self.__x_grid, np.nan, dtype=np.float64) for _ in self.__k_values
        ]
        for i in range(self.__valid_mask.shape[0]):
            for j in range(self.__valid_mask.shape[1]):
                if self.__valid_mask[i][j]:
                    val = self.__magnetic_amp_func(self.__x_grid[i][j], self.__y_grid[i][j])
                    for k in range(len(self.__k_values)):
                        mag_amp_grids[k][i, j] = val[k]

        self.__magnetic_amp_grids = []

        for grid in mag_amp_grids:
            self.__magnetic_amp_grids.append(
                GridWithInterp(
                    grid,
                    self.__x_grid[0, :],
                    self.__y_grid[:, 1]
                )
            )
    
    @property
    def x_grid(self) -> npt.NDArray[np.float64]:
        return self.__x_grid

    @property
    def y_grid(self) -> npt.NDArray[np.float64]:
        return self.__y_grid

    @property
    def potential_func(self) -> PotentialFunction:
        return self.__potential_func

    @potential_func.setter
    def potential_func(
            self,
            new_potential_func: PotentialFunction
    ) -> None:
        self.__potential_func = new_potential_func
        self.__potential_grid = None

    @property
    def potential_grid(self) -> GridWithInterp:
        return self.__potential_grid

    @property
    def magnetic_amp_grids(self) -> List[GridWithInterp]:
        return self.__magnetic_amp_grids
    
    @property 
    def hamiltonian_grids(self) -> List[npt.NDArray[np.float64]]:
        return self.__hamiltonian_grids

    def calc_cross_products_grids(self) -> None:
        self.__cross_product_grids = []
        du_y, du_x = np.gradient(self.__potential_grid.discrete)

        for grid in self.__magnetic_amp_grids:
            db_y, db_x = np.gradient(grid.discrete)
            self.__cross_product_grids.append(
                np.multiply(du_x, db_y) - np.multiply(du_y, db_x)
            )

    def calc_hamiltonian_grids(self, hamiltonian_func: HamiltonianFunction) -> None:
        self.__hamiltonian_grids = []

        for grid in self.__magnetic_amp_grids:
            self.__hamiltonian_grids.append(
                hamiltonian_func(
                    grid.discrete,
                    self.__potential_grid.discrete
                )
            )

    def __find_closed_contour(
        self,
        grid: npt.NDArray[np.float64],
        level
    ) -> None:
        contours = measure.find_contours(grid, level)

        for contour in contours:
            if _is_closed(contour) and (len(contour) > self.__x_grid.shape[0] / 10):
                if _includes_centre(contour, grid.shape[0]):
                    return contour

        return None

    def __contour_level_search_step(
        self,
        grid: npt.NDArray[np.float64],
        level_closed: float,
        level_open: float
    ) -> Tuple[float, float]:

        new_level = (level_closed + level_open) / 2

        res = self.__find_closed_contour(grid, new_level)
        
        if res is None:
            return (level_closed, new_level)
        else:
            return (new_level, level_open)

    def __find_initial_levels(
        self,
        grid: npt.NDArray[np.float64]
    ) -> Tuple[float, float]:

        valid_data = grid[~np.isnan(grid)]
        lo, hi = np.min(valid_data), np.max(valid_data)

        delta = (hi - lo) / 100

        mid = lo

        for i in range(100):
            
            mid_closed = (self.__find_closed_contour(grid, mid) is not None)

            if mid_closed is True:
                break

            mid = mid + delta

        if mid_closed is False:
            raise ValueError("Implementation error, could not find suitable starting value for H contours")
        
        return (mid, lo, hi)

    def __interp_contour(
            self,
            contour
    ) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:

        x_vals = np.interp(
            contour[:, 1],
            [0, self.__x_grid.shape[0] - 1],
            [self.__x_bounds[0], self.__x_bounds[1]]
        )
        y_vals = np.interp(
            contour[:, 0],
            [0, self.__y_grid.shape[0] - 1],
            [self.__y_bounds[0], self.__y_bounds[1]]
        )
        return x_vals, y_vals

    def find_cross_product_zeros(self) -> List[List[TurningPoint]]:

        retVal = []

        for j, grid in enumerate(self.__cross_product_grids):
            temp = []
            contours = measure.find_contours(grid, level=0)

            for contour in contours:
                if len(contour) < self.__x_grid.shape[0] / 10:
                    continue
                
                x_vals, y_vals = self.__interp_contour(contour)
                temp.append([TurningPoint(
                    type=TurningPointType.MAXIMUM,
                    x=x_vals[i],
                    y=y_vals[i],
                    B=self.__magnetic_amp_grids[j](x_vals[i], y_vals[i]),
                    U=self.__potential_grid(x_vals[i], y_vals[i])
                ) for i in range(contour.shape[0])])

            retVal.append(temp)

        return retVal

    def __find_lcds_contours(
            self,
            grid: npt.NDArray[np.float64],
            iterations: int,
            closed_level: float,
            open_level: float  
    ) -> Tuple[float, float, npt.NDArray[np.float64]]:
        closed = closed_level
        open = open_level
        for i in range(iterations):
            closed, open = self.__contour_level_search_step(grid, closed, open)
    
        return closed, open, self.__find_closed_contour(grid, closed)

    def find_lcds_contours(
            self,
            iterations: int = 50,
            check_iterations: int = 10
    ) -> List[Tuple(float, npt.NDArray, npt.NDArray)]:

        assert check_iterations < iterations

        retVal = []

        for grid in self.__hamiltonian_grids:
            mid, lo, hi = self.__find_initial_levels(grid)
            midlo = mid
            midhi = mid
            lo_closed, lo_open, lo_contour = self.__find_lcds_contours(grid, check_iterations, midlo, lo)
            hi_closed, hi_open, hi_contour = self.__find_lcds_contours(grid, check_iterations, midhi, hi)
            if _get_outer_contour_by_distance(lo_contour, hi_contour):
                closed, open, contour = self.__find_lcds_contours(grid, iterations - check_iterations, lo_closed, lo_open)
            else:
                closed, open, contour = self.__find_lcds_contours(grid, iterations - check_iterations, hi_closed, hi_open)
            
            x_vals, y_vals = self.__interp_contour(contour)
            retVal.append((closed, x_vals, y_vals))

        return retVal

    def find_drift_shell(
            self,
            hamiltonian: float,
            k_idx: int
    ) -> List[tuple(npt.NDArray, npt.NDArray)]:
        contours = measure.find_contours(self.__hamiltonian_grids[k_idx], level=hamiltonian)

        x_vals = []
        y_vals = []
        for contour in contours:
            x, y = self.__interp_contour(contour)
            x_vals.append(x)
            y_vals.append(y)

        return x_vals, y_vals

    def calc_magnetic_amp_grid_parallel(
            self, 
            num_workers: int = None, 
            use_processes: bool = True
    ) -> None:
        """
        Fill magnetic amplitude grids using parallel processing across rows.
        
        Parameters:
        -----------
        num_workers : int, optional
            Number of worker processes/threads. Defaults to CPU count.
        use_processes : bool, optional
            If True, use ProcessPoolExecutor. If False, use ThreadPoolExecutor.
        """
        if num_workers is None:
            num_workers = mp.cpu_count()
        
        rows = list(range(self.__x_grid.shape[0]))
        
        Executor = ProcessPoolExecutor if use_processes else ThreadPoolExecutor
        
        with Executor(max_workers=num_workers) as executor:
            futures = {
                executor.submit(
                    _process_row,
                    row_idx,
                    self.__x_grid,
                    self.__y_grid,
                    self.__valid_mask,
                    self.__k_values,
                    self.__magnetic_amp_func
                ): row_idx
                for row_idx in rows
            }
            
            mag_amp_grids = [
                np.full_like(self.__x_grid, np.nan, dtype=np.float64) for _ in self.__k_values
            ]
            for future in futures:
                row_idx = futures[future]
                try:
                    row_data = future.result()
                    for k_idx in range(len(self.__k_values)):
                        mag_amp_grids[k_idx][row_idx, :] = row_data[k_idx]
                except Exception as e:
                    print(f"Error processing row {row_idx}: {e}")

        self.__magnetic_amp_grids = []
        for grid in mag_amp_grids:
            self.__magnetic_amp_grids.append(
                GridWithInterp(
                    grid,
                    self.__x_grid[0, :],
                    self.__y_grid[:, 1]
                )
            )

    def save(self, filepath: str) -> None:
        """
        Save grid data to a .npz file.
        
        Parameters:
        -----------
        filepath : str
            Path to save the file (without extension, .npz will be added).
        """
        if not filepath.endswith('.npz'):
            filepath += '.npz'
        
        save_dict = {
            'x_grid': self.__x_grid,
            'y_grid': self.__y_grid,
            'potential_grid': self.__potential_grid.discrete,
            'valid_mask': self.__valid_mask,
            'k_values': np.array(self.__k_values),
            'x_bounds': np.array(self.__x_bounds),
            'y_bounds': np.array(self.__y_bounds),
        }
        
        for k_idx, grid in enumerate(self.__magnetic_amp_grids):
            save_dict[f'magnetic_amp_grid_{k_idx}'] = grid.discrete
        
        np.savez_compressed(filepath, **save_dict)

    @classmethod
    def load(cls, filepath: str) -> 'FullGrid':
        """
        Load grid data from a .npz file. Note: This restores the data grids
        but the potential_func and magnetic_amp_func must be provided separately.
        
        Parameters:
        -----------
        filepath : str
            Path to the .npz file.
            
        Returns:
        --------
        Grid
            A Grid instance with restored data but dummy functions.
        """
        if not filepath.endswith('.npz'):
            filepath += '.npz'
        
        data = np.load(filepath)
        
        grid = cls(
            k_values=data['k_values'].tolist(),
            potential_func=None,
            magnetic_amp_func=None,
            x_bounds=tuple(data['x_bounds']),
            y_bounds=tuple(data['y_bounds']),
            resolution=data['x_grid'].shape[0]
        )
        
        grid.__potential_grid = GridWithInterp(
            data['potential_grid'],
            grid.__x_grid[0, :],
            grid.__y_grid[:, 1]
        )
        
        k_idx = 0
        grid.__magnetic_amp_grids = []
        while f'magnetic_amp_grid_{k_idx}' in data:
            grid.__magnetic_amp_grids.append(GridWithInterp(
                data[f'magnetic_amp_grid_{k_idx}'],
                grid.__x_grid[0, :],
                grid.__y_grid[:, 1]
            ))
            k_idx += 1
        
        return grid
