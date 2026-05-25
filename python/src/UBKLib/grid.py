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
    TurningPointType
)


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


class Grid:
    
    __x_bounds: Tuple[float, float]
    __y_bounds: Tuple[float, float]
    __k_values: List[float]
    __potential_func: PotentialFunction
    __magnetic_amp_func: MagneticAmplitudeFunctionWithK

    __x_grid: npt.NDArray[np.float64]
    __y_grid: npt.NDArray[np.float64]
    __potential_grid: npt.NDArray[np.float64]
    __magnetic_amp_grids: List[npt.NDArray[np.float64]]
    __cross_product_grids: List[npt.NDArray[np.float64]]
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
        self.__init_field_grids(len(k_values))

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

    def __init_field_grids(
            self,
            num_k_values: int
    ) -> None:
        self.__potential_grid = np.full_like(self.__x_grid, np.nan, dtype=np.float64)
        self.__magnetic_amp_grids = [
            np.full_like(self.__x_grid, np.nan, dtype=np.float64) for _ in range(num_k_values)
        ]

    def calc_potential_grid(self) -> None:
        self.__potential_grid[self.__valid_mask] = self.__potential_func(
            self.__x_grid[self.__valid_mask],
            self.__y_grid[self.__valid_mask]
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
        self.__potential_grid.fill(np.nan)

    @property
    def potential_grid(self) -> npt.NDArray[np.float64]:
        return self.__potential_grid

    @property
    def magnetic_amp_grids(self) -> List[npt.NDArray[np.float64]]:
        return self.__magnetic_amp_grids

    def calc_cross_products(self) -> None:
        self.__cross_product_grids = []
        du_y, du_x = np.gradient(self.__potential_grid)

        for grid in self.__magnetic_amp_grids:
            db_y, db_x = np.gradient(grid)
            self.__cross_product_grids.append(
                np.multiply(du_x, db_y) - np.multiply(du_y, db_x)
            )

    def find_cross_product_zeros(self) -> List[List[TurningPoint]]:
        interp_pot_grid = RegularGridInterpolator(
            (self.__x_grid[0, :], self.__y_grid[:, 1]),
            self.__potential_grid,
            bounds_error=False, 
            fill_value=np.nan
        )

        retVal = []

        for j, grid in enumerate(self.__cross_product_grids):
            temp = []
            contours = measure.find_contours(grid, level=0)

            for contour in contours:
                if len(contour) < self.__x_grid.shape[0] / 10:
                    continue

                interp = RegularGridInterpolator(
                    (self.__x_grid[0, :], self.__y_grid[:, 1]),
                    self.magnetic_amp_grids[j],
                    bounds_error=False, 
                    fill_value=np.nan
                )

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
                temp.append([TurningPoint(
                    type=TurningPointType.MAXIMUM,
                    x=x_vals[i],
                    y=y_vals[i],
                    B=interp((y_vals[i], x_vals[i])),
                    U=interp_pot_grid((y_vals[i], x_vals[i]))
                ) for i in range(contour.shape[0])])

            retVal.append(temp)

        return retVal

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
            
            for future in futures:
                row_idx = futures[future]
                try:
                    row_data = future.result()
                    for k_idx in range(len(self.__k_values)):
                        self.__magnetic_amp_grids[k_idx][row_idx, :] = row_data[k_idx]
                except Exception as e:
                    print(f"Error processing row {row_idx}: {e}")

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
            'potential_grid': self.__potential_grid,
            'valid_mask': self.__valid_mask,
            'k_values': np.array(self.__k_values),
            'x_bounds': np.array(self.__x_bounds),
            'y_bounds': np.array(self.__y_bounds),
        }
        
        for k_idx, grid in enumerate(self.__magnetic_amp_grids):
            save_dict[f'magnetic_amp_grid_{k_idx}'] = grid
        
        np.savez_compressed(filepath, **save_dict)

    @classmethod
    def load(cls, filepath: str) -> 'Grid':
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
        
        grid.__potential_grid = data['potential_grid']
        
        k_idx = 0
        while f'magnetic_amp_grid_{k_idx}' in data:
            grid.__magnetic_amp_grids[k_idx] = data[f'magnetic_amp_grid_{k_idx}']
            k_idx += 1
        
        return grid
