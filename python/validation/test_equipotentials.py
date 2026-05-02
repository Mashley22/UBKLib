import pytest
import numpy as np
from UBKLib import generate_equipotentials


def radial_potential(x, y):
    """Simple radial potential for testing."""
    return np.hypot(x, y)

def linear_potential(x, y):
    """Linear potential for testing."""
    return x + y

def constant_potential(x, y):
    """Constant potential for testing."""
    return np.ones_like(x) * 2.0

class TestGenerateEquipotentials:
    
    @pytest.fixture(autouse=True)
    def setup(self):
        self.levels = [1.0, 2.0, 3.0]
        self.x_bounds = (-5.0, 5.0)
        self.y_bounds = (-5.0, 5.0)
        self.resolution = 100
    
    def test_returns_list_of_correct_length(self):
        result = generate_equipotentials(
            radial_potential,
            self.levels,
            self.x_bounds,
            self.y_bounds,
            self.resolution
        )
        
        assert isinstance(result, list)
        assert len(result) == len(self.levels)
    
    def test_each_level_is_list_of_arrays(self):
        result = generate_equipotentials(
            radial_potential,
            self.levels,
            self.x_bounds,
            self.y_bounds,
            self.resolution
        )
        
        for level_contours in result:
            assert isinstance(level_contours, list)
            for contour in level_contours:
                assert isinstance(contour, np.ndarray)
    
    def test_contour_arrays_have_shape_n_2(self):
        result = generate_equipotentials(
            radial_potential,
            self.levels,
            self.x_bounds,
            self.y_bounds,
            self.resolution
        )
        
        for level_contours in result:
            for contour in level_contours:
                assert contour.ndim == 2
                assert contour.shape[1] == 2
    
    def test_contours_found_for_valid_levels(self):
        result = generate_equipotentials(
            radial_potential,
            [2.0],
            self.x_bounds,
            self.y_bounds,
            self.resolution
        )
        
        assert len(result[0]) > 0
    
    def test_no_contours_for_out_of_range_level(self):
        result = generate_equipotentials(
            radial_potential,
            [100.0],
            self.x_bounds,
            self.y_bounds,
            self.resolution
        )
        
        assert len(result[0]) == 0
    
    def test_default_parameters_work(self):
        result = generate_equipotentials(
            radial_potential,
            [1.0]
        )
        
        assert isinstance(result, list)
        assert len(result) == 1
    
    def test_x_bounds_min_greater_than_max_raises_error(self):
        with pytest.raises(AssertionError):
            generate_equipotentials(
                radial_potential,
                self.levels,
                x_bounds=(10.0, -10.0)
            )
    
    def test_y_bounds_min_equal_to_max_raises_error(self):
        with pytest.raises(AssertionError):
            generate_equipotentials(
                radial_potential,
                self.levels,
                y_bounds=(15.0, 15.0)
            )
    
    def test_equal_x_bounds_raises_error(self):
        with pytest.raises(AssertionError):
            generate_equipotentials(
                radial_potential,
                self.levels,
                x_bounds=(5.0, 5.0)
            )
    
    def test_empty_levels_returns_empty_list(self):
        result = generate_equipotentials(
            radial_potential,
            [],
            self.x_bounds,
            self.y_bounds,
            self.resolution
        )
        
        assert isinstance(result, list)
        assert result == []
    
    def test_single_level_single_contour_structure(self):
        result = generate_equipotentials(
            radial_potential,
            [1.5],
            self.x_bounds,
            self.y_bounds,
            self.resolution
        )
        
        assert len(result) == 1
        contour = result[0][0]
        assert isinstance(contour, np.ndarray)
        assert contour.shape[1] == 2
    
    def test_higher_resolution_produces_smoother_contours(self):
        low_res = generate_equipotentials(
            radial_potential,
            [2.0],
            self.x_bounds,
            self.y_bounds,
            resolution=10
        )
        
        high_res = generate_equipotentials(
            radial_potential,
            [2.0],
            self.x_bounds,
            self.y_bounds,
            resolution=100
        )
        
        assert len(high_res[0][0]) > len(low_res[0][0])
    
    def test_different_levels_produce_different_contours(self):
        result = generate_equipotentials(
            radial_potential,
            [1.5, 3.0],
            self.x_bounds,
            self.y_bounds,
            self.resolution
        )
        
        contour_1 = result[0][0]
        contour_3 = result[1][0]
        
        mean_radius_1 = np.mean(np.sqrt(contour_1[:, 0]**2 + contour_1[:, 1]**2))
        mean_radius_3 = np.mean(np.sqrt(contour_3[:, 0]**2 + contour_3[:, 1]**2))
        
        assert mean_radius_3 > mean_radius_1
    
    def test_asymmetric_bounds(self):
        result = generate_equipotentials(
            radial_potential,
            [2.0],
            x_bounds=(-5.0, 10.0),
            y_bounds=(-10.0, 5.0),
            resolution=30
        )
        
        assert isinstance(result, list)
        assert len(result) == 1


class TestGenerateEquipotentialsCustomPotentials:
    
    def test_constant_potential(self):
        result = generate_equipotentials(
            constant_potential,
            [2.0],
            (-5.0, 5.0),
            (-5.0, 5.0),
            50
        )
        
        assert isinstance(result, list)
    
    def test_linear_potential(self):
        result = generate_equipotentials(
            linear_potential,
            [0.0],
            (-5.0, 5.0),
            (-5.0, 5.0),
            50
        )
        
        assert isinstance(result, list)
        assert len(result) == 1


@pytest.mark.parametrize("levels,expected_length", [
    ([1.0], 1),
    ([1.0, 2.0], 2),
    ([1.0, 2.0, 3.0, 4.0], 4),
    ([], 0),
])
def test_output_length_matches_levels_length(levels, expected_length):
    result = generate_equipotentials(
        radial_potential,
        levels,
        (-5.0, 5.0),
        (-5.0, 5.0),
        30
    )
    
    assert len(result) == expected_length


@pytest.mark.parametrize("invalid_bounds", [
    pytest.param({"x_bounds": (10.0, -10.0)}, id="x_min_greater_than_max"),
    pytest.param({"x_bounds": (5.0, 5.0)}, id="x_min_equal_to_max"),
    pytest.param({"y_bounds": (15.0, 15.0)}, id="y_min_equal_to_max"),
    pytest.param({"y_bounds": (10.0, -10.0)}, id="y_min_greater_than_max"),
])
def test_invalid_bounds_raise_error(invalid_bounds):
    with pytest.raises(AssertionError):
        generate_equipotentials(
            radial_potential,
            [1.0],
            **invalid_bounds
        )


def test_all_contours_are_finite():
    result = generate_equipotentials(
        radial_potential,
        [1.0, 2.0, 3.0],
        (-5.0, 5.0),
        (-5.0, 5.0),
        50
    )
    
    for level_contours in result:
        for contour in level_contours:
            assert np.all(np.isfinite(contour))
