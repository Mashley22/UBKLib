import pytest
import numpy as np

from UBKLib import (
    cross_tail_potential,
    volland_stern_potential
)

class TestCrossTailPotential:
    
    @pytest.mark.parametrize("x, y, expected", [
        # r=1.0. Corotation: -94.2/1. Convection: 0. Total: -94.2
        (1.0, 0.0, -94.2),         
        
        # r=1.0. Corotation: -94.2/1. Convection: 10(1). Total: -84.2
        (0.0, 1.0, -84.2),        
        
        # r=1.0. Corotation: -94.2/1. Convection: -10(-1). Total: -104.2
        (0.0, -1.0, -104.2),        
        
        # r=5.0. Corotation: -94.2/5 = -18.84. Convection: 10(4) = 40. Total: 21.16
        (3.0, 4.0, 21.16),        
    ])
    def test_scalar_coordinates(self, x, y, expected):
        result = cross_tail_potential(x, y)
        assert result == pytest.approx(expected)

    def test_numpy_arrays(self):
        x_arr = np.array([1.0, 0.0, 3.0])
        y_arr = np.array([0.0, 1.0, 4.0])
        expected = np.array([-94.2, -84.2, 21.16])
        
        result = cross_tail_potential(x_arr, y_arr)
        
        assert isinstance(result, np.ndarray)
        np.testing.assert_allclose(result, expected)

    @pytest.mark.parametrize("x, y", [
        (0.0, 0.0),
        (0.5, 0.0),
        (0.0, 0.5)
    ])
    def test_radial_dist_below_surface(self, x, y):
        with pytest.raises(ValueError):
            cross_tail_potential(0.0, 0.0)

class TestVollandSternPotential:

    @pytest.mark.parametrize("kp, expected_total", [
        # Kp = 0.0: Denominator is 1.0. 
        # e0 = 0.045. Convection = -0.045 * 5.0 * 4.0 = -0.9. 
        # Total = -18.84 + (-0.9) = -19.74
        (0.0, -19.74),
        
        # Kp = 3.0: Denom = 1.0 - 0.0159(3) + 0.0093(9) = 1.036
        # e0 = 0.045 / (1.036^3) = 0.040473
        # Convection = -0.040473 * 20.0 = -0.80946
        # Total = -18.84 + (-0.80946) = -19.64946
        (3.0, -19.6494676),
        
        # Kp = 5.0: Denom = 1.0 - 0.0159(5) + 0.0093(25) = 1.153
        # e0 = 0.045 / (1.153^3) = 0.029358
        # Convection = -0.029358 * 20.0 = -0.58715
        # Total = -18.84 + (-0.58715) = -19.42715
        (5.0, -19.4271516),
    ])
    def test_kp_scalar_values(self, kp, expected_total):
        result = volland_stern_potential(x=3.0, y=4.0, kp=kp)
        
        assert result == pytest.approx(expected_total, rel=1e-5)

    def test_numpy_arrays(self):
        x_arr = np.array([3.0, 0.0, -4.0])
        y_arr = np.array([4.0, 5.0,  3.0])
        
        result = volland_stern_potential(x_arr, y_arr, kp=3.0)
        
        assert isinstance(result, np.ndarray)
        assert result.shape == (3,)
        
        for i in range(x_arr.shape[0]):
            scalar_result = volland_stern_potential(x=x_arr[i], y=y_arr[i], kp=3.0)
            assert result[i] == pytest.approx(scalar_result)

    def test_dawn_dusk_asymmetry(self):
        pot_dusk = volland_stern_potential(x=0.0, y=5.0)
        pot_dawn = volland_stern_potential(x=0.0, y=-5.0)
        
        assert pot_dusk < pot_dawn
        
        corotation = -18.84
        convection_dusk = pot_dusk - corotation
        convection_dawn = pot_dawn - corotation
        
        assert convection_dusk == pytest.approx(-convection_dawn)

    @pytest.mark.parametrize("x, y", [
        (0.0, 0.0),    # Origin (r=0)
        (0.5, 0.5),    # Inside Earth (r ≈ 0.707)
        (0.8, 0.0),    # Inside Earth on X-axis (r=0.8)
    ])
    def test_invalid_radial_distance_scalars(self, x, y):
        with pytest.raises(ValueError, match="Coordinates must supply a radial distance larger than 1 Re"):
            volland_stern_potential(x, y)
