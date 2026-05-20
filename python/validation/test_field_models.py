import pytest
import numpy as np

from UBKLib import (
    volland_stern_potential, cross_tail_potential, equatorial_dipole_amplitude
)

from UBKLib.field_models import cross_tail, volland_stern


class TestCrossTailPotential:
    
    @pytest.mark.parametrize("x, y, expected", [
        # r=1.0. Corotation: -92.4/1. Convection: 0. 
        (1.0, 0.0, -92.4),         
        
        # r=1.0. Corotation: -92.4/1. Convection: -10(-1). 
        (0.0, 1.0, -102.4),        
        
        # r=1.0. Corotation: -92.4/1. Convection: 10(1). 
        (0.0, -1.0, -82.4),        
        
        # r=5.0. Corotation: -92.4/5 = -18.44. Convection: 10(-4) = -40. 
        (3.0, 4.0, -58.48),        
    ])
    def test_scalar_coordinates(self, x, y, expected):
        result = cross_tail_potential(x, y)
        assert result == pytest.approx(expected)

    def test_numpy_arrays(self):
        x_arr = np.array([1.0, 0.0, 3.0])
        y_arr = np.array([0.0, 1.0, 4.0])
        expected = np.array([-92.4, -102.4, -58.48])
        
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
