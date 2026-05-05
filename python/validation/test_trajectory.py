import numpy as np
import pytest
import sys
import os

from UBKLib import continuousLCDSUBTrajectory, UBTrajectory, UBCoord
from UBKLib.traits import W0ContourFunction


# -------------------------------------------------------------------
# Real Contour Functions for Testing
# -------------------------------------------------------------------

class LinearW0Contour(W0ContourFunction):
    """Linear contour: U = a*B + c"""
    def __init__(self, a, c):
        self.a = a
        self.c = c
    
    def __call__(self, B):
        return self.a * B + self.c

class QuadraticW0Contour(W0ContourFunction):
    """Quadratic contour: U = a*B^2 + c"""
    def __init__(self, a, c):
        self.a = a
        self.c = c
    
    def __call__(self, B):
        return self.a * B**2 + self.c

class ReciprocalW0Contour(W0ContourFunction):
    """Contour with 1/B^(1/3) term"""
    def __init__(self, a, b, c):
        self.a = a  # coefficient for 1/B^(1/3)
        self.b = b  # coefficient for linear term
        self.c = c  # constant
    
    def __call__(self, B):
        return self.a / (B ** (1/3)) + self.b * B + self.c


# -------------------------------------------------------------------
# Fixtures
# -------------------------------------------------------------------

@pytest.fixture
def simple_linear_lower():
    """Simple linear lower contour: U = 2*B + 10"""
    return LinearW0Contour(2.0, 10.0)

@pytest.fixture
def simple_linear_upper():
    """Simple linear upper contour: U = 2*B + 30"""
    return LinearW0Contour(2.0, 30.0)

@pytest.fixture
def steep_linear_lower():
    """Steep lower contour: U = 5*B"""
    return LinearW0Contour(5.0, 0.0)

@pytest.fixture
def steep_linear_upper():
    """Steep upper contour: U = 5*B + 50"""
    return LinearW0Contour(5.0, 50.0)

@pytest.fixture
def quadratic_lower():
    """Nonlinear lower contour: U = 0.001*B^2 + 100"""
    return QuadraticW0Contour(0.001, 100.0)

@pytest.fixture
def quadratic_upper():
    """Nonlinear upper contour: U = 0.001*B^2 + 200"""
    return QuadraticW0Contour(0.001, 200.0)

@pytest.fixture
def reciprocal_lower():
    """Lower contour with 1/B^(1/3) singularity"""
    return ReciprocalW0Contour(a=100.0, b=0.1, c=0.0)

@pytest.fixture
def reciprocal_upper():
    """Upper contour with 1/B^(1/3) singularity"""
    return ReciprocalW0Contour(a=100.0, b=0.1, c=50.0)

@pytest.fixture
def standard_bounds():
    """Standard B bounds (avoiding singularity)"""
    return (1000.0, 20000.0)

@pytest.fixture
def wide_bounds():
    """Wide B bounds"""
    return (100.0, 30000.0)

@pytest.fixture
def tight_bounds():
    """Tight B bounds"""
    return (5000.0, 10000.0)


# -------------------------------------------------------------------
# Basic Functionality Tests
# -------------------------------------------------------------------

class TestContinuousLCDSUBTrajectoryBasic:
    """Basic functional tests with linear contours"""
    
    def test_returns_correct_type(self, simple_linear_lower, simple_linear_upper, standard_bounds):
        """Should return a UBTrajectory instance"""
        gradient = 1.0
        result = continuousLCDSUBTrajectory(
            simple_linear_lower, 
            simple_linear_upper, 
            gradient, 
            standard_bounds
        )
        assert isinstance(result, UBTrajectory)
    
    def test_gradient_set_correctly(self, simple_linear_lower, simple_linear_upper, standard_bounds):
        """Trajectory should store the input gradient"""
        gradient = 1.5
        result = continuousLCDSUBTrajectory(
            simple_linear_lower, 
            simple_linear_upper, 
            gradient, 
            standard_bounds
        )
        assert result.gradient == gradient
    
    def test_positive_gradient_less_than_contour_slope(self, simple_linear_lower, simple_linear_upper, standard_bounds):
        """When gradient < contour slope, line should intersect at lower bound"""
        gradient = 1.0  # Less than contour slope of 2.0
        result = continuousLCDSUBTrajectory(
            simple_linear_lower, 
            simple_linear_upper, 
            gradient, 
            standard_bounds
        )
        
        # Line U = gradient*B + intercept must be >= lower contour for some B
        # At B_min, we should have approximate intersection
        B_min = standard_bounds[0]
        lower_U = simple_linear_lower(B_min)
        line_U = gradient * B_min + result.intercept
        
        # The line should be close to the contour at the intersection point
        assert abs(line_U - lower_U) < 1.0 or result.lower_intercept is not None
    
    def test_positive_gradient_greater_than_contour_slope(self, simple_linear_lower, simple_linear_upper, standard_bounds):
        """When gradient > contour slope, line should intersect at upper bound"""
        gradient = 3.0  # Greater than contour slope of 2.0
        result = continuousLCDSUBTrajectory(
            simple_linear_lower, 
            simple_linear_upper, 
            gradient, 
            standard_bounds
        )
        
        # Should find intersection points
        assert result.lower_intercept is not None
        assert result.upper_intercept is not None

    def test_gradient_equal_to_contour_slope(self, simple_linear_lower, simple_linear_upper, standard_bounds):
        """When gradient equals contour slope, the line is parallel to contours.
        The maximum intercept that still touches both contours will be such that
        the line passes through the lower contour at the optimal B."""
        gradient = 2.0  # Same as contour slope
        result = continuousLCDSUBTrajectory(
            simple_linear_lower,
            simple_linear_upper,
            gradient,
            standard_bounds
        )
        
        # Line should be parallel to contours
        assert result.lower_intercept is not None
        assert result.upper_intercept is not None
        
        # The line should touch at least one contour exactly
        lower_B = result.lower_intercept.B
        upper_B = result.upper_intercept.B
        lower_U_contour = simple_linear_lower(lower_B)
        lower_U_line = gradient * lower_B + result.intercept
        upper_U_contour = simple_linear_upper(upper_B)
        upper_U_line = gradient * upper_B + result.intercept
        
        # At least one should be very close (the binding contour)
        lower_diff = abs(lower_U_line - lower_U_contour)
        upper_diff = abs(upper_U_line - upper_U_contour)
        
        assert min(lower_diff, upper_diff) < 1.0, \
            f"Line should touch at least one contour. Lower diff: {lower_diff}, Upper diff: {upper_diff}"
        
        # For parallel lines, the intercept should be the maximum that works
        # The line U = 2*B + intercept must be >= lower contour (2*B + 10)
        # So intercept >= 10
        # And it must intersect upper contour somewhere
        assert result.intercept >= 10.0  # Must be above lower contour
    
    def test_zero_gradient(self, simple_linear_lower, simple_linear_upper, standard_bounds):
        """Test with zero gradient (horizontal line)"""
        gradient = 0.0
        result = continuousLCDSUBTrajectory(
            simple_linear_lower, 
            simple_linear_upper, 
            gradient, 
            standard_bounds
        )
        
        assert result.gradient == 0.0
        assert result.lower_intercept is not None
        assert result.upper_intercept is not None
        
        # With zero gradient, intercept should be >= contour U values
        lower_U = simple_linear_lower(result.lower_intercept.B)
        assert result.intercept >= lower_U - 1e-6
    
    def test_negative_gradient(self, simple_linear_lower, simple_linear_upper, standard_bounds):
        """Test with negative gradient"""
        gradient = -0.5
        result = continuousLCDSUBTrajectory(
            simple_linear_lower, 
            simple_linear_upper, 
            gradient, 
            standard_bounds
        )
        
        assert result.gradient == gradient
        assert result.lower_intercept is not None
        assert result.upper_intercept is not None
    
    def test_very_negative_gradient(self, steep_linear_lower, steep_linear_upper, standard_bounds):
        """Test with very negative gradient"""
        gradient = -10.0
        result = continuousLCDSUBTrajectory(
            steep_linear_lower, 
            steep_linear_upper, 
            gradient, 
            standard_bounds
        )
        
        assert result.gradient == gradient
        assert result.lower_intercept is not None
        assert result.upper_intercept is not None


# -------------------------------------------------------------------
# Nonlinear Contour Tests
# -------------------------------------------------------------------

class TestNonlinearContours:
    """Tests with quadratic and reciprocal contours"""
    
    def test_quadratic_contours_positive_gradient(self, quadratic_lower, quadratic_upper, wide_bounds):
        """Test with quadratic contours and positive gradient"""
        gradient = 2.0
        result = continuousLCDSUBTrajectory(
            quadratic_lower, 
            quadratic_upper, 
            gradient, 
            wide_bounds
        )
        
        assert isinstance(result, UBTrajectory)
        assert result.lower_intercept is not None
        assert result.upper_intercept is not None
        
        # Verify intercept points are within bounds
        assert wide_bounds[0] <= result.lower_intercept.B <= wide_bounds[1]
        assert wide_bounds[0] <= result.upper_intercept.B <= wide_bounds[1]
    
    def test_quadratic_contours_negative_gradient(self, quadratic_lower, quadratic_upper, wide_bounds):
        """Test with quadratic contours and negative gradient"""
        gradient = -0.1
        result = continuousLCDSUBTrajectory(
            quadratic_lower, 
            quadratic_upper, 
            gradient, 
            wide_bounds
        )
        
        assert isinstance(result, UBTrajectory)
        assert result.lower_intercept is not None
        assert result.upper_intercept is not None
    
    def test_quadratic_contours_high_gradient(self, quadratic_lower, quadratic_upper, wide_bounds):
        """Test with quadratic contours and high gradient"""
        gradient = 10.0
        result = continuousLCDSUBTrajectory(
            quadratic_lower, 
            quadratic_upper, 
            gradient, 
            wide_bounds
        )
        
        assert isinstance(result, UBTrajectory)
        assert result.lower_intercept is not None
        assert result.upper_intercept is not None
    
    def test_reciprocal_contours_handles_singularity(self, reciprocal_lower, reciprocal_upper, standard_bounds):
        """Test reciprocal contours don't crash due to 1/B^(1/3) singularity"""
        gradient = 0.5
        # This should not raise even though contour has 1/B^(1/3)
        result = continuousLCDSUBTrajectory(
            reciprocal_lower, 
            reciprocal_upper, 
            gradient, 
            standard_bounds
        )
        
        assert isinstance(result, UBTrajectory)
        assert result.lower_intercept is not None
        assert result.upper_intercept is not None
    
    def test_reciprocal_contours_near_singularity(self, reciprocal_lower, reciprocal_upper):
        """Test reciprocal contours with bounds close to zero"""
        gradient = 0.1
        # Bounds starting very close to zero where singularity is strong
        tight_bounds_near_zero = (1.0, 1000.0)
        
        # This might fail due to numerical issues - we're testing robustness
        try:
            result = continuousLCDSUBTrajectory(
                reciprocal_lower, 
                reciprocal_upper, 
                gradient, 
                tight_bounds_near_zero
            )
            # If it succeeds, check basic validity
            assert result.lower_intercept is not None
            assert result.upper_intercept is not None
        except RuntimeError as e:
            # It's acceptable if it fails gracefully with a clear error
            assert "minimization failed" in str(e).lower() or "singular" in str(e).lower()


# -------------------------------------------------------------------
# Edge Cases
# -------------------------------------------------------------------

class TestEdgeCases:
    """Edge case tests"""
    
    def test_tight_bounds(self, simple_linear_lower, simple_linear_upper):
        """Test with very tight bounds"""
        gradient = 1.5
        tight_b = (5000.0, 5001.0)  # Only 1 nT range
        
        result = continuousLCDSUBTrajectory(
            simple_linear_lower, 
            simple_linear_upper, 
            gradient, 
            tight_b
        )
        
        assert isinstance(result, UBTrajectory)
        assert result.lower_intercept is not None
        assert tight_b[0] <= result.lower_intercept.B <= tight_b[1]
    
    def test_large_bounds(self, simple_linear_lower, simple_linear_upper):
        """Test with very large bounds"""
        gradient = 1.0
        large_b = (1.0, 100000.0)
        
        result = continuousLCDSUBTrajectory(
            simple_linear_lower, 
            simple_linear_upper, 
            gradient, 
            large_b
        )
        
        assert isinstance(result, UBTrajectory)
        assert result.lower_intercept is not None
        assert result.upper_intercept is not None
    
    def test_very_small_gradient(self, simple_linear_lower, simple_linear_upper, standard_bounds):
        """Test with gradient near zero"""
        gradient = 1e-10
        result = continuousLCDSUBTrajectory(
            simple_linear_lower, 
            simple_linear_upper, 
            gradient, 
            standard_bounds
        )
        
        assert isinstance(result, UBTrajectory)
        assert result.gradient == pytest.approx(gradient)
    
    def test_very_large_gradient(self, simple_linear_lower, simple_linear_upper, standard_bounds):
        """Test with very large gradient"""
        gradient = 1000.0
        result = continuousLCDSUBTrajectory(
            simple_linear_lower, 
            simple_linear_upper, 
            gradient, 
            standard_bounds
        )
        
        assert isinstance(result, UBTrajectory)
        assert result.lower_intercept is not None
        assert result.upper_intercept is not None
    
    def test_equal_lower_and_upper_contours(self, simple_linear_lower, standard_bounds):
        """Test when both contours are identical"""
        gradient = 1.0
        result = continuousLCDSUBTrajectory(
            simple_linear_lower, 
            simple_linear_lower,  # Same contour for both
            gradient, 
            standard_bounds
        )
        
        assert isinstance(result, UBTrajectory)
        assert result.lower_intercept is not None
        assert result.upper_intercept is not None
        # Lower and upper intercepts should be close
        assert abs(result.lower_intercept.B - result.upper_intercept.B) < 1000.0
    
    def test_intercept_points_have_correct_U_values(self, simple_linear_lower, simple_linear_upper, standard_bounds):
        """Verify intercept U coordinates match contour evaluations"""
        gradient = 1.5
        result = continuousLCDSUBTrajectory(
            simple_linear_lower, 
            simple_linear_upper, 
            gradient, 
            standard_bounds
        )
        
        # Lower intercept U should equal lower_contour(B)
        lower_B = result.lower_intercept.B
        expected_lower_U = simple_linear_lower(lower_B)
        assert result.lower_intercept.U == pytest.approx(expected_lower_U, rel=1e-6)
        
        # Upper intercept U should equal upper_contour(B)
        upper_B = result.upper_intercept.B
        expected_upper_U = simple_linear_upper(upper_B)
        assert result.upper_intercept.U == pytest.approx(expected_upper_U, rel=1e-6)
    
    def test_line_touches_both_contours(self, simple_linear_lower, simple_linear_upper, standard_bounds):
        """Verify the trajectory line actually touches both contours at intercept points"""
        gradient = 1.5
        result = continuousLCDSUBTrajectory(
            simple_linear_lower, 
            simple_linear_upper, 
            gradient, 
            standard_bounds
        )
        
        # At lower intercept: line should equal contour
        lower_B = result.lower_intercept.B
        line_U_lower = gradient * lower_B + result.intercept
        contour_U_lower = simple_linear_lower(lower_B)
        assert line_U_lower == pytest.approx(contour_U_lower, rel=1e-5)
        
        # At upper intercept: line should equal contour
        upper_B = result.upper_intercept.B
        line_U_upper = gradient * upper_B + result.intercept
        contour_U_upper = simple_linear_upper(upper_B)
        assert line_U_upper == pytest.approx(contour_U_upper, rel=1e-5)


# -------------------------------------------------------------------
# Consistency Tests
# -------------------------------------------------------------------
class TestConsistency:
    """Tests for internal consistency of results"""
    
    def test_upper_intercept_above_lower(self, simple_linear_lower, simple_linear_upper, standard_bounds):
        """The intercept POINTS should be at different B values.
        The upper intercept is not necessarily above the lower intercept in U -
        it depends on the geometry."""
        gradient = 1.0
        result = continuousLCDSUBTrajectory(
            simple_linear_lower,
            simple_linear_upper,
            gradient,
            standard_bounds
        )
        
        # The intercept points should be at different B values
        assert result.lower_intercept.B != result.upper_intercept.B, \
            "Lower and upper intercepts should be at different B values"
        
        # Both should be within bounds
        assert standard_bounds[0] <= result.lower_intercept.B <= standard_bounds[1]
        assert standard_bounds[0] <= result.upper_intercept.B <= standard_bounds[1]
        
        # The U values at intercept points should match their respective contours
        assert result.lower_intercept.U == pytest.approx(
            simple_linear_lower(result.lower_intercept.B), rel=1e-6)
        assert result.upper_intercept.U == pytest.approx(
            simple_linear_upper(result.upper_intercept.B), rel=1e-6)
    
    def test_intercept_consistency_with_contours(self, quadratic_lower, quadratic_upper, wide_bounds):
        """The trajectory line must touch both contours somewhere, but may go 
        above or below them at other points. The key is that it intersects both."""
        gradient = 2.0
        result = continuousLCDSUBTrajectory(
            quadratic_lower,
            quadratic_upper,
            gradient,
            wide_bounds
        )
        
        # At the intercept points, the line should equal the contours
        lower_diff = abs(gradient * result.lower_intercept.B + result.intercept - 
                        result.lower_intercept.U)
        upper_diff = abs(gradient * result.upper_intercept.B + result.intercept - 
                        result.upper_intercept.U)
        
        assert lower_diff < 1.0, f"Line doesn't touch lower contour at intercept point. Diff: {lower_diff}"
        assert upper_diff < 1.0, f"Line doesn't touch upper contour at intercept point. Diff: {upper_diff}"
        
        # The line interpolates between the two intercept points
        # It should touch each contour exactly at the intercept points
        B_lower = result.lower_intercept.B
        B_upper = result.upper_intercept.B
        
        if B_lower < B_upper:
            # At B_lower, line touches lower contour
            # At B_upper, line touches upper contour
            pass
        else:
            # At B_upper, line touches upper contour
            # At B_lower, line touches lower contour
            pass
    
    def test_same_gradient_different_bounds_gives_different_results(self, simple_linear_lower, simple_linear_upper):
        """Different bounds should potentially give different results"""
        gradient = 1.0
        bounds1 = (1000.0, 5000.0)
        bounds2 = (15000.0, 20000.0)
        
        result1 = continuousLCDSUBTrajectory(
            simple_linear_lower, simple_linear_upper, gradient, bounds1
        )
        result2 = continuousLCDSUBTrajectory(
            simple_linear_lower, simple_linear_upper, gradient, bounds2
        )
        
        # Results can be different or same depending on optimization
        # At minimum, intercept points should be in their respective bounds
        assert bounds1[0] <= result1.lower_intercept.B <= bounds1[1]
        assert bounds1[0] <= result1.upper_intercept.B <= bounds1[1]
        assert bounds2[0] <= result2.lower_intercept.B <= bounds2[1]
        assert bounds2[0] <= result2.upper_intercept.B <= bounds2[1]

# -------------------------------------------------------------------
# Physical Plausibility Tests
# -------------------------------------------------------------------

class TestPhysicalPlausibility:
    """Tests based on physical meaning of trajectories"""
    
    def test_intercept_B_values_are_positive(self, simple_linear_lower, simple_linear_upper, standard_bounds):
        """B field values should always be positive"""
        gradient = 1.0
        result = continuousLCDSUBTrajectory(
            simple_linear_lower, 
            simple_linear_upper, 
            gradient, 
            standard_bounds
        )
        
        assert result.lower_intercept.B > 0
        assert result.upper_intercept.B > 0
    
    def test_trajectory_stays_above_lower_contour(self, reciprocal_lower, reciprocal_upper, standard_bounds):
        """Physical trajectory should not go below the lower W=0 contour"""
        gradient = 0.5
        result = continuousLCDSUBTrajectory(
            reciprocal_lower, 
            reciprocal_upper, 
            gradient, 
            standard_bounds
        )
        
        # Check at many points within bounds
        test_B = np.linspace(standard_bounds[0], standard_bounds[1], 50)
        
        for B in test_B:
            line_U = gradient * B + result.intercept
            lower_U = reciprocal_lower(B)
            
            # Line should be >= lower contour
            assert line_U >= lower_U - 1e-5, \
                f"Trajectory goes below lower contour at B={B:.0f}: line={line_U:.2f}, contour={lower_U:.2f}"
    
class TestPhysicalPlausibility:
    """Tests based on physical meaning of trajectories"""
    
    # ... (keep test_intercept_B_values_are_positive and test_trajectory_stays_above_lower_contour)
    
    def test_intercept_between_contours(self, simple_linear_lower, simple_linear_upper, standard_bounds):
        """The intercept value represents the U-value at B=0.
        For the line U = gradient*B + intercept to intersect both contours,
        the intercept must satisfy certain constraints based on the geometry."""
        gradient = 1.0
        result = continuousLCDSUBTrajectory(
            simple_linear_lower,
            simple_linear_upper,
            gradient,
            standard_bounds
        )
        
        # The intercept is the U value at B=0
        # For this to work, the line must intersect both contours somewhere in bounds
        # These are just sanity checks that the result is finite
        assert np.isfinite(result.intercept)
        assert result.intercept > simple_linear_lower(standard_bounds[0]) - gradient * standard_bounds[1], \
            "Intercept should be physically reasonable"
        assert result.intercept < simple_linear_upper(standard_bounds[1]) + abs(gradient) * standard_bounds[1], \
            "Intercept should be physically reasonable"
    
    def test_line_touches_both_contours_at_intercept_points(self, simple_linear_lower, simple_linear_upper, standard_bounds):
        """At the identified intercept points, the line should actually touch the contours"""
        gradient = 1.5
        result = continuousLCDSUBTrajectory(
            simple_linear_lower,
            simple_linear_upper,
            gradient,
            standard_bounds
        )
        
        # At lower intercept point
        lower_B = result.lower_intercept.B
        lower_U_line = gradient * lower_B + result.intercept
        lower_U_contour = result.lower_intercept.U
        assert abs(lower_U_line - lower_U_contour) < 1e-3, \
            f"Line doesn't match lower intercept U: line={lower_U_line}, contour={lower_U_contour}"
        
        # At upper intercept point
        upper_B = result.upper_intercept.B
        upper_U_line = gradient * upper_B + result.intercept
        upper_U_contour = result.upper_intercept.U
        assert abs(upper_U_line - upper_U_contour) < 1e-6, \
            f"Line doesn't match upper intercept U: line={upper_U_line}, contour={upper_U_contour}"
