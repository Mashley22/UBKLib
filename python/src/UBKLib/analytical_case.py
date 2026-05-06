import numpy as np
from numpy.polynomial import Polynomial
import matplotlib.pyplot as plt

from .field_models import cross_tail, dipole
from .types import Vectorizable, UBCoord, UBTrajectory
from .trajectory import continuous_lcds_ub_trajectory

'''
This case is for the cross tail potential combined with a dipole 
magnetic field at K=0, the same case as supplied in Whipple.
In this case the U-B curves can be found analytically and as such the 
trajectories in U-B space can also be found analytically, the extrema points 
in U-B space can be found analytically as well.
'''
class UBCurve:
    surface_potential: float
    cross_tail_strength: float
    dipole_strength: float
    inner_lim: float
    outer_lim: float

    def __init__(self,
                 surface_potential: float = cross_tail.SURFACE_POTENTIAL,
                 cross_tail_strength: float = cross_tail.CROSS_TAIL_STRENGTH,
                 dipole_strength: float = dipole.SURFACE_STRENGTH,
                 inner_lim: float = 1.05,
                 outer_lim: float = 15.0):
                         
        self.surface_potential = surface_potential
        self.cross_tail_strength = cross_tail_strength
        self.dipole_strength = dipole_strength
        self.inner_lim = inner_lim
        self.outer_lim = outer_lim

    def upper_contour(self,
                      magnetic_amplitude: Vectorizable) -> Vectorizable:
        '''
        Calculates the value of the electric potential upper value for a given
        magnetic field amplitude
        
        Parameters:
            magnetic_amplitude: in nT

        Returns:
            electric potential (kV)
        '''
        return self.cross_tail_strength * np.pow(self.dipole_strength / magnetic_amplitude, 1/3) \
            - self.surface_potential * np.pow(magnetic_amplitude / self.dipole_strength, 1/3)

    def lower_contour(self,
                      magnetic_amplitude: Vectorizable) -> Vectorizable:
        '''
        Calculates the value of the electric potential lower value for a given
        magnetic field amplitude
        
        Parameters:
            magnetic_amplitude (nT)

        Returns:
            electric potential (kV)
        '''
        return -self.cross_tail_strength * np.pow(self.dipole_strength / magnetic_amplitude, 1/3) \
            - self.surface_potential * np.pow(magnetic_amplitude / self.dipole_strength, 1/3)

    def radius(self,
               magnetic_amplitude: Vectorizable) -> Vectorizable:
        '''
        Calculates the value of the radius for a given value of the magnetic field amplitude

        Parameters:
            magnetic_amplitude (nT)

        Returns:
            radius (Re)
        '''

        return np.pow(self.dipole_strength / magnetic_amplitude, 1/3)

    def lcds_line(self, gradient: float) -> UBTrajectory:
        '''
        Calculates the UB coordinates of the LCDS trajectory in the UB space where this trajectory
        meets the W=0 contours

        Parameters:
            gradient (kV/nT): the value of -mu/q

        Returns:
            B1, U1, B2, U2 where B1, U1 refer to the point where this LCDS trajectory meets the lower curve.
            Likewise B2, U2, refer to upper curve meeting point.
        '''

        if gradient >= 0:
            return self.__lower_tangent_lcds_line(gradient)
        else:
            return continuous_lcds_ub_trajectory(
                self.lower_contour,
                self.upper_contour,
                gradient,
                (
                    dipole.equitorial_dipole_amplitude(self.outer_lim, 0, self.dipole_strength),
                    dipole.equitorial_dipole_amplitude(self.inner_lim, 0, self.dipole_strength),
                )
            )

    def __lower_tangent_lcds_line(self, gradient: float) -> UBTrajectory:
        '''
        Compute the lcds meeting points with the W=0 contour if the lcds trajectory in UB space
        is taken to be tangent to the lower curve
        '''
        trajectory = UBTrajectory()
        trajectory.gradient = gradient
        trajectory.lower_intercept = self.__lcds_lower_intercept(gradient)

        trajectory.intercept = U1 - gradient * B1

        upper_sol = Polynomial([
            -self.cross_tail_strength * np.pow(self.dipole_strength, 1/3),
            trajectory.intercept,
            self.surface_potential / np.pow(self.dipole_strength, 1/3),
            0,
            gradient
        ])

        roots = upper_sol.roots() 
        valid_roots = roots[(roots.imag == 0) & (roots.real > 0)].real

        assert valid_roots.size == 1, "Implementation error"

        B2 = np.pow(valid_roots[0], 3)
        trajectory.upper_intercept = UBCoord(
                B=B2,
                U=self.upper_contour(B2)
        )

        return trajectory 

    def __lcds_lower_intercept(self, gradient: float) -> UBCoord:

        a = self.cross_tail_strength * np.pow(self.dipole_strength, 1/3) / 3
        b = - self.surface_potential /(3 * np.pow(self.dipole_strength, 1/3))
        c = -gradient

        x = (-b + np.sqrt(b * b - 4 * a * c)) / (2 * a)

        B1 = np.pow(x, -3/2)
        return UBCoord(B=B1,
            U=self.lower_contour(B1)
        )
