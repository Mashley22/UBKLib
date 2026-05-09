import numpy as np

from .field_models import cross_tail, dipole
from .types import Vectorizable, UBTrajectory
from .trajectory import continuous_lcds_ub_trajectory


class UBCurve:
    '''
    This case is for the cross tail potential combined with a dipole 
    magnetic field at K=0, the same case as supplied in Whipple.
    In this case the U-B curves can be found analytically and as such the 
    trajectories in U-B space can also be found analytically, the extrema points 
    in U-B space can be found analytically as well.
    '''

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
        return self.cross_tail_strength * np.pow(self.dipole_strength / magnetic_amplitude, 1 / 3) \
            - self.surface_potential * np.pow(magnetic_amplitude / self.dipole_strength, 1 / 3)

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
        return -self.cross_tail_strength * np.pow(self.dipole_strength / magnetic_amplitude, 1 / 3) \
            - self.surface_potential * np.pow(magnetic_amplitude / self.dipole_strength, 1 / 3)

    def radius(self,
               magnetic_amplitude: Vectorizable) -> Vectorizable:
        '''
        Calculates the value of the radius for a given value of the magnetic field amplitude

        Parameters:
            magnetic_amplitude (nT)

        Returns:
            radius (Re)
        '''

        return np.pow(self.dipole_strength / magnetic_amplitude, 1 / 3)

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

        return continuous_lcds_ub_trajectory(
            self.lower_contour,
            self.upper_contour,
            gradient,
            (
                dipole.equatorial_dipole_amplitude(self.outer_lim, 0, self.dipole_strength),
                dipole.equatorial_dipole_amplitude(self.inner_lim, 0, self.dipole_strength),
            )
        )
