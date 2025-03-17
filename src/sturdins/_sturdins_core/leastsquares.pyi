"""

Least Squares
=============

Least squares algorithms for initializing the GNSS-InertialNav filter.
"""

from __future__ import annotations
import numpy

__all__ = ["GnssPVT", "MUSIC", "PhasedArrayAttitude", "RangeAndRate", "Wahba"]

def GnssPVT(
    x: numpy.ndarray[numpy.float64[m, 1], numpy.ndarray.flags.writeable],
    P: numpy.ndarray[
        numpy.float64[m, n], numpy.ndarray.flags.writeable, numpy.ndarray.flags.f_contiguous
    ],
    sv_pos: numpy.ndarray[numpy.float64[3, n], numpy.ndarray.flags.f_contiguous],
    sv_vel: numpy.ndarray[numpy.float64[3, n], numpy.ndarray.flags.f_contiguous],
    psr: numpy.ndarray[numpy.float64[m, 1]],
    psrdot: numpy.ndarray[numpy.float64[m, 1]],
    psrvar: numpy.ndarray[numpy.float64[m, 1]],
    psrdotvar: numpy.ndarray[numpy.float64[m, 1]],
) -> bool:
    """
    GnssPVT
    =======

    Least Squares solver for GNSS position, velocity, and timing terms

    Parameters
    ----------

    x : np.ndarray

        Initial state estimate

    P : np.ndarray

        Initial covariance estimate

    sv_pos : np.ndarray

        Satellite ECEF positions [m]

    sv_vel : np.ndarray

        Satellite ECEF velocities [m/s]

    psr : np.ndarray

        Pseudorange measurements [m]

    psrdot : np.ndarray

        Pseudorange-rate measurements [m/s]

    psrvar : np.ndarray

        Pseudorange measurement variance [m^2]

    psrdotvar : np.ndarray

        Pseudorange-rate measurement variance [(m/s)^2]

    Returns
    -------

    status : bool

        Convergence success
    """

def MUSIC(
    az_mean: float,
    el_mean: float,
    P: numpy.ndarray[numpy.complex128[m, 1]],
    ant_xyz: numpy.ndarray[numpy.float64[3, n], numpy.ndarray.flags.f_contiguous],
    n_ant: int,
    lamb: float,
    thresh: float = 0.0001,
) -> None:
    """
    MUSIC
    ===================

    MUSIC estimator using Prompt correlators (I & Q)

    Parameters
    ----------

    az_mean : double

        Azimuth estimates [rad]

    el_mean : double

        Elevation estimates [rad]

    P : np.ndarray

        Measured prompt correlators

    ant_xyz : np.ndarray

        Known antenna positions in the body frame

    n_ant : np.ndarray

        Known number of antennas in the array

    lamb : np.ndarray

        Wavelength for the signal of interest [m/rad]

    thresh : np.ndarray

        Desired threshold of convergence

    Returns
    -------

    status : bool

        Convergence success
    """

def PhasedArrayAttitude(
    C_b_l: numpy.ndarray[
        numpy.float64[3, 3], numpy.ndarray.flags.writeable, numpy.ndarray.flags.f_contiguous
    ],
    u_ned: numpy.ndarray[numpy.float64[3, n], numpy.ndarray.flags.f_contiguous],
    phase: numpy.ndarray[numpy.float64[m, n], numpy.ndarray.flags.f_contiguous],
    phasevar: numpy.ndarray[numpy.float64[m, n], numpy.ndarray.flags.f_contiguous],
    ant_xyz: numpy.ndarray[numpy.float64[m, n], numpy.ndarray.flags.f_contiguous],
    n_ant: int,
    lamb: float,
    thresh: float = 1e-06,
) -> bool:
    """
    PhasedArrayAttitude
    ===================

    Iterative attitude estimate based on the known spatial phase of an antenna array

    Parameters
    ----------

    C_b_l : np.ndarray

        Initial estimate of the body to local-nav frame attitude dcm

    u_ned : np.ndarray

        3 x n_sv Ephemeris based unit vectors in the local-nav frame

    phase : np.ndarray

        n_ant x n_sv matrix of measured differential gnss phase values

    phasevar : np.ndarray

        Variance of each phase measurement

    ant_xyz : np.ndarray

        Known antenna positions in the body frame

    n_ant : np.ndarray

        Known number of antennas in the array

    lamb : np.ndarray

        Wavelength for the signal of interest [m/rad]

    thresh : np.ndarray

        Desired threshold of convergence

    Returns
    -------

    status : bool

        Convergence success
    """

def RangeAndRate(
    pos: numpy.ndarray[numpy.float64[3, 1]],
    vel: numpy.ndarray[numpy.float64[3, 1]],
    cb: float,
    cd: float,
    sv_pos: numpy.ndarray[numpy.float64[3, 1]],
    sv_vel: numpy.ndarray[numpy.float64[3, 1]],
    pred_u: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable],
    pred_udot: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable],
    pred_psr: float,
    pred_psrdot: float,
) -> None:
    """
    RangeAndRate
    ============

    Predicts a range and rate based on a satellite location and velocity

    Parameters
    ----------

    pos : np.ndarray

        3x1 User ECEF position [m]

    vel : np.ndarray

        3x1 User ECEF velocity [m/s]

    cb : double

        User clock bias [m]

    cd : double

        User clock drift [m/s]

    sv_pos : np.ndarray

        3x1 Satellite ECEF positions [m]

    sv_vel : np.ndarray

        3x1 Satellite ECEF velocities [m/s]

    pred_u : np.ndarray

        3x1 reference to unit vector to satellite

    pred_udot : np.ndarray

        3x1 reference to unit vector rate of change to satellite

    pred_psr : double

        Reference to pseudorange prediction

    pred_psrdot :double

        Reference to pseudorange-rate prediction
    """

def Wahba(
    C_b_l: numpy.ndarray[
        numpy.float64[3, 3], numpy.ndarray.flags.writeable, numpy.ndarray.flags.f_contiguous
    ],
    u_body: numpy.ndarray[numpy.float64[3, n], numpy.ndarray.flags.f_contiguous],
    u_ned: numpy.ndarray[numpy.float64[3, n], numpy.ndarray.flags.f_contiguous],
    u_body_var: numpy.ndarray[numpy.float64[m, 1]],
) -> None:
    """
    Whaba
    =====

    Wahba's problem solver

    Parameters
    ----------

    C_b_l : np.ndarray

        Attitude DCM estimate (local-nav to body)

    u_body : np.ndarray

        Measured unit vectors in the body frame

    u_ned : np.ndarray

        Ephemeris based unit vectors in the local-nav frame

    u_body_var : np.ndarray

        Variance of the measured unit vectors
    """
