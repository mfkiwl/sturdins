"""

SturdINS
========

Daniel Sturdivant's GNSS-INS Sofware.

Contains the following classes:

1. `InertialNav`
2. `KinematicNav`
3. `Strapdown`

Contains the following modules:

1. `leastsquares`
2. `navsense`


"""

from __future__ import annotations
import numpy
import typing
from . import leastsquares
from . import navsense

__all__ = ["InertialNav", "KinematicNav", "Strapdown", "leastsquares", "navsense"]

class InertialNav:
    """

    InertialNav
    ===

    Inertial navigation Kalman Filter equations.

    """

    C_b_l_: numpy.ndarray[numpy.float64[3, 3]]
    P_: numpy.ndarray[numpy.float64[m, n]]
    ba_: numpy.ndarray[numpy.float64[3, 1]]
    bg_: numpy.ndarray[numpy.float64[3, 1]]
    cb_: float
    cd_: float
    h_: float
    lam_: float
    phi_: float
    q_b_l_: numpy.ndarray[numpy.float64[4, 1]]
    vd_: float
    ve_: float
    vn_: float
    def GnssUpdate(
        self,
        sv_pos: numpy.ndarray[numpy.float64[m, n], numpy.ndarray.flags.f_contiguous],
        sv_vel: numpy.ndarray[numpy.float64[m, n], numpy.ndarray.flags.f_contiguous],
        psr: numpy.ndarray[numpy.float64[m, 1]],
        psrdot: numpy.ndarray[numpy.float64[m, 1]],
        psrvar: numpy.ndarray[numpy.float64[m, 1]],
        psrdotvar: numpy.ndarray[numpy.float64[m, 1]],
    ) -> None:
        """
        GnssUpdate
        ==========

        Correct state with GPS measurements

        Parameters
        ----------

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
        """

    def Mechanize(
        self,
        w_ib_b: numpy.ndarray[numpy.float64[3, 1]],
        f_ib_b: numpy.ndarray[numpy.float64[3, 1]],
        dt: float,
    ) -> None:
        """
        Mechanize
        =========

        Integrate measured angular rates (delta thetas) and specific forces (delta velocities)

        Parameters
        ----------

        w_ib_b : np.ndarray

            Measured angular rates (delta thetas) in the body frame [rad/s]

        f_ib_b : np.ndarray

            Measured specific forces (delta velocities) in the body frame [m/s^2]

        dt : double

            Integration time [s]
        """

    def PhasedArrayUpdate(
        self,
        sv_pos: numpy.ndarray[numpy.float64[3, n], numpy.ndarray.flags.f_contiguous],
        sv_vel: numpy.ndarray[numpy.float64[3, n], numpy.ndarray.flags.f_contiguous],
        psr: numpy.ndarray[numpy.float64[m, 1]],
        psrdot: numpy.ndarray[numpy.float64[m, 1]],
        phase: numpy.ndarray[numpy.float64[m, n], numpy.ndarray.flags.f_contiguous],
        psrvar: numpy.ndarray[numpy.float64[m, 1]],
        psrdotvar: numpy.ndarray[numpy.float64[m, 1]],
        phasevar: numpy.ndarray[numpy.float64[m, n], numpy.ndarray.flags.f_contiguous],
        ant_xyz: numpy.ndarray[numpy.float64[3, n], numpy.ndarray.flags.f_contiguous],
        n_ant: int,
        lamb: float,
    ) -> None:
        """
        PhasedArrayUpdate
        =================

        Correct state with GPS measurements

        Parameters
        ----------

        sv_pos : np.ndarray

            Satellite ECEF positions [m]

        sv_vel : np.ndarray

            Satellite ECEF velocities [m/s]

        psr : np.ndarray

            Pseudorange measurements [m]

        psrdot : np.ndarray

            Pseudorange-rate measurements [m/s]

        phase : np.ndarray

            Delta phase measurements [m/s]

        psrvar : np.ndarray

            Pseudorange measurement variance [m^2]

        psrdotvar : np.ndarray

            Pseudorange-rate measurement variance [(m/s)^2]

        phasevar : np.ndarray

            Variance of each phase measurement

        ant_xyz : np.ndarray

            Known antenna positions in the body frame

        n_ant : np.ndarray

            Known number of antennas in the array

        lamb : np.ndarray

            Wavelength for the signal of interest [m/rad]
        """

    def Propagate(
        self,
        w_ib_b: numpy.ndarray[numpy.float64[3, 1]],
        f_ib_b: numpy.ndarray[numpy.float64[3, 1]],
        dt: float,
    ) -> None:
        """
        Propagate
        =========

        Propagate the error state matrices

        Parameters
        ----------

        w_ib_b : np.ndarray

            Measured angular rates (delta thetas) in the body frame [rad/s]

        f_ib_b : np.ndarray

            Measured specific forces (delta velocities) in the body frame [m/s^2]

        dt : double

            Integration time [s]
        """

    @typing.overload
    def SetAttitude(
        self, C: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.f_contiguous]
    ) -> None:
        """
        SetAttitude
        ===========

        Set the attitude of the navigator

        Parameters
        ----------

        C : np.ndarray

            body-to-ned rotation matrix
        """

    @typing.overload
    def SetAttitude(self, r: float, p: float, y: float) -> None:
        """
        SetAttitude
        ===========

        Set the attitude of the navigator

        Parameters
        ----------

        r : double

            roll angle [rad]

        p : double

            pitch angle [rad]

        y : double

            yaw angle [rad]
        """

    def SetClock(self, cb: float, cd: float) -> None:
        """
        SetClock
        ========

        Set the clock states of the navigator

        Parameters
        ----------

        cb : double

            Clock bias [m]

        cd : double

            Clock drift [m/s]
        """

    def SetClockSpec(self, h0: float, h1: float, h2: float) -> None:
        """
        SetClockSpec
        ============

        Set the noise parameters of the Clock

        Parameters
        ----------

        h0 : double

            white frequency modulation

        h1 : double

            flicker frequency modulation

        h2 : double

            random walk frequency modulation
        """

    def SetImuSpec(self, Ba: float, Na: float, Bg: float, Ng: float) -> None:
        """
        SetImuSpec
        ==========

        Set the noise parameters of the IMU

        Parameters
        ----------

        Ba : double

            Accelerometer bias InertialNavtability [mg]

        Na : double

            Accelerometer random walk [m/s/sqrt(hr)]

        Bg : double

            Gyroscope bias InertialNavtability [deg/hr]

        Ng : double

            Gyroscope random walk [deg/sqrt(hr)]
        """

    def SetPosition(self, lat: float, lon: float, alt: float) -> None:
        """
        SetPosition
        ===========

        Set the position of the navigator

        Parameters
        ----------

        lat : double

            Latitude [rad]

        lon : double

            Longitude [rad]

        alt : double

            Altitude/Height [m]
        """

    def SetVelocity(self, vn: float, ve: float, vd: float) -> None:
        """
        SetVelocity
        ===========

        Set the velocity of the navigator

        Parameters
        ----------

        vn : double

            north velocity [m/s]

        ve : double

            east velocity [m/s]

        vd : double

            down velocity [m/s]
        """

    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(
        self,
        lat: float,
        lon: float,
        alt: float,
        vn: float,
        ve: float,
        vd: float,
        r: float,
        p: float,
        y: float,
        cb: float,
        cd: float,
    ) -> None: ...

class KinematicNav:
    """

    KinematicNav
    ===

    KinematicNav navigation Kalman Filter equations.

    """

    C_b_l_: numpy.ndarray[numpy.float64[3, 3]]
    P_: numpy.ndarray[numpy.float64[m, n]]
    cb_: float
    cd_: float
    h_: float
    lam_: float
    phi_: float
    q_b_l_: numpy.ndarray[numpy.float64[4, 1]]
    vd_: float
    ve_: float
    vn_: float
    def GnssUpdate(
        self,
        sv_pos: numpy.ndarray[numpy.float64[3, n], numpy.ndarray.flags.f_contiguous],
        sv_vel: numpy.ndarray[numpy.float64[3, n], numpy.ndarray.flags.f_contiguous],
        psr: numpy.ndarray[numpy.float64[m, 1]],
        psrdot: numpy.ndarray[numpy.float64[m, 1]],
        psrvar: numpy.ndarray[numpy.float64[m, 1]],
        psrdotvar: numpy.ndarray[numpy.float64[m, 1]],
    ) -> None:
        """
        GnssUpdate
        ==========

        Correct state with GPS measurements

        Parameters
        ----------

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
        """

    def PhasedArrayUpdate(
        self,
        sv_pos: numpy.ndarray[numpy.float64[3, n], numpy.ndarray.flags.f_contiguous],
        sv_vel: numpy.ndarray[numpy.float64[3, n], numpy.ndarray.flags.f_contiguous],
        psr: numpy.ndarray[numpy.float64[m, 1]],
        psrdot: numpy.ndarray[numpy.float64[m, 1]],
        phase: numpy.ndarray[numpy.float64[m, n], numpy.ndarray.flags.f_contiguous],
        psrvar: numpy.ndarray[numpy.float64[m, 1]],
        psrdotvar: numpy.ndarray[numpy.float64[m, 1]],
        phasevar: numpy.ndarray[numpy.float64[m, n], numpy.ndarray.flags.f_contiguous],
        ant_xyz: numpy.ndarray[numpy.float64[3, n], numpy.ndarray.flags.f_contiguous],
        n_ant: int,
        lamb: float,
    ) -> None:
        """
        PhasedArrayUpdate
        =================

        Correct state with GPS measurements

        Parameters
        ----------

        sv_pos : np.ndarray

            Satellite ECEF positions [m]

        sv_vel : np.ndarray

            Satellite ECEF velocities [m/s]

        psr : np.ndarray

            Pseudorange measurements [m]

        psrdot : np.ndarray

            Pseudorange-rate measurements [m/s]

        phase : np.ndarray

            Delta phase measurements [m/s]

        psrvar : np.ndarray

            Pseudorange measurement variance [m^2]

        psrdotvar : np.ndarray

            Pseudorange-rate measurement variance [(m/s)^2]

        phasevar : np.ndarray

            Variance of each phase measurement

        ant_xyz : np.ndarray

            Known antenna positions in the body frame

        n_ant : np.ndarray

            Known number of antennas in the array

        lamb : np.ndarray

            Wavelength for the signal of interest [m/rad]
        """

    def Propagate(self, dt: float) -> None:
        """
        Propagate
        =========

        Propagate the error state matrices

        Parameters
        ----------

        dt : double

            Integration time [s]
        """

    @typing.overload
    def SetAttitude(
        self, C: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.f_contiguous]
    ) -> None:
        """
        SetAttitude
        ===========

        Set the attitude of the navigator

        Parameters
        ----------

        C : np.ndarray

            body-to-ned rotation matrix
        """

    @typing.overload
    def SetAttitude(self, r: float, p: float, y: float) -> None:
        """
        SetAttitude
        ===========

        Set the attitude of the navigator

        Parameters
        ----------

        r : double

            roll angle [rad]

        p : double

            pitch angle [rad]

        y : double

            yaw angle [rad]
        """

    def SetClock(self, cb: float, cd: float) -> None:
        """
        SetClock
        ========

        Set the clock states of the navigator

        Parameters
        ----------

        cb : double

            Clock bias [m]

        cd : double

            Clock drift [m/s]
        """

    def SetClockSpec(self, h0: float, h1: float, h2: float) -> None:
        """
        SetClockSpec
        ============

        Set the noise parameters of the Clock

        Parameters
        ----------

        h0 : double

            white frequency modulation

        h1 : double

            flicker frequency modulation

        h2 : double

            random walk frequency modulation
        """

    def SetPosition(self, lat: float, lon: float, alt: float) -> None:
        """
        SetPosition
        ===========

        Set the position of the navigator

        Parameters
        ----------

        lat : double

            Latitude [rad]

        lon : double

            Longitude [rad]

        alt : double

            Altitude/Height [m]
        """

    def SetProcessNoise(self, Svel: float, Satt: float) -> None:
        """
        SetProcessNoise
        ===============

        Set the noise parameters of the KinematicNav (constant velocity) model

        Parameters
        ----------

        Svel : double

            PSD of expected acceleration white noise [(m/s^2)^2]

        Satt : double

            PSD of expected angular rate white noise [(rad/s)^2]
        """

    def SetVelocity(self, vn: float, ve: float, vd: float) -> None:
        """
        SetVelocity
        ===========

        Set the velocity of the navigator

        Parameters
        ----------

        vn : double

            north velocity [m/s]

        ve : double

            east velocity [m/s]

        vd : double

            down velocity [m/s]
        """

    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(
        self,
        lat: float,
        lon: float,
        alt: float,
        vn: float,
        ve: float,
        vd: float,
        cb: float,
        cd: float,
    ) -> None: ...
    @typing.overload
    def __init__(
        self,
        lat: float,
        lon: float,
        alt: float,
        vn: float,
        ve: float,
        vd: float,
        r: float,
        p: float,
        y: float,
        cb: float,
        cd: float,
    ) -> None: ...

class Strapdown:
    """

    Strapdown
    =========

    Inertial navigation strapdown integration equations.

    """

    C_b_l_: numpy.ndarray[numpy.float64[3, 3]]
    h_: float
    lam_: float
    phi_: float
    q_b_l_: numpy.ndarray[numpy.float64[4, 1]]
    vd_: float
    ve_: float
    vn_: float
    def Mechanize(
        self,
        w_ib_b: numpy.ndarray[numpy.float64[3, 1]],
        f_ib_b: numpy.ndarray[numpy.float64[3, 1]],
        dt: float,
    ) -> None:
        """
        Mechanize
        =========

        Integrate measured angular rates (delta thetas) and specific forces (delta velocities)

        Parameters
        ----------

        w_ib_b : np.ndarray

            Measured angular rates (delta thetas) in the body frame [rad/s]

        f_ib_b : np.ndarray

            Measured specific forces (delta velocities) in the body frame [m/s^2]

        dt : double

            Integration time [s]
        """

    @typing.overload
    def SetAttitude(
        self, C: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.f_contiguous]
    ) -> None:
        """
        SetAttitude
        ===========

        Set the attitude of the navigator

        Parameters
        ----------

        C : np.ndarray

            body-to-ned rotation matrix
        """

    @typing.overload
    def SetAttitude(self, r: float, p: float, y: float) -> None:
        """
        SetAttitude
        ===========

        Set the attitude of the navigator

        Parameters
        ----------

        r : double

            roll angle [rad]

        p : double

            pitch angle [rad]

        y : double

            yaw angle [rad]
        """

    def SetPosition(self, lat: float, lon: float, alt: float) -> None:
        """
        SetPosition
        ===========

        Set the position of the navigator

        Parameters
        ----------

        lat : double

            Latitude [rad]

        lon : double

            Longitude [rad]

        alt : double

            Altitude/Height [m]
        """

    def SetVelocity(self, vn: float, ve: float, vd: float) -> None:
        """
        SetVelocity
        ===========

        Set the velocity of the navigator

        Parameters
        ----------

        vn : double

            north velocity [m/s]

        ve : double

            east velocity [m/s]

        vd : double

            down velocity [m/s]
        """

    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(
        self,
        lat: float,
        lon: float,
        alt: float,
        vn: float,
        ve: float,
        vd: float,
        r: float,
        p: float,
        y: float,
    ) -> None: ...

__version__: str = "1.0.0"
