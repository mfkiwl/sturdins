"""

Nav Sensors
===========

Generators for navigation sensor error parameters.

Currently Contains the following sensors:

1. `NavigationIMU`
2. `NavigationClock`

"""

from __future__ import annotations

__all__ = [
    "AUTOMOTIVE",
    "CESIUM",
    "CONSUMER",
    "GetNavClock",
    "GetNavImu",
    "HIGH_QUALITY_TCXO",
    "LOW_QUALTIY_TCXO",
    "NavImuToSiUnits",
    "NavigationClock",
    "NavigationIMU",
    "OCXO",
    "RUBIDIUM",
    "TACTICAL",
]

class NavigationClock:
    """

    NavigationClock
    ===============

    Struct of navigation clock Allan variance values

    """

    h0: float
    h1: float
    h2: float
    def __init__(self) -> None: ...

class NavigationIMU:
    """

    NavigationIMU
    =============

    Struct of navigation IMU Allan variance values

    """

    Ba: float
    Bg: float
    Ka: float
    Kg: float
    Na: float
    Ng: float
    Ta: float
    Tg: float
    def __init__(self) -> None: ...

def GetNavClock(clock_name: str) -> NavigationClock:
    """
    GetNavClock
    ===========

    Generator for the default navigation clocks

    Parameters
    ----------

    clock_name : str

        Name of the desired default clock

    Returns
    -------

    imu : NavigationClock

        Allan variance parameters for the desired clock
    """

def GetNavImu(imu_name: str) -> NavigationIMU:
    """
    GetNavImu
    =========

    Generator for the default navigation IMUs

    Parameters
    ----------

    imu_name : str

        Name of the desired default IMU

    Returns
    -------

    imu : NavigationIMU

        Allan variance parameters for the desired IMU
    """

def NavImuToSiUnits(imu: NavigationIMU) -> None:
    """
    NavImuToSiUnits
    ===============

    Converts the NavigationIMU parameters to SI units of a navigation filter

    Parameters
    ----------

    imu : NavigationIMU

        Allan variance parameters for the IMU
    """

AUTOMOTIVE: NavigationIMU  # value = <sturdins._sturdins_core.navsense.NavigationIMU object>
CESIUM: NavigationClock  # value = <sturdins._sturdins_core.navsense.NavigationClock object>
CONSUMER: NavigationIMU  # value = <sturdins._sturdins_core.navsense.NavigationIMU object>
HIGH_QUALITY_TCXO: (
    NavigationClock  # value = <sturdins._sturdins_core.navsense.NavigationClock object>
)
LOW_QUALTIY_TCXO: (
    NavigationClock  # value = <sturdins._sturdins_core.navsense.NavigationClock object>
)
OCXO: NavigationClock  # value = <sturdins._sturdins_core.navsense.NavigationClock object>
RUBIDIUM: NavigationClock  # value = <sturdins._sturdins_core.navsense.NavigationClock object>
TACTICAL: NavigationIMU  # value = <sturdins._sturdins_core.navsense.NavigationIMU object>
