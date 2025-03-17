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
from sturdins._sturdins_core import InertialNav
from sturdins._sturdins_core import KinematicNav
from sturdins._sturdins_core import Strapdown
from sturdins._sturdins_core import leastsquares
from sturdins._sturdins_core import navsense
from . import _sturdins_core

__all__: list = [
    "__doc__",
    "__version__",
    "InertialNav",
    "KinematicNav",
    "Strapdown",
    "leastsquares",
    "navsense",
]
__version__: str = "1.0.0"
