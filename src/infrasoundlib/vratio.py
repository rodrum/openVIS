"""
Vratio class.
"""

from dataclasses import dataclass
from datetime import datetime


@dataclass
class Vratio:
    """
    A class to represent a Vratio computation.

    Attributes
    ----------
    dt : datetime
        Time
    azimuth : int
        Azimtuth for the vratio
    value : float
        Value

    """

    dt: datetime
    azimuth: int
    value: float
