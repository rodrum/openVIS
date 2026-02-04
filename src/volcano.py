"""
This module defines the Volcano class.
"""


class Volcano:
    """
    A class to represent a volcano.

    Attributes
    ----------
    code : int
        Official code of the volcano.
    name : str
        Name of the volcano.
    lat : float
        Latitude of the volcano.
    lon : float
        Longitude of the volcano.
    alt : float
        Altitude of the summit of the volcano.

    """

    def __init__(self, code: int, name: str, lat: float, lon: float, alt: float):
        """
        Constructor of :class:`.Volcano`.

        Parameters
        ----------
        code : int
            Official code of the volcano.
        name : str
            Name of the volcano.
        lat : float
            Latitude of the volcano.
        lon : float
            Longitude of the volcano.
        alt : float
            Altitude of the summit of the volcano.

        """

        self.code = int(code)
        self.name = name
        self.lat = float(lat)
        self.lon = float(lon)
        self.alt = float(alt)
