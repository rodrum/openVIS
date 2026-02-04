"""
This module defines the Station class.
"""
from src.infrasoundlib.detection import Detection
from src.infrasoundlib.vratio import Vratio


class Station:
    """
    A class to represent a station. It contains its position information,
    its detections and vratio.
    Parent of class :class:`.StationVolc`.

    Attributes
    ----------
    name : str
        Name of the station.
    lat : float
        Latitude of the station.
    lon : float
        Longitude of the station.
    v_ratio : :obj:`list`
        List containing the vratio data.
    detections : :obj:`list` of :class:`.Detection`
        List of detections of the station.
    nb_dets : int
        Number of detections of the station.
    az_avg : :obj:`list` of :obj:`float`
        List of average azimuth (unused)

    """

    # pylint: disable=too-many-arguments

    def __init__(self, name: str, lat: float, lon: float):
        """
        Constructor of :class:`.Station`. Initialize the name and position and set other attributes
        to None or empty :obj:`list`.

        Parameters
        ----------
        name : str
            Name of the station.
        lat : float
            Latitude of the station.
        lon : float
            Longitude of the station.

        """

        self.name = name
        self.lat = lat
        self.lon = lon
        self.v_ratio: list[Vratio] | None = None
        self.detections: list[Detection] = []
        self.nb_dets: int | None = None
