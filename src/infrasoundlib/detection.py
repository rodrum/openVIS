"""
This module defines the Detection class.
"""


from datetime import datetime, timezone
from typing import Any


class Detection:
    """
    A class to represent a Detection. It contains all the quantities that can be read from pmcc
    files. It can be a family (pmcc analysis) but also directly a pixel.

    Attributes
    ----------
    data_list : :obj:`list`
        Usually a splitted line of a bulletin.
    t_start : datetime
        Start time of detection.
    t_end : datetime
        End time of detection.
    f_min : float
        Minimum frequency of detection.
    f_max : float
        Maximum frequency of detection.
    amp : float
        Amplitude of detection.
    azi : float
        Azimuth of the station at which the detection was recorded.
    speed : float
        Speed of the detection.
    std_azi
    std_speed
    nb_pts : integer
        Number of pixels or number of families per nominal detection.
    cons : float
        Consistency of detection (family).
    cor : float
        Correlation of detection (family).
    f_mean : float
        Mean frequency of detection.
    q_tau
    delays
    sensors
    picking
    fisc_stat
    p_max
    pseudo_period
    ampmax : float
        Max amplitude of the detection (family)
    nb_sensors
    attenuation : float
        Attenuation of the detection.


    """

    # pylint: disable=too-many-instance-attributes
    # All those attributes can be used (or at least recorded) for a detection.

    def __init__(self, _data: list[Any]):
        """
        Constructor of the :class:`.Detection`. It loads the data and initialize
        all the quantities to None.

        Parameters
        ----------
        _data : :obj:`list`
             List in which the detection data is contained.

        """

        self.data_list = _data
        self.t_start: datetime | None = None
        self.t_end: datetime | None = None
        self.f_min: float | None = None
        self.f_max: float | None = None
        self.amp: float | None = None
        self.azi: float | None = None
        self.speed: float | None = None
        self.std_azi: float | None = None
        self.std_speed: float | None = None
        self.nb_pts: int | None = None
        self.cons: float | None = None
        self.cor: float | None = None
        self.f_mean: float | None = None
        self.q_tau: float | None = None
        self.delays: list[Any] = []
        self.sensors: list[Any] = []
        self.picking: float | None = None
        self.fisc_stat: float | None = None
        self.p_max: float | None = None
        self.pseudo_period: float | None = None
        self.ampmax: float | None = None
        self.nb_sensors: int | None = None
        self.attenuation: float | None = None

    @classmethod
    def from_db(cls, det_db) -> "Detection":  # type: ignore
        """
        Create a :class:`.Detection` object and load its information from a
        Pandas DataFrame object.

        Parameters
        ----------
        det_db : Panadas DataFrame
            Loaded from the local binary. Each 'row' of the form (e.g.):
                station_id                  I26DE
                t_start       2021-01-01 09:40:00
                t_end         2021-01-01 09:45:00
                f_mean                   2.202104
                nb_pts                       46.0
                amp                      0.001084
                amp_max                  0.001084
                azi                     292.55856
                speed                    0.334075
                Name: 0, dtype: object

        Returns
        -------
        :class:`.Detection`
            New detection object.

        """

        det = cls([])
        det.t_start = det_db['t_start'].replace(tzinfo=timezone.utc)
        det.t_end = det_db['t_end'].replace(tzinfo=timezone.utc)
        det.f_min = None
        det.f_max = None
        det.f_mean = det_db['f_mean']
        det.amp = det_db['amp']
        det.ampmax = det_db['amp_max']
        det.azi = det_db['azi']
        det.speed = det_db['speed']
        det.nb_pts = det_db['nb_pts']
        return det

    def get_tstart(self) -> datetime:
        """
        Get t_start of detection.

        Returns
        -------
        datetime
            Start time of detection.

        """
        if self.t_start is None:
            raise TypeError("Start time is not known (None)")
        return self.t_start
