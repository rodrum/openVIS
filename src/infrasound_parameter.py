"""
This module defines the IP class.
"""

from datetime import datetime


class IP:
    """
    A class to represent IP values.

    The Infrasound Parameter is used to characterize eruptions.
    If its value exceed a certain threshold; the VIS will emit a new eruption notification.

    Attributes
    ----------
    dt : datetime
        Date at which this IP was recorded.
    value : float
        Value of the IP.
    nb_det : int
        Number of detections that were recorded during this IP time interval.
    avg_amp : float
        mean amplitude of the detections during this IP time interval.
    fmaean : float
        Mean frequency of the detections during the IP time interval.
    vratio : float
        Value of the vratio for this IP calculation.

    """

    def __init__(
        self,
        dt: datetime,
        val: float = 0,
        nb_det: int = 0,
        avg_amp: float | None = None,
        avg_source_amp: float | None = None,
        persistency: float | None = None,
        fmean: float | None = None,
        vratio: float | None = None,
    ):
        """
        Constructor for the class :class:`.IP`.

        Parameters
        ----------
        dt : datetime
            Date at which this IP was recorded.
        value : float
            Value of the IP.
        nb_det : int
            Number of detections that were recorded during this IP time interval.
        avg_amp : float
            Average amplitude of the detections during this IP time interval.
        fmaean : float
            Mean frequency of the detections during the IP time interval.
        vratio : float
            Value of the vratio for this IP calculation.

        """

        self.dt = dt
        self.value = val
        self.nb_det = nb_det
        self.avg_amp = avg_amp
        self.fmean = fmean
        self.vratio = vratio
        self.persistency = persistency
        self.avg_source_amp = avg_source_amp
