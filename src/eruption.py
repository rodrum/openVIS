"""
This module defines the Eruption class.
"""

from datetime import datetime, timedelta

import numpy as np

from src.region import Region
from src.settings import STA_CLOSE_CONFIDENCE
from src.station_volc import StationVolc


class Eruption:
    """
    A class to represent an Eruption. Eruption are linked to a volcano and also have the data
    from the stations detecting.

    Attributes
    ----------
    volcanoes : :obj:`list` of :class:`.Volcano`
        Volcanoes at which the eruption takes place.
    t_start : datetime
        Start date of the eruption.
    t_end : datetime
        End date of the eruption.
    er_code : str
        Code of the eruption.
    data_stations : :obj:`list` of :class:`.Station`
        Data of the stations that recorded the eruption.
    amp_source : float
        Estimated source amplitude of the eruption.
    status : str
        Status of the eruption (ONGOING or ENDED or UNCONFIRMED).
    last_notification : datetime
        Date at which the last notification was emmited.
    revision : int
        Number of revision of the eruption.
    confidence_level : int
        Confidence level of the eruption.
        - 1: the eruption is uncertain, it means it was detected by only one
             station at more than 500 km.
        - 2: the eruption is probable, it was detected by two station at more
             that 500 km.
        - 3: the eruption is certain, it was either detected by three stations
             or by one station at less than 500 km.

    """

    # pylint: disable=too-many-instance-attributes

    def __init__(
        self,
        region: Region,
        stations: list[StationVolc],
        thresh_ip: float = 0,
        ip_interval: int = 0,
    ):
        """
        Init :class:`.Eruption` class from region and stations. Also perform the first update of
        the eruption.

        Parameters
        ----------
        region : :class:`.Region`
            Region in which the eruption takes place.
        stations : :obj:`list` of :class:`.Station`
            List of station recording the eruption.
        thresh_ip : float
            IP threshold.
        ip_interval : int
            Number of seconds of the IP calculation window.

        """

        self.region = region
        self.t_start: datetime | None = None
        self.t_end: datetime | None = None
        self.er_code = ""
        self.data_stations: list[StationVolc] = []
        self.amp_source = 0.0  # Pa
        self.status = "ONGOING"
        self.last_notification: datetime | None = None
        self.detection_date: datetime | None = None
        self.revision = 0
        self.confidence_level = 0
        self.update_eruption(stations, thresh_ip, ip_interval)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Eruption):
            raise NotImplementedError
        return (
            self.t_start == other.t_start
            and self.t_end == other.t_end
            and self.detection_date == other.detection_date
            and self.last_notification == other.last_notification
            and self.revision == other.revision
            and self.amp_source == other.amp_source
            and self.confidence_level == other.confidence_level
        )

    def generate_eruption_code(self) -> str:
        """
        Generate a code for the eruption.

        """

        if self.t_start:
            return (
                "Erupt"
                + self.create_julian_date(self.t_start)
                + "".join(sorted([v.name[:4].upper() for v in self.region.volcanoes]))
            )
        return ""

    def set_t_start(self, t_start: datetime) -> None:
        """
        Set the start time and generate a code fur the eruption.

        Parameters
        ----------
        t_start : datetime
            Start date of the eruption.

        """

        self.t_start = t_start
        self.er_code = self.generate_eruption_code()

    @staticmethod
    def create_julian_date(_date: datetime) -> str:
        """
        Static method creating a julian date from a datetime object.

        Parameters
        ----------
        _date : datetime
            Date used to generate the julian date.

        """

        time_tup = _date.timetuple()
        julian_date = str(time_tup.tm_year)
        if time_tup.tm_yday < 100:
            julian_date += "0"
            if time_tup.tm_yday < 10:
                julian_date += "0"

        return julian_date + str(time_tup.tm_yday)

    # TODO : Improve source term estimation
    def estimate_source_term(self) -> None:
        """
        Estimate the source amplitude of the eruption.

        The estimated amplitude is the maximum amplitude ever recorded at one of
        the detecting stations.

        """

        for sta in self.data_stations:
            if sta.detections:
                lst_amp = [d.amp for d in sta.detections]
                max_amp: float = np.max(lst_amp)
                if sta.max_amplitude < max_amp:
                    max_amp_det = sta.detections[lst_amp.index(max_amp)]
                    sta.max_amplitude = max_amp
                    if max_amp_det.attenuation is None:
                        raise TypeError("Attenuation of max amp det is unknown (None)")
                    sta.estimated_amplitude = (
                        sta.max_amplitude / max_amp_det.attenuation
                    )

                if sta.estimated_amplitude > self.amp_source:
                    self.amp_source = sta.estimated_amplitude

    def calculate_confidence_level(self) -> None:
        """
        Compute the confidence level of the eruption.

        An eruption is certain if at least 3 stations see it or one station sees it at less than a
        given distance (500 km by defaults).
        An eruption is probable if at least 2 stations see it at more than the given distance.
        An eruption is uncertain if only 1 station sees it at more than the given distance.

        """

        stations_detecting = [
            sta for sta in self.data_stations if sta.detecting is True
        ]
        sta_close = [
            sta
            for sta in stations_detecting
            if sta.distaz["dkm"] < STA_CLOSE_CONFIDENCE
        ]
        if sta_close or len(stations_detecting) >= 3:
            self.confidence_level = 3
        elif len(stations_detecting) == 2:
            self.confidence_level = 2
        elif len(stations_detecting) == 1:
            self.confidence_level = 1
        else:
            self.confidence_level = 0

    def update_eruption(
        self, stations: list[StationVolc], thresh_ip: float, ip_inter: int
    ) -> None:
        """
        Update the eruption infos based on the new data from the stations.

        Eruption update is done every analysis with the new data acquired. The end date is obvioulsy
        modified but also the start date in case there is a new station (located further away) that
        gets detection corresponding to an earlier start time that what was assessed from previous
        data. Also check which stations are detecting the eruption and estimate the source
        amplitude and confidence level.

        Parameters
        ----------
        stations : :obj:`list` of :class:`.Station`
            List of stations used for analysis.
        thresh_ip : float
            IP threshold.
        ip_inter : int
            Number of seconds of the IP calculation window.

        """

        for sta in stations:
            if sta.detections:
                if not sta.detecting:
                    tmp_ip = [
                        ip
                        for ip in sta.list_ip
                        if ip.value >= thresh_ip
                        and ip.dt  # type: ignore # mypy think comparison is datetime and timedelta
                        >= (
                            sta.detections[0].t_start - timedelta(seconds=ip_inter)  # type: ignore
                        )
                    ]
                    if tmp_ip:
                        sta.detecting = True
                if sta.detecting:
                    if (
                        sta.detections[0].t_start is None
                        or sta.detections[-1].t_end is None
                    ):
                        raise TypeError("Detection time is unknown (None)")
                    t_start = sta.detections[0].t_start - sta.mean_travel_time
                    t_end = sta.detections[-1].t_end - sta.mean_travel_time
                    if self.t_start is None or self.t_start > t_start:
                        self.set_t_start(t_start)
                    if self.t_end is None or self.t_end < t_end:
                        self.t_end = t_end

        self.data_stations = stations
        self.estimate_source_term()
        self.calculate_confidence_level()
