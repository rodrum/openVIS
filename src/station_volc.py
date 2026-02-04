"""
This module defines the StationVolc class
"""

import copy
from datetime import datetime, timedelta
from typing import Any

import numpy as np

from src.infrasound_parameter import IP
from src.infrasoundlib import attenuation, util
from src.infrasoundlib.detection import Detection
from src.infrasoundlib.station import Station
from src.logger import logger
from src.volcano import Volcano


class StationVolc(Station):
    """
    A class to represent a :class:`.Station` monitoring a :class:`.Volcano`. Children of
    :class:`Station`.

    Attributes
    ----------
    distaz : :obj:`dict`
        A dictionnary containing information about the volcano monitored. The azimuth, back-azimuth
        and distance in km and degrees.
    baz_min : float
        Minimum back azimuth used to filter detections.
    baz_max : float
        Maximum back azimuth used to filter detections.
    mean_travel_time : datetime.timedelta
        Mean travel time of the infrasound wave from the volcano to the station.
    list_ip : :obj:`list` of :class:`.IP`
        List of IP object.
    max_amplitude : float
        Max amplitude of the detections.
    estimated_amplitude : float
        Estimated amplitude of the eruption at the station.
    detecting : bool
        Indicate if the station is detecting the eruption (True if IP exceed the threshold one time
        during the analysis).
    nb_dets_previous_analysis : int
        Number of detections from the previous analysis. Only used when there is an ongoing eruption
        saved in the database at the start analysis time of VIS.

    """

    def __init__(self, name: str, lat: float, lon: float):
        """
        Constructor of :class:`.StationVolc`.

        Parameters
        ----------
        name : str
            Name of the station.
        lat : float
            Latitude of the station.
        lon : float
            Longitude of the station.

        """

        super().__init__(name, lat, lon)
        self.distaz: dict[str, float] = {}
        self.mean_travel_time = timedelta(seconds=0)
        self.list_ip: list[IP] = []
        self.max_amplitude = 0.0
        self.estimated_amplitude = 0.0
        self.detecting = False
        self.nb_dets_previous_analysis = 0
        self.baz_max: int = 360
        self.baz_min: int = 0

    def __deepcopy__(self, memo: dict[int, Any]) -> "StationVolc":
        """
        Only keep relevant info to write in summary. Override deepcopy.
        """

        result = StationVolc(self.name, self.lat, self.lon)
        result.detecting = self.detecting
        result.nb_dets = len(self.detections)
        result.nb_dets_previous_analysis = self.nb_dets_previous_analysis
        result.max_amplitude = self.max_amplitude
        result.estimated_amplitude = self.estimated_amplitude
        result.mean_travel_time = self.mean_travel_time
        return result

    def compute_baz_filter(self, volcanoes: list[Volcano], dazim: int) -> None:
        """
        For a given set of volcanoes, compute the min and max azimuth for the 
        station to see every volcano.

        Parameters
        ----------
        volcanoes : :obj:`list` of :class:`.Volcano`
            List of volcanoes.
        dazim : int
            Delta azimuth to account for wind deflection.
        end_d: datetime
            End date.
        """
        baz_min = self.distaz["baz"]
        baz_max = self.distaz["baz"]  # NOTE: why are both the same?
        for vol in volcanoes:
            # NOTE: azimuth from station to volcano
            new_baz = util.dist_az(self.lat, self.lon, vol.lat, vol.lon)["baz"]
            if not np.mod(new_baz - baz_min, 360) <= np.mod(baz_max - baz_min, 360):
                diff_min = np.mod(baz_min - new_baz, 360)
                diff_max = np.mod(new_baz - baz_max, 360)
                if diff_min < diff_max:
                    baz_min = new_baz
                elif diff_max < diff_min:
                    baz_max = new_baz
                else:
                    baz_min = 0
                    baz_max = 360
        if baz_min < dazim and baz_max > 360 - dazim:
            self.baz_min = 0
            self.baz_max = 360
        else:
            self.baz_min = np.mod(baz_min - dazim, 360)
            self.baz_max = np.mod(baz_max + dazim, 360)

    def deepcopy(self) -> "StationVolc":
        """
        Deepcopy is overrided to only keep relevant info to write in summary.
        Before this patch, the RAM usage of VIS was enormous,
        especially on long period of analysis.

        This function redefines the deepcopy (basic usage) and should be used
        when wanting to perform a real deepcopy. Note that the detections are
        not deepcopied (we don't need to in its current usage).

        This was implemented this way because the __deepcopy__ function of
        StationVolc is the one called when __deepcopy__ is called on an object
        containing StationVolc objects (a list for example).


        """
        memo = {}
        cls = self.__class__  # Extract the class of the object
        result = cls.__new__(
            cls
        )  # Create a new instance of the object based on extracted class
        memo[id(self)] = result
        for key, value in self.__dict__.items():
            if key == "detections":
                pass
            else:
                # Copy over attributes by copying directly or in case of complex
                # objects like lists for exaample calling the `__deepcopy()__`
                # method defined by them. Thus recursively copying the whole
                # tree of objects.
                setattr(result, key, copy.deepcopy(value, memo))
        result.detections = []
        return result

    @classmethod
    def from_station(cls, station: Station) -> "StationVolc":
        """
        Create a :class:`.StationVolc` object from a :class:`.Station` object.

        Parameters
        ----------
        station : :class:`.Station`
            Station.

        Returns
        -------
        :class:`.StationVolc`
            StationVolc.

        """

        sta = cls(station.name, station.lat, station.lon)
        sta.detections = station.detections
        sta.v_ratio = station.v_ratio
        return sta

    def calculate_ip(self,
                     ed: datetime,
                     ip_interval: int,
                     max_amp: float = 500.0,
                     persistency_threshold: int = 20) -> int:
        """
        Calculate the IP at the date `ed` using detections from the last `ip_interval` number
        of seconds and filtering every detections with an amplitude higher than `max_amp`.

        The IP calculation corresponds to two filters: 
        1) The first filter removes every detection which has an amplitude higher than `max_amp`. 
        2) The second filter is on the persistency of the detections: the ratio of the time on 
        which the station receives detections and the `ip_interval` number of seconds must be 
        higher than the IP threshold.

        If detections have a too high amplitude, they are completely removed from the detection
        list of the station.

        Parameters
        ----------
        ed : datetime
            Date at which the IP calculation is done.
        ip_interval : int
            Number of seconds if the ip time interval.
        max_amp : float, optional
            Maximum amplitude for a detection to be considered trustworthy.

        Returns
        -------
        int
            Number of detections removed from the station.
        """
        det_removed = 0
        sd = ed - timedelta(seconds=ip_interval)
        if not ip_interval or ip_interval < 0:
            logger.critical("IP_interval parameter must be > 0")
            raise ValueError
        detections: list[Detection] = [
            d for d in self.detections if sd <= d.t_start < ed  # type: ignore
        ]
        vratio = attenuation.interpolate_vratio(self.v_ratio, sd)
        if detections:
            for det in detections:
                if det.f_mean is None:
                    raise TypeError("Mean freq of detection is unknown (None)")
                det.attenuation = attenuation.calculate_att_coeff(
                    det.f_mean, vratio, self.distaz["dkm"]
                )
            dets_to_remove = [
                d for d in detections if d.amp / d.attenuation > max_amp  # type: ignore
            ]
            det_removed += len(dets_to_remove)
            self.detections: list[Detection] = [
                d for d in self.detections if d not in dets_to_remove
            ]
            detections = [d for d in detections if d not in dets_to_remove]
            if detections:
                intervals: list[datetime] = []
                #time_det = timedelta()
                for det in detections:
                    if det.t_end is None or det.t_start is None:
                        raise TypeError("Time of detection is unknown (None)")
                    t_end = det.t_end if det.t_end < ed else ed
                    intervals = util.interval_processing(
                        intervals, [t_end, det.t_start]
                    )
                time_det = np.sum(
                    [
                        intervals[i + 1] - intervals[i]
                        for i in range(0, len(intervals)-1, 2)
                    ]
                )
                # convert to seconds
                persistency = time_det.seconds / ip_interval * 100 # pylint: disable=E1101
                avg_source_amp = np.mean(
                    [d.amp / d.attenuation for d in detections]  # type: ignore
                )
                avg_amp = np.mean([d.amp for d in detections])
                mean_freq = np.mean([d.f_mean for d in detections])

                if persistency > persistency_threshold:
                    ip_val = persistency * avg_source_amp
                    self.list_ip.append(
                        IP(
                            dt=ed,
                            val=ip_val,
                            nb_det=len(detections),
                            avg_amp=avg_amp,
                            avg_source_amp=avg_source_amp,
                            persistency=persistency,
                            fmean=mean_freq,
                            vratio=vratio,
                        )
                    )
                else:
                    self.list_ip.append(
                        IP(
                            dt=ed,
                            val=0,
                            nb_det=len(detections),
                            avg_amp=avg_amp,
                            avg_source_amp=avg_source_amp,
                            persistency=persistency,
                            fmean=mean_freq,
                            vratio=vratio,
                        )
                    )
            else:
                self.list_ip.append(IP(ed, vratio=vratio))
        else:
            self.list_ip.append(IP(ed, vratio=vratio))

        return det_removed

    def clear(self) -> None:
        """
        Clear the station info.
        """
        self.detections = []
        self.detecting = False
        self.max_amplitude = 0
        self.estimated_amplitude = 0
        self.nb_dets_previous_analysis = 0
