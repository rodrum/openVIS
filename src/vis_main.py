"""
Main module of the Volcanic Information System. It contains the main algorithm
and should be run to perform an analysis.
"""

import copy
import sys
import pandas as pd

from datetime import datetime, timedelta, timezone
from math import ceil
from typing import cast

sys.path.append(
    "."
)  # TODO: Bad practice but necessary for the modules to work -- see distutil


from src.db_vis import (
    read_volcanoes_from_db,
    get_volcanoes_from_gvp_database,
    load_data,
    get_filtered_detections_from_db,
    get_vratio,
    read_stations_from_db,
    save_eruptions,
    save_ip_list_to_db,
)
from src.eruption import Eruption
from src.logger import logger
from src.infrasoundlib.station import Station
from src.station_volc import StationVolc
from src.infrasoundlib import util
from src.infrasoundlib.detection import Detection
from src.progress import ProgressBar
from src.region import Region
from src.settings import (
    STATIONS_TABLE,
    VOLCANO_TABLE,
    BULL_PATH,
    DATA_PATH,
    START_DATE,
    END_DATE,
    STATIONS,
    MONITORING_AREA,
    REGIONS,
    VOLCANOES,
    ANALYSIS_TIME_INTERVAL,
    DAZIM,
    DELTA_CLOSE,
    IP_THRESHOLD,
    IP_TIME_INTERVAL,
    MAX_AMP,
    MAX_DIST,
    NOTIFICATION_INTERVAL,
    PERSISTENCY_THRESHOLD,
    REF_SPEED,
)


def main() -> None:
    """
    Main function of the VIS. Performs analysis on given period.
    """

    start_d = START_DATE.replace(tzinfo=timezone.utc)
    end_d = END_DATE.replace(tzinfo=timezone.utc)

    volc_tables = get_volcanoes_from_gvp_database(VOLCANO_TABLE)
    stat_tables = pd.read_csv(STATIONS_TABLE)

    lst_eruptions: list[Eruption] = []
    volcanoes, coords = read_volcanoes_from_db(volc_tables)
    stations = read_stations_from_db(stat_tables)
    regions = []

    # Unpack region from volcanoes and coords
    idx = 0
    for volcs_aux, coord in zip(volcanoes, coords):
        if volcs_aux:
            idx += 1
            if coord != []:
                regions.append(Region(coord[0], coord[1], volcs_aux, idx))
            else:
                regions.append(Region(None, None, volcs_aux, idx))

    # Start main analysis
    for region in regions:
        logger.info(
            f"Monitoring {region.name}: {[vol.name for vol in region.volcanoes]} - "
            f"{regions.index(region) + 1}/{len(regions)}"
        )

        # Select nearest stations
        selected_stations = get_stations_for_region(region, stations, start_d, end_d)

        stations_detecting = [s.deepcopy() for s in selected_stations]

        load_data(stat_tables, selected_stations, volc_tables, start_d, end_d)

        cur_erupt = None

        analysis_sd = start_d
        analysis_ed = min(
            analysis_sd + timedelta(seconds=ANALYSIS_TIME_INTERVAL), end_d
        )

        p_bar = ProgressBar(
            ceil((end_d - start_d) / timedelta(seconds=ANALYSIS_TIME_INTERVAL))
        )
        cnt_progress = 1

        start_d_filter = start_d

        while analysis_ed <= end_d:
            # Get detections in bulk of one month or less
            # This is done to reduce the computation time in case of long
            # period analysis.

            if analysis_ed > start_d_filter:
                end_d_filter = start_d_filter + timedelta(days=50)
                if end_d_filter > end_d:
                    end_d_filter = end_d

                selected_stations = get_stations_data(
                    selected_stations, start_d_filter, end_d_filter, region
                )

                for sta_d, sta_s in zip(stations_detecting, selected_stations):
                    sta_d.v_ratio = sta_s.v_ratio
                start_d_filter = end_d_filter

            nb_new_dets = get_new_detections(
                region,
                analysis_ed,
                selected_stations,
                stations_detecting,
            )

            # Ongoing eruption ?
            if cur_erupt:
                # New detections ?
                if nb_new_dets > 0:
                    # YES => Update eruption
                    update_eruption(cur_erupt, stations_detecting)
                else:
                    # NO => Test DELTA_CLOSE ?
                    cur_erupt = test_delta_close(
                        analysis_ed, cur_erupt, lst_eruptions, stations_detecting
                    )
            else:
                # Enough detection to create eruption ?
                if evaluate_ip_threshold(
                    stations_detecting, IP_THRESHOLD, IP_TIME_INTERVAL, analysis_ed
                    ):
                    # OK => New eruption, notification
                    cur_erupt = generate_new_eruption(
                        region, stations_detecting, lst_eruptions, analysis_ed
                        )
                    lst_eruptions.append(cur_erupt)
                else:
                    # Remove isolated detections
                    for sta in stations_detecting:
                        sta.detections = []

            analysis_sd = analysis_ed
            analysis_ed += timedelta(seconds=ANALYSIS_TIME_INTERVAL)
            p_bar.set_and_plot(cnt_progress)
            cnt_progress += 1

        del p_bar

        # In case there is still an ongoing eruption when reaching end date
        if lst_eruptions and lst_eruptions[-1].data_stations[0].nb_dets is None:
            lst_eruptions[-1] = cast(Eruption, copy.deepcopy(cur_erupt))
        
        logger.info('Saving IP values...')
        save_ip_list_to_db(region, stations_detecting)

    if lst_eruptions:
        save_eruptions(lst_eruptions)


def remove_isolated_detections(stations_detecting: list[StationVolc]) -> None:
    """
    Remove detections from stations.

    Parameters
    ----------
    stations_detecting : :obj:`list` of :class:`.Station`
        List of detecting Station.

    """
    for sta in stations_detecting:
        sta.detections = []


def generate_new_eruption(
    region: Region,
    stations_detecting: list[StationVolc],
    lst_eruptions: list[Eruption],
    analysis_ed: datetime,
    ) -> Eruption:
    """
    Generate an eruption code if another eruption already exist with the same id.

    Parameters
    ----------
    region : :class:`.Region`
        Monitored region in eruption.
    stations_detecting : :obj:`list` of :class:`.StationVolc`
        List of Station detecting.
    lst_eruptions : :obj:`list` of :class:`.Eruption`
        List of past eruptions.
    analysis_ed : datetime
        End date of analysis.

    Returns
    -------
    :class:`.Eruption`
        New Eruption.

    """
    cur_eruption = Eruption(
        region,
        stations_detecting,
        thresh_ip=IP_THRESHOLD,
        ip_interval=IP_TIME_INTERVAL
    )
    er_same_id = [
        e for e in lst_eruptions if e.er_code.startswith(cur_eruption.er_code)
    ]
    if er_same_id:
        cur_eruption.er_code = cur_eruption.er_code + "-" + str(len(er_same_id))

    cur_eruption.last_notification = datetime.now(timezone.utc)
    cur_eruption.detection_date = analysis_ed
    return cur_eruption


def evaluate_ip_threshold(
    stations: list[StationVolc],
    thresh_ip_val: float,
    ip_interval: int,
    ed: datetime,
    nb_ip: int = 1,
    ) -> bool:
    """
    Evaluate if the IP values exceed the threshold `thresh_ip_val`.

    Evaluate the IP for each given `stations` on the interval `ip_interval`.
    The `nb_ip` parameter can also be set to specify the minimum number of value
    that needs to exceed the threshold.

    Parameters
    ----------
    stations : :obj:`list` of :class:`.StationVolc`
        List of Station.
    thresh_ip_val : float
        Threshold IP value.
    ip_interval : int
        IP interval time (in number of seconds).
    ed : datetime
        Date of analysis.
    nb_ip : int, optional
        Minimum number of IP values that needs to exceed the threshold, defaults to 1.

    Returns
    -------
    bool
        Threshold exceeded for at least `nb_ip` times.

    """
    tmp_ip = 0
    for sta in stations:
        if sta.detections:
            tmp_ip += len(
                [
                    ip
                    for ip in sta.list_ip
                    if ip.value >= thresh_ip_val
                    and ip.dt >= ed - timedelta(seconds=ip_interval)
                ]
            )
    return bool(tmp_ip >= nb_ip)


def update_eruption(cur_erupt: Eruption, stations_detecting: list[StationVolc]) -> None:
    """
    Update Eruption and generate notification if needed.

    Also take care of removing unwanted detection in real time mode.

    Parameters
    ----------
    cur_erupt : :class:`.Eruption`
        Current Eruption
    stations_detecting : :obj:`list` of :class:`.StationVolc`
        List of Station detecting.
    """

    cur_erupt.update_eruption(
        stations_detecting, thresh_ip=IP_THRESHOLD, ip_inter=IP_TIME_INTERVAL
    )

    if cur_erupt.t_start is None:
        raise TypeError("Eruption has no start time (None)")

    nxt_not_time = cur_erupt.t_start + cur_erupt.revision * timedelta(
        seconds=NOTIFICATION_INTERVAL
    )

    if cur_erupt.t_end is None:
        raise TypeError("Eruption has no end time (None)")

    if nxt_not_time < cur_erupt.t_end:
        cur_erupt.revision += 1
        cur_erupt.last_notification = datetime.now(timezone.utc)


def test_delta_close(
    analysis_ed: datetime,
    cur_erupt: Eruption,
    lst_eruptions: list[Eruption],
    stations_detecting: list[StationVolc],
) -> Eruption | None:
    """
    Test delta close and end eruption if needed.

    Parameters
    ----------
    analysis_ed : datetime
        End date of analysis.
    cur_erupt : :class:`.Eruption`
        Current Eruption.
    lst_eruptions : :obj:`list` of :class:`.Eruption`
        List of Eruption.
    stations_detecting : :obj:`list` of :class:`.StationVolc`
        List of Station detecting.

    Returns
    -------
    :class:`.Eruption` or None
        Eruption or None if the eruption ended.
    """

    max_avg_tt = max((s.mean_travel_time for s in stations_detecting))

    if cur_erupt.t_end is None:
        raise TypeError("Eruption has no end time (None)")
    if analysis_ed - cur_erupt.t_end - max_avg_tt > timedelta(seconds=DELTA_CLOSE):
        cur_erupt.update_eruption(
            stations_detecting, thresh_ip=IP_THRESHOLD, ip_inter=IP_TIME_INTERVAL
        )
        cur_erupt.status = "ENDED"
        cur_erupt.revision += 1
        cur_erupt.last_notification = datetime.now(timezone.utc)

        lst_eruptions[-1] = copy.deepcopy(cur_erupt)
        for sta in stations_detecting:
            sta.clear()
        return None

    return cur_erupt


def get_new_detections(
    region: Region,
    end_t: datetime,
    selected_stations: list[StationVolc],
    stations_detecting: list[StationVolc],
    ) -> int:
    """
    For each station, add the filtered detections to their detections,
    calculate the IP value and count the number of new detections on
    the given period.

    Parameters
    ----------
    end_t : datetime
        End date.
    selected_stations : :obj:`list` of :class:`.StationVolc`
        List of stations (with all detections).
    stations_detecting : :obj:`list` of :class:`.StationVolc`
        List of stations (with filtered detections).

    Returns
    -------
    int
        Number of new detections

    """
    nb_new_dets = 0
    start_t = end_t - timedelta(seconds=IP_TIME_INTERVAL)
    for idx_sta, sta in enumerate(selected_stations):
        # Select and filter detections by azimuth and frequency
        det_filtered = filter_detections(
            sta, start_t, end_t, sta.distaz["ddeg"]
        )
        if det_filtered:
            # Do not include the same detection twice 
            if stations_detecting[idx_sta].detections:
                det_filtered = [
                    d
                    for d in det_filtered
                    if (
                        d.t_start
                        > stations_detecting[idx_sta].detections[-1].t_start
                    )
                ]
            stations_detecting[idx_sta].detections.extend(det_filtered)
            nb_new_dets += len(det_filtered)
        nb_new_dets -= stations_detecting[idx_sta].calculate_ip(
            end_t, IP_TIME_INTERVAL, max_amp=MAX_AMP,
            persistency_threshold=PERSISTENCY_THRESHOLD
        )

    return nb_new_dets


def filter_detections(
    station: StationVolc,
    t_min: datetime,
    t_max: datetime,
    ddeg: float,
    ) -> list[Detection]:
    """
    Filter `detections` with a start time between `t_min` and `t_max`,
    and a max frequency lower than the theoretical max frequency
    (calculated with `ddeg`).

    Parameters
    ----------
    station : :class:`.Station`
        Station.
    t_min : datetime
        Minimum time for detection to start.
    t_max : datetime
        Maximum time for detection to end.
    ddeg : float
        Great circle arc distance in degrees.

    Returns
    -------
    :obj:`list` of :class:`.Detection`
        Sorted list of Detection (by start time).

    """
    # NOTE-1: We made a rough filter on azimuth when querying from the db
    # in order to limit the number of detections
    # (as it slows down the program significantly).
    # NOTE-2: Does not take into account the unpractical case where a station
    # is located inside a region.
    # NOTE-3: When referring to 'baz' it means azimuth from station to volcano
    # min_baz = station.baz_min  # min azimuth from station to volcano(es)
    # max_baz = station.baz_max  # max ''
    # IMPORTANT: Assumption below is from -------------
    # Brachet et al. 2010, Section 3.2.3 (if ddeg < 60)
    max_frequency = 4 - 0.055 * ddeg
    # -------------------------------------------------

    filtered_dets = [
        d
        for d in station.detections
        if (t_min <= d.t_start < t_max and d.f_mean <= max_frequency)
    ]
    #if min_baz < max_baz:
    #    filtered_dets = [d for d in filtered_dets if min_baz < d.azi < max_baz]
    #else: # wrap-around case (e.g., min=355, max=5)
    #    filtered_dets = [
    #        d for d in filtered_dets if (max_baz > d.azi or min_baz < d.azi)
    #    ]
    return sorted(filtered_dets, key=Detection.get_tstart)


def get_stations_data(
    selected_stations: list[StationVolc],
    start_d: datetime,
    end_d: datetime,
    region: Region,
) -> list[StationVolc]:
    """
    Get the detections and vratios for the stations selected on the analysis
    period and filter stations that don't have detections.

    Parameters
    ----------
    region : :class:`.Region`
        Region monitored by the stations.
    start_d : datetime
        Start date.
    end_d : datetime
        End date.
    selected_stations : :obj:`list` of :class:`.StationVolc`
        List of stations selected for the analysis.

    Returns
    -------
    :obj:`list` of :class:`.Station`
        List of updated Station

    """
    for sta in selected_stations:
        # ===================
        # Retrieve detections
        # ===================
        sta.detections = get_filtered_detections_from_db(
            sta,
            start_d,
            end_d,
            region,
        )
        # ====================
        # Retrieve Veff-ratios
        # ====================
        if get_vratio(sta, start_d, end_d) == 1:
            sta.v_ratio = None
        else:
            sta.v_ratio = get_vratio(sta, start_d, end_d)
    return selected_stations


def get_current_eruption_data(
    cur_erupt: Eruption, stations_detecting: list[StationVolc]
) -> None:
    """
    Retrieve the station data from the database that is related with the
    current eruption.

    Parameters
    ----------
    cur_erupt : :class:`.Eruption`
        Current Eruption.
    stations_detecting : :obj:`list` of :class:`.Station`
        List of detecting stations.
    """
    for sta in stations_detecting:
        sta_from_db = [s for s in cur_erupt.data_stations if s.name == sta.name]
        if len(sta_from_db) == 1:
            sta.estimated_amplitude = sta_from_db[0].estimated_amplitude
            sta.max_amplitude = sta_from_db[0].max_amplitude
            sta.detecting = sta_from_db[0].detecting
            sta.nb_dets_previous_analysis = sta_from_db[0].nb_dets_previous_analysis


def get_stations_for_region(
    region: Region,
    stations: list[Station],
    start_d: datetime,
    end_d: datetime,
) -> list[StationVolc]:
    """
    Select closest stations to region and load their detections and vratios.

    Parameters
    ----------
    region : :class:`.Region`
        Monitored Region.
    stations : :obj:`list` of :class:`.Station`
        List of Station given for the analysis.
    start_d : datetime
        Start date.
    end_d : datetime
        End date.

    Returns
    -------
    :obj:`list` of :class:`.Station`
        List of Station for the analysis.

    """
    selected_stations = select_closest_stations(region, stations)
    if len(selected_stations) == 0:
        logger.warning(
            "No station found for monitored region - "
            "Skipping to next volcano"
        )
        return []
    logger.info(f"Stations selected: {' - '.join([s.name for s in selected_stations])}")

    for sta in selected_stations:
        sta.compute_baz_filter(region.volcanoes, DAZIM)

    if len(selected_stations) == 0:
        logger.warning("No detections found for selected - Skipping to next region")
        return []

    return selected_stations


def select_closest_stations(
    region: Region,
    stations: list[Station]
) -> list[StationVolc]:
    """
    Select `stations` that are located within a maximum distance of
    the `volcano` for the analysis.

    Parameters
    ----------
    region : :class:`.Region`
        Monitored Region.
    stations : :obj:`list` of :class:`.Station`
        List of Station given for the analysis.

    Returns
    -------
    :obj:`list` of :class:`.StationVolc`
        List of Station.

    """
    selected_stations = []
    for sta in stations:
        distaz = util.dist_az(sta.lat, sta.lon, region.lat, region.lon)
        if distaz["dkm"] <= MAX_DIST:
            sta_volc = StationVolc.from_station(sta)
            sta_volc.distaz = distaz
            sta_volc.mean_travel_time = timedelta(
                seconds=(sta_volc.distaz["dkm"] / REF_SPEED)
            )
            selected_stations.append(sta_volc)
    return selected_stations


if __name__ == '__main__':
    main()
