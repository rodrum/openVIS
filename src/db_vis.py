"""
This module deals with everything related to the database.
"""

import os
from os.path import join, isfile, isdir
from fnmatch import fnmatch
import sys
from datetime import datetime, timedelta, timezone
import scipy.io as sio
import pickle
import netCDF4
from cftime import num2date
from itertools import product

import numpy as np
import pandas as pd
import xarray as xr

from tqdm import tqdm

sys.path.append(".")

from src.logger import logger
from src.infrasoundlib.station import Station
from src.eruption import Eruption
from src.infrasoundlib.detection import Detection
from src.infrasoundlib.vratio import Vratio
from src.region import Region
from src.settings import (
    ANALYSIS_TIME_INTERVAL,
    BULL_PATH,
    BAZDEV_PATH,
    IP_TIME_INTERVAL,
    MIN_MEAN_FREQ,
    MAX_MEAN_FREQ,
    MONITORING_AREA,
    REGIONS,
    VOLCANOES,
    DAZIM,
    FORCE_DAZIM,
    NUM_BAZ_STD,
    RESULTS_PATH,
    VEFF_PATH,
    VEFF_FORMAT,
    DATA_PATH,
    STATIONS,
    START_DATE,
    END_DATE
)
from src.station_volc import StationVolc
from src.volcano import Volcano
from src.models import create_output_dataframes
create_output_dataframes()


def get_stat_name(stat_name):
    """
    Temporary patch to find interpolations following database names of stations
    Used for back-azimuth deviation based tolerances.
    NOTE: fix this
    """
    stat_dic = {
        'IS01': 'I01AR',
        'IS02': 'I02AR',
        'IS03': 'I03AU',
        'IS04': 'I04AU',
        'IS05': 'I05AU',
        'IS06': 'I06AU',
        'IS07': 'I07AU',
        'IS08': 'I08BO',
        'IS09': 'I09BR',
        'IS10': 'I10CA',
        'IS11': 'I11CV',
        'IS12': 'I12CR',
        'IS13': 'I13CL',
        'IS14': 'I14CL',
        'IS15': 'I15CN',
        'IS16': 'I16CN',
        'IS17': 'I17CI',
        'IS18': 'I18DK',
        'IS19': 'I19DJ',
        'IS20': 'I20EQ',
        'IS21': 'I21FR',
        'IS22': 'I22FR',
        'IS23': 'I23FR',
        'IS24': 'I24FR',
        'IS25': 'I25FR',
        'IS26': 'I26DE',
        'IS27': 'I27DE',
        'IS29': 'I29IR',
        'IS30': 'I30JP',
        'IS31': 'I31KZ',
        'IS32': 'I32KE',
        'IS33': 'I33MG',
        'IS34': 'I34MN',
        'IS35': 'I35NA',
        'IS36': 'I36NZ',
        'IS37': 'I37NO',
        'IS38': 'I38PA',
        'IS39': 'I39PW',
        'IS40': 'I40PN',
        'IS41': 'I41PY',
        'IS42': 'I42PT',
        'IS43': 'I43RU',
        'IS44': 'I44RU',
        'IS45': 'I45RU',
        'IS46': 'I46RU',
        'IS47': 'I47ZA',
        'IS48': 'I48TN',
        'IS49': 'I49GB',
        'IS50': 'I50GB',
        'IS51': 'I51GB',
        'IS52': 'I52GB',
        'IS53': 'I53US',
        'IS54': 'I54US',
        'IS55': 'I55US',
        'IS56': 'I56US',
        'IS57': 'I57US',
        'IS58': 'I58US',
        'IS59': 'I59US',
        'IS60': 'I60US',
        'DBN': 'DBNI',
        'JAM': 'JAMT',
        'KIR': 'KIRU',
        'LYC': 'LYCK',
        'SOD': 'SODA',
    }
    if stat_name in stat_dic:
        return stat_dic[stat_name]
    else:
        return stat_name


def get_volcanoes_region(volcs, top_left_corner, bot_right_corner):
    """
    Get the volcanoes located inside the region. Compute the barycenter of the region.

    Parameters
    ----------
    session : Session
        SQLAlchemy session object linked to local database
    top_left_corner : :obj:`tuple` of float
        Coordinates of the top left corner of the region.
    bot_right_corner : :obj:`tuple` of float
        Coordinates of the bottom right corner of the region.

    Returns
    -------
    :obj:`list` of :class:`.Volcano`
        List of volcanoes in region.
    """
    lst_volcanoes = []
    lat_max = top_left_corner[0]
    lon_min = top_left_corner[1]
    lat_min = bot_right_corner[0]
    lon_max = bot_right_corner[1]
    volcs_aux = None
    if lat_max < lat_min:
        volcs_aux_1 = volcs[volcs['Latitude'] >= lat_max]
        volcs_aux_2 = volcs[volcs['Latitude'] <= lat_min]
        volcs_aux = pd.merge(volcs_aux_1, volcs_aux_2, 'outer')
    else:
        volcs_aux_1 = volcs[volcs['Latitude'] <= lat_max]
        volcs_aux_2 = volcs[volcs['Latitude'] >= lat_min]
        volcs_aux = pd.merge(volcs_aux_1, volcs_aux_2, 'inner')

    if lon_max < lon_min:
        volcs_aux_1 = volcs_aux[volcs_aux['Longitude'] >= lon_max]
        volcs_aux_2 = volcs_aux[volcs_aux['Longitude'] <= lon_min]
        volcs_aux = pd.merge(volcs_aux_1, volcs_aux_2, 'outer')
    else:
        volcs_aux_1 = volcs_aux[volcs_aux['Longitude'] <= lon_max]
        volcs_aux_2 = volcs_aux[volcs_aux['Longitude'] >= lon_min]
        volcs_aux = pd.merge(volcs_aux_1, volcs_aux_2, 'inner')

    for row in range(len(volcs_aux)):
        lst_volcanoes.append(
            [
                volcs_aux['Volcano Number'].iloc[row],
                volcs_aux['Volcano Name'].iloc[row],
                volcs_aux['Latitude'].iloc[row],
                volcs_aux['Longitude'].iloc[row],
                volcs_aux['Elevation (m)'].iloc[row],
            ]
        )
    return lst_volcanoes


def input_is_tuple_coordinate(region):
    if (
        len(region) == 2
        and len(region[0]) == 2
        and len(region[1]) == 2
        and isinstance(region[0][0], (float, int))
        and isinstance(region[0][1], (float, int))
        and isinstance(region[1][0], (float, int))
        and isinstance(region[1][1], (float, int))
    ):
        return True
    return False


def get_volcano_name(volcs, volc: str):
    """
    Get a volcano from its name in the database.

    Parameters
    ----------
    volcs : Pandas DataFrame with volcano list
    volc : str
        Name of the volcano.

    Returns
    -------
        Volcano row from the dataframe.
    """

    volc_db = volcs[volcs['Volcano Name'] == volc]  # NOTE: case sensitive
    if len(volc_db) == 1:
        return volc_db
    else:
        logger.warning(f"No volcano found in database for {volc} - Skipping")
        return None


def volcano_from_dbo(volc_db) -> Volcano:
    """
    Create :class:`.Volcano` from volcano row in table.

    Parameters
    ----------
    volc_db : selected 'row' from dataframe

    Returns
    -------
    :class:`.Volcano`
        Volcano object used for calculations.

    """

    return Volcano(
        volc_db['Volcano Number'].item(),
        volc_db['Volcano Name'].item(),
        volc_db['Latitude'].item(),
        volc_db['Longitude'].item(),
        volc_db['Elevation (m)'].item()
    )


def get_volcano_code(volcs, volc: int):
    """
    Get a volcano from its name in the database.

    Parameters
    ----------
    volcs : pandas dataframe
    volc : int
        Code of the volcano.

    Returns
    -------
    pandas dataframe
        Volcano from the database.
    """

    volc_db = volcs[volcs['Volcano Number'] == volc]
    if len(volc_db) == 1:
        return volc_db
    else:
        logger.warning(f"No volcano found in database for {volc} - Skipping")
        return None


def get_volcanoes_from_gvp_database(volcano_path: str) -> pd.DataFrame:
    """
    Get a list of interesting volcanoes from the very complete GVP database.

    Only volcanoes that have a recorded eruption since 2003 are selected in the
    VIS database.

    Returns
    -------
    :class:`pandas.DataFrame`
        Filtered dataframe containing the volcanoes.

    Notes
    -----
    The GVP database can be downloaded at:
    https://volcano.si.edu/volcanolist_holocene.cfm.
    """

    volcanoes_all_data = pd.read_csv(
        volcano_path,
        sep=";",
        encoding="latin-1",
        dtype={'Volcano Number': int,
               'Volcano Name': str,
               'Country': str,
               'Last Known Eruption': str,
               'Latitude': float,
               'Longitude': float,
               'Elevation (m)': float,
               },
        decimal=','
    )

    volcanoes = volcanoes_all_data[
       [
            "Volcano Number",
            "Volcano Name",
            "Country",
            "Last Known Eruption",
            "Latitude",
            "Longitude",
            "Elevation (m)",
        ]
    ].copy()
    lke = "Last Known Eruption"
    # Filter all volcanoes with lke unknown and BCE (BC equivalent)
    volcanoes.drop(
        volcanoes[
            (volcanoes[lke] == "Unknown") | (volcanoes[lke].str.contains("BCE"))
        ].index,
        inplace=True,
    )
    # Filter all volcanoes with lke before 2003
    # (year where we have most of the IMS data)
    volcanoes.drop(
        volcanoes[
            volcanoes[lke].str.extract(r"([\d]+)", expand=False).astype(int) < 2003  
        ].index,
        inplace=True,
    )
    # Create escape character for sql request
    volcanoes.replace("'", "''", regex=True, inplace=True)
    # replace coma by points in number
    volcanoes.replace(",", ".", regex=True, inplace=True)

    volcanoes = pd.DataFrame(
        {
            "Volcano Number": volcanoes['Volcano Number'].to_numpy(),
            "Volcano Name": volcanoes['Volcano Name'].to_numpy(),
            "Country": volcanoes['Country'].to_numpy(),
            "Last Known Eruption": volcanoes['Last Known Eruption'].to_numpy(),
            "Latitude": volcanoes['Latitude'].to_numpy(),
            "Longitude": volcanoes['Longitude'].to_numpy(),
            "Elevation (m)": volcanoes['Elevation (m)'].to_numpy(),
        },
    )

    return volcanoes


def read_volcanoes_from_db(volcs_tab):
    """
    Get volcanoes info from the database. The list of :class:`.Volcano` is
    parsed from the config file.

    Parameters
    ----------
    volcs_tab : pandas dataframe
        Contains table volcanoes defined by csv file (<volcanoes.csv>).

    """

    lst_volcanoes = []
    lst_region_coords = []

    if input_is_tuple_coordinate(MONITORING_AREA):
        top_left_corner: list[float, float] = MONITORING_AREA[0]
        bot_right_corner: list[float, float] = MONITORING_AREA[1]
        lst_region_coords = [top_left_corner, bot_right_corner]
        lst_volcanoes = get_volcanoes_region(
                volcs_tab, top_left_corner, bot_right_corner
            )
    else:
        for region in REGIONS:
            lst_vol_region = []
            lst_region_coords_aux = []
            logger.info(f"Coniguring volcanic region: {region}")
            for vol in region:
                volc_db = None
                if isinstance(vol, int):
                    volc_db = get_volcano_code(volcs_tab, vol)
                elif isinstance(vol, str):
                    volc_db = get_volcano_name(volcs_tab, vol)
                if volc_db is not None:
                    lst_vol_region.append(volcano_from_dbo(volc_db))
                    lst_region_coords_aux.append([])
            lst_region_coords.append(lst_region_coords_aux)
            lst_volcanoes.append(lst_vol_region)

        for volc in VOLCANOES:
            if isinstance(volc, int):
                logger.info(f"Volcano by code: {volc}")
                volc_db = get_volcano_code(volcs_tab, volc)
                if volc_db is not None:
                    lst_volcanoes.append([volcano_from_dbo(volc_db)])
                    lst_region_coords.append([])
            elif isinstance(volc, str):
                logger.info(f"Volcano by name: {volc}")
                volc_db = get_volcano_name(volcs_tab, volc)
                if volc_db is not None:
                    lst_volcanoes.append([volcano_from_dbo(volc_db)])
                    lst_region_coords.append([])

    return lst_volcanoes, lst_region_coords


def read_stations_from_db(stats_tab):
    """
    Get stations infos from the database. The list of :class:`.Station` is
    parsed from the config file.

    Parameters
    ----------
    stats_tab : Pandas DataFrame
        Table with stations in network.
    """

    lst_stations = []
    if len(STATIONS) == 0:
        logger.critical("No station selected for the analysis")
        sys.exit()

    for sta_nam in STATIONS:
        sta_db = stats_tab[stats_tab['Station Name'] == sta_nam]
        lst_stations.append(station_from_dbo(sta_db))

    if len(lst_stations) == 0:
        logger.critical("No station found for the analysis")
        sys.exit()

    return lst_stations


def station_from_dbo(sta_db) -> Station:
    """
    Create :class:`.Station` from :class:`.StationDB`.

    Parameters
    ----------
    sta_db : dataframe 'stats'

    Returns
    -------
    :class:`.Station`
        Station object used for calculations.

    """

    return Station(sta_db['Station Name'].values[0],
                   sta_db['Latitude'].values[0],
                   sta_db['Longitude'].values[0])


def get_vratio(sta: StationVolc, t_min: datetime, t_max: datetime,
               veff_path=join(DATA_PATH, 'veff_ratios.pkl')):
    """
    Assuming vratio is a table with one record per time of day (usually every
    6 hours) and angle (0 to 359 in steps of 1 degree).
    Hence dimension Datetimes x 360, whatever Datatimes is.
    In 365 days, four times a day, this would be [365x4 x 360].
    The headers are: 'dt', 'name', 'azimuth' and 'value' (see 'load_data')
    """
    veff_sta = pd.read_pickle(veff_path)
    if type(veff_sta) == int:  # a terrible way of doing this, FIXME
        return 1
    vratios_out = []
    vratios = veff_sta[veff_sta['dt'] <= t_max]
    vratios = vratios[vratios['dt'] >= t_min]
    vratios = vratios[vratios['azimuth'] == int(round(sta.distaz["baz"]))%360]
    for row in range(len(vratios)):
        vratios_out.append(
            Vratio(
                vratios['dt'].iloc[row].replace(tzinfo=timezone.utc),
                vratios['azimuth'].iloc[row],
                vratios['value'].iloc[row]
                )
        )
    return vratios_out


def save_eruptions(eruptions: list[Eruption]) -> None:
    """
    Save a list of eruption in the database

    Parameters
    ----------
    eruptions : :obj:`list` of :class:`.Eruption`

    """

    results_filename = join(RESULTS_PATH, 'eruption_results.pkl')
    assoc_sta_er_filename = join(RESULTS_PATH, 'assoc_sta_er.pkl')
    assoc_volc_er_filename = join(RESULTS_PATH, 'assoc_volc_er.pkl')
    eruption_results = pd.read_pickle(results_filename)
    assoc_sta_er = pd.read_pickle(assoc_sta_er_filename)
    assoc_volc_er = pd.read_pickle(assoc_volc_er_filename)

    for eruption in eruptions:
        # Get data volcanoes in eruption, first one suffices

        er_db = eruption_results[eruption_results['Eruption Code'] == eruption.er_code].index

        if len(er_db) == 0:
            for volc in eruption.region.volcanoes:
                this_volc_er = pd.DataFrame(
                    [
                        [
                            eruption.er_code,
                            volc.code,
                            volc.name,
                        ]
                    ],
                    columns=assoc_volc_er.columns
                )

                assoc_volc_er = pd.concat(
                    [
                        this_volc_er,
                        assoc_volc_er
                    ],
                    ignore_index=True
                )

            for sta in eruption.data_stations:
                this_sta_er = pd.DataFrame(
                    [
                        [
                            sta.name,
                            eruption.er_code,
                            sta.nb_dets + sta.nb_dets_previous_analysis,
                            round(sta.max_amplitude, 5),
                            round(sta.estimated_amplitude),
                            sta.detecting,
                        ]
                    ],
                    columns=assoc_sta_er.columns
                )

                assoc_sta_er = pd.concat(
                    [
                        this_sta_er,
                        assoc_sta_er
                    ],
                    ignore_index=True
                )

            this_eruption = pd.DataFrame(
                [
                    [
                        eruption.er_code,
                        eruption.t_start,
                        eruption.t_end,
                        eruption.detection_date,
                        eruption.confidence_level,
                        eruption.last_notification,
                        eruption.revision,
                        str.lower(eruption.status),
                        round(eruption.amp_source),
                    ]
                ],
                columns=eruption_results.columns
            )

            eruption_results = pd.concat(
                [
                    this_eruption,
                    eruption_results
                ],
                ignore_index=True
            )
        else:
            # Update every eruption
            if (
                eruption.t_start is None
                or eruption.t_end is None
                or eruption.last_notification is None
            ):
                raise TypeError("Eruption from DB has no time info (None)")

            for ind in er_db:  # NOTE: should be just a list with one element
                eruption_results.at[ind, 'Start Date (UTC)'] = \
                    eruption.t_start.replace(tzinfo=timezone.utc)
                eruption_results.at[ind, 'End Date (UTC)'] = \
                    eruption.t_end.replace(tzinfo=timezone.utc)
                eruption_results.at[ind, 'Last Notification (UTC)'] = \
                    eruption.last_notification.replace(tzinfo=timezone.utc)
                eruption_results.at[ind, 'Status'] = \
                    str.lower(eruption.status)
                eruption_results.at[ind, 'Estimated Amplitude [Pa]'] = \
                    round(eruption.amp_source)
                eruption_results.at[ind, 'Confidence Level'] = \
                    eruption.confidence_level
                eruption_results.at[ind, 'Revision'] = eruption.revision

            for sta in eruption.data_stations:
                ind_sta_er = assoc_sta_er[assoc_sta_er['Eruption Code'] == eruption.er_code].index
                for ind in ind_sta_er:
                    assoc_sta_er.at[ind, 'Station Name'] = sta.name
                    assoc_sta_er.at[ind, 'Eruption Code'] = eruption.er_code
                    assoc_sta_er.at[ind, 'Num. Detections'] = \
                        sta.nb_dets + sta.nb_dets_previous_analysis
                    assoc_sta_er.at[ind, 'Max. Amp. [Pa]'] = \
                        round(sta.max_amplitude, 5)
                    assoc_sta_er.at[ind, 'Estimated Amp. [Pa]'] = \
                        round(sta.estimated_amplitude)
                    assoc_sta_er.at[ind, 'Detecting'] = sta.detecting

            for volc in eruption.region.volcanoes:
                ind_volc_er = assoc_volc_er[assoc_volc_er['Eruption Code'] == eruption.er_code].index
                for ind in ind_volc_er:
                    assoc_volc_er.at[ind, 'Eruption Code'] = eruption.er_code
                    assoc_volc_er.at[ind, 'Volcano Code'] = volc.code
                    assoc_volc_er.at[ind, 'Volcano Name'] = volc.name

    eruption_results = drop_none_row_from_dataframe(eruption_results, 'Eruption Code')
    eruption_results.to_pickle(results_filename)
    logger.info(f"Saved {results_filename}")

    assoc_sta_er = drop_none_row_from_dataframe(assoc_sta_er, 'Eruption Code')
    assoc_sta_er.to_pickle(assoc_sta_er_filename)
    logger.info(f"Saved {assoc_sta_er_filename}")

    assoc_volc_er = drop_none_row_from_dataframe(assoc_volc_er, 'Eruption Code')
    assoc_volc_er.to_pickle(assoc_volc_er_filename)
    logger.info(f"Saved {assoc_volc_er_filename}")


def save_ip_list_to_db(region: Region, lst_station: list[StationVolc]) -> None:
    """
    Save the IP values for all station monitoring a region to the database.

    Parameters
    ----------
    region : :class:`.Region`
        region monitored.
    lst_station : :obj:`list` of :class:`.StationVolc`
        List of station for which to save IP values.

    """

    ip_results = pd.read_pickle(join(RESULTS_PATH, 'ip_results.pkl'))
    
    tqdm_num = []
    for vol in region.volcanoes:
        for sta in lst_station:
            for ip in sta.list_ip:
                if ip.nb_det > 0:  # NOTE: basically save all 
                    tqdm_num.append([
                                sta.name,
                                vol.code,
                                ip.value,
                                ip.avg_amp,
                                ip.nb_det,
                                ip.vratio,
                                sta.distaz["baz"],
                                sta.distaz["dkm"],
                                ip.dt,
                                ip.avg_source_amp,
                                ip.persistency,
                                ip.fmean
                            ])

    for item in tqdm(tqdm_num, total=len(tqdm_num)):
        this_ip = pd.DataFrame(
            [
                [
                    item[0],  #sta.name,
                    item[1],  #vol.code,
                    item[2],  #ip.value,
                    item[3],  #ip.avg_amp,
                    item[4],  #ip.nb_det,
                    item[5],  #ip.vratio,
                    item[6],  #sta.distaz["baz"],
                    item[7],  #sta.distaz["dkm"],
                    item[8],  #ip.dt,
                    item[9],  #ip.avg_source_amp,
                    item[10], # ip.persistency,
                    item[11], # ip.fmean
                ]
            ],
            columns=ip_results.columns
        )

        ip_results = pd.concat(
            [
                this_ip,
                ip_results
            ],
            ignore_index=True
        )

    #for vol in region.volcanoes:
    #    for sta in lst_station:
    #        for ip in sta.list_ip:
    #            if ip.nb_det >= 0:  # NOTE: basically save all 
    #                this_ip = pd.DataFrame(
    #                    [
    #                        [
    #                            sta.name,
    #                            vol.code,
    #                            ip.value,
    #                            ip.avg_amp,
    #                            ip.nb_det,
    #                            ip.vratio,
    #                            sta.distaz["baz"],
    #                            sta.distaz["dkm"],
    #                            ip.dt,
    #                            ip.avg_source_amp,
    #                            ip.persistency,
    #                            ip.fmean
    #                        ]
    #                    ],
    #                    columns=ip_results.columns
    #                )

    #                ip_results = pd.concat(
    #                    [
    #                        this_ip,
    #                        ip_results
    #                    ],
    #                    ignore_index=True
    #                )

    ip_results = drop_none_row_from_dataframe(ip_results, 'Station Name')

    ip_results.to_pickle(join(RESULTS_PATH, 'ip_results.pkl'))
    logger.info(f"Saved {join(RESULTS_PATH, 'ip_results.pkl')}")


def drop_none_row_from_dataframe(dataframe_input, col_name):
    none_rows = dataframe_input[dataframe_input[col_name] == 'None'].index
    for ind in none_rows:
        dataframe_input = dataframe_input.drop(ind)
    return dataframe_input


def get_filtered_detections_from_db(
    sta: StationVolc,
    t_min: datetime,
    t_max: datetime,
    region: Region,
) -> list[Detection]:
    """
    Get detections from the database for a given :class:`.Station` name.

    Detections are filtered by time and by azimuth.

    Parameters
    ----------
    session : Session
        Session for sqlalchemy local database.
    sta : :class:`.StationVolc`.
        Station at which we get detections.
    t_min : datetime
        Filter detection received after this time.
    t_max : datetime
        Filter detection received before this time.

    Returns
    -------
    :obj:`list` of :class:`.Detection`
        List of detections meeting the criteria.

    """

    detections = pd.read_pickle(join(DATA_PATH, 'detections.pkl'))

    t_min -= timedelta(seconds=IP_TIME_INTERVAL - ANALYSIS_TIME_INTERVAL)

    lst_det = []

    if BAZDEV_PATH is False:
        logger.info(f"Using static back-azimuth tolerance with DAZIM={DAZIM}")
        # Select detections at the given station
        cond1 = detections['station_id'] == sta.name
        # Select by date
        cond2 = (detections['t_start'] >= t_min) & (detections['t_start'] < t_max)
        # Select by frequencies
        #cond3 = (detections['f_mean'] >= MIN_MEAN_FREQ) & (detections['f_mean'] <= 4.0)
        cond3 = (detections['f_mean'] >= MIN_MEAN_FREQ) & (detections['f_mean'] <= MAX_MEAN_FREQ)

        cond4 = None
        if sta.baz_min > sta.baz_max:
            cond4 = (detections['azi'] >= sta.baz_min) | (detections['azi'] <= sta.baz_max)
        else:
            cond4 = (detections['azi'] >= sta.baz_min) & (detections['azi'] <= sta.baz_max)

        lst_det_db = detections[cond1 & cond2 & cond3 & cond4]
        for i in range(len(lst_det_db)):
            det_db = lst_det_db.iloc[i]
            lst_det.append(detection_from_dbo(det_db))
    else:
        logger.info("Using back-azimuth deviations to define back-azimuth tolerance")
        # Get interpolations ===================================================
        volc_name = ''
        if len(region.volcanoes) == 1:
            volc_name = region.volcanoes[0].name
        else:
            # NOTE: this case will not be tested yet
            volc_name = f"Region-{region.name}"
        this_sta_name = get_stat_name(sta.name)  # replaces name if necessary
        this_interp = f"{volc_name[0:4]}_{this_sta_name}.pkl"
        # NOTE: the interpolation format is not final. For now, it's a list
        #       containing the inerpolation data, plus extra information that
        #       helps deciding if using the interpolation values or not.
        bazdevs = pickle.load(open(join(BAZDEV_PATH, this_interp), 'rb'))
        bazdevs_days = bazdevs[0]  # doys from interpolation
        bazdevs_vals = bazdevs[1]  # interp. of average daily back-azimuth deviation
        std_days = bazdevs[7]  # doys from interps. with at least four bazdevs
        std_av = bazdevs[8]  # std of bazdevs with at least four (every 6 hour)
        # strato_arriv = bazdevs[5] # list of booleans if there are strato arrivals
                                    # per day
        n_times = bazdevs[6] # number of arrivals per day

        # Iterate by day =======================================================
        delta_time = t_max - t_min
        this_time = t_min
        for day_i in range(delta_time.days):
            # --------------------------------------
            # Start filter from name and frequencies
            # --------------------------------------
            cond1 = detections['station_id'] == sta.name
            cond3 = (detections['f_mean'] >= MIN_MEAN_FREQ) &\
                    (detections['f_mean'] <= MAX_MEAN_FREQ)
                    #(detections['f_mean'] <= 4.0)

            # ----------------------------
            # Set filter by date (per DOY)
            # ----------------------------
            cond2 = (detections['t_start'] >= this_time) & \
                    (detections['t_start'] < this_time + timedelta(days=1))

            # ----------------------
            # Get DOY from this_time
            # ----------------------
            this_doy = (datetime(this_time.year, this_time.month, this_time.day,
                                 tzinfo=timezone.utc)
                        - datetime(this_time.year, 1, 1, tzinfo=timezone.utc)
                        + timedelta(days=1)).days

            # -------------------------
            # get indexes corresponding
            # to +-1 of the DOY
            # -------------------------
            # Special case below
            if this_doy == 366:  # is a leap year
                this_doy = 365  # use 365, they are similar  # FIXME
            inds = np.where(bazdevs_days < this_doy+1)[0]
            inds2 = np.where(bazdevs_days > this_doy-1)[0]
            lists = [inds, inds2]
            lists = [set(li) for li in lists]
            inds3 = list(lists[0].intersection(*lists))
            this_bazdev = np.mean(bazdevs_vals[inds3]) if len(bazdevs_vals[inds3])>0 else 0
            # below modified (corrected) back-azimuth considering deviation ('expected')
            calc_baz = sta.distaz['baz']-this_bazdev
            # ----------------------------------
            # use STD of bazdev if there
            # are four estimates (every 6 hours)
            # ----------------------------------
            inds_3 = []
            if len(std_days)>0:
                inds_1 = np.where(std_days < this_doy+1)[0]
                inds_2 = np.where(std_days > this_doy-1)[0]
                lists = [inds_1, inds_2]
                lists = [set(li) for li in lists]
                inds_3 = list(lists[0].intersection(*lists))
            this_std = 0
            if n_times[this_doy-1] >= 3:
                this_std = np.mean(std_av[inds_3]) if len(std_av[inds_3]) > 0 else 0

            # ----------------------------------------------
            # get minimum and maximum possible back-azimuths
            # ----------------------------------------------
            # back to +/- DAZIM around corrected back-azimuth
            this_dazim = float(NUM_BAZ_STD)*this_std if this_std > 0 else DAZIM
            # NOTE: testing foced +/- fixed value instead
            if FORCE_DAZIM is True:
                this_dazim = DAZIM
            baz_min = calc_baz - this_dazim
            baz_max = calc_baz + this_dazim
            # wraparound if necessary
            baz_min = baz_min if baz_min > 0 else baz_min+360
            baz_max = baz_max if baz_max < 360 else baz_max-360

            # ---------------------------
            # Set filter by back-azimuths
            # ---------------------------
            cond4 = None
            if baz_min > baz_max:
                cond4 = (detections['azi'] >= baz_min) |\
                        (detections['azi'] <= baz_max)
            else:
                cond4 = (detections['azi'] >= baz_min) &\
                        (detections['azi'] <= baz_max)

            # --------------------------------------------------------------
            # Filter all
            # --------------------------------------------------------------
            lst_det_db = detections[cond1 & cond2 & cond3 & cond4]
            for i in range(len(lst_det_db)):
                det_db = lst_det_db.iloc[i]
                lst_det.append(detection_from_dbo(det_db))


            # --------------------
            # Continue to next day
            # --------------------
            this_time = this_time + timedelta(days=1)

            #if len(lst_det)>0 and sta.name=='I02AR' and this_time.month == 6:
            #    breakpoint()

    return lst_det


def detection_from_dbo(det_db) -> Detection:
    """
    Create :class:`.Detection` from Pandas DataFrame

    Parameters
    ----------
    det_db : DataFrame table

    Returns
    -------
    :class:`.Detection`
        Detection object used for calculations.

    """

    return Detection.from_db(det_db)


def fill_vratio_clim(vratio_file):
    """
    FILLME
    """
    if os.path.exists(vratio_file):
        vratio_data = netCDF4.Dataset(vratio_file)
        num_times = vratio_data.dimensions['t'].size
        num_az = vratio_data.dimensions['az'].size
        all_times = []
        all_veff = []
        # NOTE: assumes resuglar steps (6 h)
        t_0 = num2date(vratio_data['time'][0],
                       units=vratio_data['time'].units,
                       calendar=vratio_data['time'].calendar)
        t_1 = num2date(vratio_data['time'][1],
                       units=vratio_data['time'].units,
                       calendar=vratio_data['time'].calendar)
        t_n = num2date(vratio_data['time'][-1],
                       units=vratio_data['time'].units,
                       calendar=vratio_data['time'].calendar)
        t_step = t_1 - t_0

        t_0 = datetime(START_DATE.year, t_0.month, t_0.day,
                       t_0.hour, t_0.minute, t_0.second)

        t_end = datetime(END_DATE.year, END_DATE.month, END_DATE.day,
                         t_n.hour, t_n.minute, t_n.second)
        # NOTE: year should not be important, as it is about climatologies,
        #       but this is a workaround for a multi-year calculation
        t_ii = 0
        logger.info(f"Importing {vratio_file}...")
        while t_0 + t_step*t_ii < t_end:
            all_veff_circle = []
            for (t_i, az_i) in tqdm(product(range(num_times), range(1, num_az)),
                                    total=num_times*(num_az-1)):
                dt = t_0 + t_step*t_i
                all_veff_circle.append(
                    float(vratio_data['veff_ratio'][t_i, az_i].data)
                )
                if int(az_i) % 360 == 0:
                    all_times.append(dt.replace(tzinfo=timezone.utc))
                    all_veff.append(all_veff_circle)
                    all_veff_circle = []

                if dt >= t_end:  # Don't get more than needed, this is expensive
                    if int(az_i) % 360 == 0:
                        break

            t_0 = datetime(t_0.year+1, t_0.month, t_0.day,
                           t_0.hour, t_0.minute, t_0.second)
            t_ii += 1

        return all_times, all_veff
    else:
        logger.warning(
            f"File {vratio_file} not found!"
        )


def read_OA_BGR_bull(bull_path: str) -> None:
    """
    Load the BGR products from netCDF format.

    See Hupe et al., 2022, Table 1
    [https://doi.org/10.5194/essd-14-4201-2022]

    Parameters
    ----------
    bull_path : str
        Path of the bulletin of the specific station-year.
    """
    # breakpoint()
    if isfile(bull_path):
        data = xr.open_dataset(bull_path, engine='netcdf4')
        if 'N_avail' not in data.sizes:
            logger.error("Empty dataset, skipping")
            return 1
        else:
            data_size = data.sizes['N_avail']
        logger.info(f"Available number of points: {data_size}")
        time_p_aux = data.time_p.to_numpy()
        start_time = []
        end_time = []
        for i in range(time_p_aux.shape[1]):
            t0 = datetime.strptime(f"{netCDF4.chartostring(time_p_aux[:, i])}",
                                   '%Y%m%dT%H%M%S')
            t0 = t0.replace(tzinfo=timezone.utc)
            start_time.append(t0)
            end_time.append(t0+timedelta(minutes=5))  # HF prods
        azim_mean = data.azim.to_numpy()[0, :]
        vapp_mean = data.vapp.to_numpy()[0, :]/1000  # km/s
        amp_mean = data.a_rms.to_numpy()[0, :]
        amp_max = data.a_rms.to_numpy()[2, :]
        freq_mean = data.freq.to_numpy()[0, :]
        fam_size_sum = data.f_size.to_numpy()[2, :]

        data = pd.DataFrame(
            {
                "station_id": get_stat_name(data.instrument.split()[0]),  # NOTE: assumes IMS 4-letter
                "t_start": start_time,
                "t_end": end_time,
                "f_mean": freq_mean,
                "nb_pts": fam_size_sum,
                "amp": amp_mean,
                "amp_max": amp_max,
                "azi": azim_mean,
                "speed": vapp_mean },
            )
        return data
    else:
        logger.error(f"Can't find netCDF file <{bull_path}>")


def load_unifi_data(sta_name: str, bull_path: str) -> None:
    """
    Load UNIFI PMCC data.

    The '.mat' binary contains the following fields.

    ['data']
    - txx: time
    - prs: pressure
    - cc: coherency
    - rs: time residual
    - Fp: peak frequency
    - snr: signal-to-noise ratio
    - azz: azimuth
    - azzsd: azimuth standard deviation
    - slw: apparent velocity
    - slwsd: apparent velocity standard deviation
    - ttfly:
    - sensor_number: sensor ID
    ['info']
    - maxl: 20 — max lag for cross-correlation (in samples)
    - freq: [1 3] — frequency band
    - win: 60 — time window of analysis (seconds)
    - shift: 10 — time shift (seconds)
    - minR: 2 — time residual
    - nmux: 1 — resampling factor (not used)
    - sensors: [1 1 1 1 0 1 1 1] — active sensors '1' ; non active sensor '0'

    Parameters
    ----------
    bull_path : str
        Path of the bulletin of the specific station-year.
    """
    if isfile(bull_path):
        data = sio.loadmat(bull_path)
        data_size = data['data']['txx'][0][0].shape[1]
        times = data['data']['txx'][0][0][0, :]
        print(f"Available number of points: {data_size}")
        start_time = []
        end_time = []
        for i in range(data_size):
            t0 = convert_time_mat2py(times[i])
            start_time.append(t0)
            end_time.append(t0+timedelta(minutes=1))  # time window is 60 s
        azim = data['data']['azz'][0][0][0, :]
        vapp = data['data']['slw'][0][0][0, :]/1000  # km/s
        amp = data['data']['prs'][0][0][0, :]
        freq = data['data']['Fp'][0][0][0, :]
        fam_size_sum = np.ones(shape=freq.shape)
        data = pd.DataFrame(
            {
                "station_id": get_stat_name(sta_name),
                "t_start": start_time,
                "t_end": end_time,
                "f_mean": freq,
                "nb_pts": fam_size_sum,
                "amp": amp,
                "amp_max": -1,
                "azi": azim,
                "speed": vapp
            },
            )
        return data
    else:
        logger.error(f"Can't find file <{bull_path}>")


def load_arise_data(sta_name: str, bull_path: str) -> None:
    """
    Load ARISE-like PMCC data.
    Parse bulletins from 2016 Etna and 2010 Eyjaf. eruptions from
    IMS and local arrays.
    Input is a 18-column text file, with spaces as separators for Etna.
    For Eyjaf., it's a 19-column text file.

    Parameters
    ----------
    bull_path : str
        Path of the bulletin of the specific station-year.
    """

    if isfile(bull_path):
        # format YYYY-MM-DD HH-mm-SS
        data = []
        with open(bull_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.split()
                data.append(line)

        start_time = []
        end_time = []
        freq = []
        fam_size_sum = []
        amp = []
        azim = []
        vapp = []

        for line in data:
            # times
            s1 = str(int(float(line[5])))
            s1 = s1 if s1!='60' else '59'
            s2 = str(int(float(line[8])))
            s2 = s2 if s2!='60' else '59'
            str_tstart = line[0]+'-'+line[1]+'-'+line[2] \
                +' '+line[3]+':'+line[4]+':'+s1
            str_tend = line[0]+'-'+line[1]+'-'+line[2] \
                +' '+line[6]+':'+line[7]+':'+s2
            t_start = datetime.strptime(
                str_tstart, "%Y-%m-%d %H:%M:%S"
                ).replace(tzinfo=timezone.utc)
            t_end = datetime.strptime(
                str_tend, "%Y-%m-%d %H:%M:%S"
                ).replace(tzinfo=timezone.utc)
            start_time.append(t_start)
            end_time.append(t_end)

            freq.append(float(line[10]))
            fam_size_sum.append(float(line[16]))
            amp.append(float(line[11]))
            azim.append(float(line[12]))
            vapp.append(float(line[13]))

        data = pd.DataFrame(
            {
                "station_id": get_stat_name(sta_name),
                "t_start": np.asarray(start_time),
                "t_end": np.asarray(end_time),
                "f_mean": np.asarray(freq),
                "nb_pts": np.asarray(fam_size_sum),
                "amp": np.asarray(amp),
                "amp_max": -1,
                "azi": np.asarray(azim),
                "speed": np.asarray(vapp)
            },
            )
        return data
    else:
        logger.error(f"Can't find file <{bull_path}>")


def load_bgr_data(sta_name, bull_path):
    """
    Load BGR 26-band PMCC data.

    The '.mat' binary contains the following columns:

    a) Year-long (16 columns):
    - 1: start time (MATLAB date number)
    - 2: duration (s)
    - 3: back-azimuth (°)
    - 4: apparent phase velocity (km/s)
    - 5: RMS amplitude (Pa)
    - 6: centre frequency (Hz)
    - 7: min freq. (Hz)
    - 8: max freq. (Hz)
    - 9: family size
    - 10: period at the maximum amplitude (s)
    - 11: pseudo period (s)
    - 12: max. amplitude (Pa)
    - 13: Fisher statistic
    - 14: number of contributing sensors
    - 15: consistency (s)
    - 16: correlation coeﬃcient

    b) Multi-year (10 columns):

    - t_start = convert_time_mat2py(self.data_list[0])
    - t_end = self.t_start + timedelta(seconds=self.data_list[1])
    - azi = self.data_list[2]
    - speed = self.data_list[3]
    - amp = self.data_list[4]
    - f_mean = self.data_list[5]
    - f_min = self.data_list[6]
    - f_max = self.data_list[7]
    - ampmax = self.data_list[9]

    Parameters
    ----------
    bull_path : str
        Path of the bulletin of the specific station-year.
    """
    if isfile(bull_path):
        data = sio.loadmat(bull_path)[f"bull_{sta_name}"]
        num_cols = data.shape[1]
        num_rows = data.shape[0]
        # Col ind 11 corresponds to 'amp_max' in 16-column format
        amp_max_ind = 11
        if num_cols == 10:  # else, it corresponds to 9 in 10-columns format
            amp_max_ind = 9
        start_time = []
        end_time = []
        for i in range(num_rows):
            t0 = convert_time_mat2py(data[i, 0])
            start_time.append(t0)
            end_time.append(t0+timedelta(seconds=data[i, 1]))  # duration
        azim = data[:, 2]
        vapp = data[:, 3]
        amp = data[:, 4]
        amp_max = data[:, amp_max_ind]
        freq = data[:, 5]
        fam_size_sum = data[:, 8]

        data = pd.DataFrame(
            {
                "station_id": sta_name,
                "t_start": start_time,
                "t_end": end_time,
                "f_mean": freq,
                "nb_pts": fam_size_sum,
                "amp": amp,
                "amp_max": amp_max,
                "azi": azim,
                "speed": vapp
            },
            )
        return data
    else:
        print(f"can't find file <{bull_path}>")


def convert_time_mat2py(time: float) -> datetime:
    """
    convert matlab time to datetime.

    parameters
    ----------
    time : float
        matlab time.

    returns
    -------
    datetime
        python time.

    notes
    -----
    see: https://stackoverflow.com/questions/13965740/converting-matlabs-datenum-format-to-python.

    """
    date = datetime.fromordinal(int(time)).replace(
        tzinfo=timezone.utc
    )  # number of full days
    date += timedelta(days=time % 1)  # decimal number of days
    date -= timedelta(days=366)  # substract 366 days
                                # (difference in calendar, see link above)
    return date


def read_veffrat_bgr(vratio_file: str, start_time, end_time):
    mat = sio.loadmat(vratio_file)
    time_ecmwf = []
    veffrat = mat["veffr"]

    for dt in mat["timeperiod"]:
        time_ecmwf.append(convert_time_mat2py(dt[0]))

    time_ecmwf = np.asarray(time_ecmwf)

    time_ind = np.where(time_ecmwf >= start_time)[0]
    time_ecmwf_filt = time_ecmwf[time_ind]
    veffrat_filt = veffrat[:, time_ind]

    time_ind = np.where(time_ecmwf_filt < end_time)[0]
    time_ecmwf_filt = time_ecmwf_filt[time_ind]
    veffrat_filt = veffrat_filt[:, time_ind]

    return time_ecmwf_filt, veffrat_filt


def fill_vratio(vratio_file: str, start_time, end_time):
    if isfile(vratio_file):
        print(f"reading atmospheric data {vratio_file}")
        time_ecmwf, veffrat = read_veffrat_bgr(vratio_file, start_time, end_time)
    else:
        print(
            f"no atmospheric data found for station {vratio_file}."
        )

    return time_ecmwf, veffrat


def load_file(file_to_load, already_processed, file_type):
    """
    Checks if there is a station in list that matches the file and loads it
    if that's the case.
    """

    # Load pre-processed detections if they exist
    detections_in = None
    detections_file = join(DATA_PATH, 'detections.pkl')
    if isfile(detections_file):
        detections_in = pd.read_pickle(detections_file)

    read_func = None
    sta_name = None
    file_name = None
    if file_type == 'OA_BGR':
        read_func = read_OA_BGR_bull
        file_name = file_to_load[1]
    elif file_type == 'BGR':
        read_func = load_bgr_data
        sta_name = file_to_load[0]
        file_name = file_to_load[1]
    elif file_type == 'UNIFI':
        read_func = load_unifi_data
        sta_name = file_to_load[0]
        file_name = file_to_load[1]
    elif file_type == 'ARISE':
        read_func = load_arise_data
        sta_name = file_to_load[0]
        file_name = file_to_load[1]

    # Load detections
    logger.info(f"Loading {file_to_load}...")
    if detections_in is not None:
        if sta_name is None:
            if type(read_func(join(BULL_PATH, file_name))) is not int:
                detections_in = pd.merge(
                        detections_in,
                        read_func(join(BULL_PATH, file_name)),
                        how='outer'
                    )
            else:
                logger.warning(f"Empty dataset for {file_to_load}")
        else:
            detections_in = pd.merge(
                    detections_in,
                    read_func(sta_name, join(BULL_PATH, file_name)),
                    how='outer'
                )
        logger.info("-> Merged into detections.")
        detections_in = detections_in.drop_duplicates()
        detections_in.to_pickle(detections_file)
        logger.info(f"Saved detections in {detections_file}.")
    else:
        logger.info("-> Starting a new table of detections.")
        if sta_name is None:
            detections_in = read_func(join(BULL_PATH, file_name))
        else:
            detections_in = read_func(sta_name, join(BULL_PATH, file_name))

        detections_in = detections_in.drop_duplicates()
        detections_in.to_pickle(detections_file)
        logger.info(f"Saved detections in {detections_file}.")


def load_bulletins(stats_detecting, already_processed):
    # Scan folder with bulletins
    all_files = os.listdir(join(BULL_PATH))
    all_files = list(filter(lambda x: x not in already_processed, all_files))
    logger.info(f"All files: {all_files}")
    sta_names = [f"{stat.name}" for stat in stats_detecting]

    for file_i in all_files:
        file_name_i = file_i.split('.')

        # In case file name has more than one '.'
        if len(file_name_i)>2 and file_name_i[-2] != 'bulletin':
            logger.error("Bulletin file name doesn't follow convention.")
            logger.error("-> File name should have only a '.' for the exension.")
            sys.exit(1)

        if fnmatch(file_i, "IS??_????_hf_1-3Hz_5min.nc"):  # OA BGR
            
            # Get only IMS
            ims_stats = list(filter(lambda x: fnmatch(x, 'I????'), sta_names))
            # Change names from IXYAB -> ISXY
            ims_names = [f"IS{stat[1:3]}" for stat in ims_stats]  # note: only ims

            this_sta = file_i.split('_')[0]
            load = len(list(filter(lambda x: x == this_sta, ims_names))) == 1
            if load is True:
                load_file([this_sta, file_i], already_processed, 'OA_BGR')
                already_processed.append(file_i)

        elif fnmatch(file_i, "bull57_*.mat"):  # 26-band BGR
            # "bull57_I????_????.mat") or "bull57_I????_????-????_VIS.mat"

            this_sta = file_i.split('_')[1]
            load = len(list(filter(lambda x: x == this_sta, sta_names))) == 1
            if load is True:
                load_file([this_sta, file_i], already_processed, 'BGR')
                already_processed.append(file_i)

        elif (fnmatch(file_i, "???_????????_????????.mat") or \
            fnmatch(file_i, "IS??_????????_????????.mat")) is True: # AMT, IS06

            this_sta = file_i.split('_')[0]
            load = len(list(filter(lambda x: x == this_sta, sta_names))) == 1

            if len(this_sta) == 4:  # e.g., IS06
                # Get only IMS
                ims_stats = list(filter(lambda x: fnmatch(x, 'I????'), sta_names))
                ims_names = [f"IS{stat[1:3]}" for stat in ims_stats]
                load = len(list(filter(lambda x: x == this_sta, ims_names))) == 1

            if load is True:
                load_file([this_sta, file_i], already_processed, 'UNIFI')
                already_processed.append(file_i)

        elif fnmatch(file_i, "*.????.bulletin.txt"):

            this_sta = file_i.split('.')[0]
            load = len(list(filter(lambda x: x == this_sta, sta_names))) == 1

            if load is True:
                load_file([this_sta, file_i], already_processed, 'ARISE')
                already_processed.append(file_i)

    return already_processed


def load_data(stats_tab, stats_detecting, volcs_tab, start_time, end_time) -> None:
    """
    uses selected stations to load data into pickle binary from pandas dataframes.
    the dataframes will be:
        - detections -> saved in $data_path/detections.pkl
        - veff-ratios -> saved in $data_path/veff_ratios.pkl
        - stations -> saved in $data_path/stations.pkl
        - volcanoes -> saved in $data_path/volcanoes.pkl
    """
    logger.info("Loading stations data in the database")

    if not isdir(DATA_PATH):
        os.mkdir(DATA_PATH)
        logger.info(f"-> Created directory {DATA_PATH}.")

    # load list of already loaded files to skip
    already_processed = []
    processed_file_path = join(DATA_PATH, 'processed_files.txt')
    try:
        with open(processed_file_path, 'r') as f:
            for line in f.readlines():
                already_processed.append(line.strip())
    except FileNotFoundError:
        logger.warning(f"File {processed_file_path} not found.")

    # ======================================================================
    # detections
    # ======================================================================
    logger.info("Loading detections...")
    load_bulletins(stats_detecting, already_processed)


    # ======================================================================
    # Veff-ratio
    # ======================================================================
    logger.info(f"Loading Veff-ratios [Format: {VEFF_FORMAT}]")
    veff_in = None
    if VEFF_FORMAT == 'BGR':
        all_files = os.listdir(join(VEFF_PATH))
        logger.info(f"All files: {all_files}")
        all_files = list(filter(lambda x: x not in already_processed, all_files))
        sta_names = [f"{stat.name}" for stat in stats_detecting]  # NOTE: only IMS
        files_to_load = []
        sta_names_to_load = []
        for file_i in all_files:
            for sta_name in sta_names:
                if fnmatch(file_i, f"veff50_{sta_name}_????_VIS.mat"):
                    logger.info(f"-> Appended {file_i} for loading.")
                    files_to_load.append(file_i)
                    already_processed.append(file_i)
                    sta_names_to_load.append(sta_name)
                elif fnmatch(file_i, f"veff50_{sta_name}_????-????_VIS.mat"):
                    logger.info(f"-> Appended {file_i} for loading.")
                    files_to_load.append(file_i)
                    already_processed.append(file_i)
                    sta_names_to_load.append(sta_name)
                if files_to_load == []:
                    logger.warning(f"Skipping {file_i} for station {sta_name}")

        veff_file = join(DATA_PATH, 'veff_ratios.pkl')
        if len(files_to_load) > 0:
            if isfile(veff_file):
                logger.info(f"Veff-ratio file {veff_file} exists.")
                veff_in = pd.read_pickle(veff_file)
                logger.info(f"-> Loaded {veff_file}.")

        for file_i, sta_name_i in zip(files_to_load, sta_names_to_load):
            logger.info(f"-> Loading {file_i}...")
            veff_time, veff_ratio = fill_vratio(join(VEFF_PATH, file_i),
                                                start_time,
                                                end_time)
            veff_data = pd.DataFrame(
                    {
                        'dt': veff_time[0],
                        'name': sta_name_i,
                        'azimuth': np.arange(0, 361, 1),
                        'value': veff_ratio[:, 0]
                    }
                )

            for i in tqdm(range(1, len(veff_time)), total=len(veff_time)-1):
                veff_data_aux = pd.DataFrame(
                    {
                        'dt': veff_time[i],
                        'name': sta_name_i,
                        'azimuth': np.arange(0, 361, 1),
                        'value': veff_ratio[:, i]
                    }
                )
                veff_data = pd.merge(veff_data, veff_data_aux, how='outer')

            if veff_in is not None:
                veff_in = pd.merge(veff_in, veff_data, how='outer')
                logger.info("-> Merged into veff-ratio file.")
            else:
                veff_in = veff_data.copy()
    elif VEFF_FORMAT == 'CLIM':
        all_files = os.listdir(join(VEFF_PATH))
        all_files = list(filter(lambda x: x not in already_processed, all_files))
        sta_names = [f"{stat.name}" for stat in stats_detecting]  # NOTE: only IMS
        files_to_load = []
        sta_names_to_load = []
        for file_i in all_files:
            for sta_name in sta_names:
                if fnmatch(file_i, f"????_{sta_name}_veff-ratios.nc"):
                    logger.info(f"-> Appended {file_i} for loading.")
                    files_to_load.append(file_i)
                    already_processed.append(file_i)
                    sta_names_to_load.append(sta_name)
                elif fnmatch(file_i, f"????-????_{sta_name}_veff-ratios.nc"):
                    logger.info(f"-> Appended {file_i} for loading.")
                    files_to_load.append(file_i)
                    already_processed.append[file_i]
                    sta_names_to_load.append(sta_name)

        veff_file = join(DATA_PATH, 'veff_ratios.pkl')
        if len(files_to_load) > 0:
            if isfile(veff_file):
                logger.info(f"Veff-ratio file {veff_file} exists.")
                veff_in = pd.read_pickle(veff_file)
                logger.info(f"-> Loaded {veff_file}.")
        else:
            logger.warning(f"Veff-ratios already loaded. Skipping...")

        for file_i, sta_name_i in zip(files_to_load, sta_names_to_load):
            veff_time, veff_ratio = fill_vratio_clim(join(VEFF_PATH, file_i))
            logger.info(f"-> Saving {file_i}...")
            veff_data = pd.DataFrame(
                    {
                        'dt': veff_time[0],
                        'name': sta_name_i,
                        'azimuth': np.arange(1, 361, 1),
                        'value': veff_ratio[0]
                    }
                )

            for i in tqdm(range(1, len(veff_time)), total=len(veff_time)-1):
                veff_data_aux = pd.DataFrame(
                    {
                        'dt': veff_time[i],
                        'name': sta_name_i,
                        'azimuth': np.arange(1, 361, 1),
                        'value': veff_ratio[i]
                    }
                )
                veff_data = pd.merge(veff_data, veff_data_aux, how='outer')

            if veff_in is not None:
                veff_in = pd.merge(veff_in, veff_data, how='outer')
                logger.info("-> Merged into veff-ratio file.")
            else:
                veff_in = veff_data.copy()
    elif VEFF_FORMAT is False:
        logger.warning("Saved veff-ratio file as 1")
        veff_file = join(DATA_PATH, 'veff_ratios.pkl')
        pickle.dump(1, open(veff_file, 'wb'))
    else:
        logger.error('Wrong VEFF format.')
        sys.exit(0)

    if veff_in is not None:
        veff_in.to_pickle(veff_file)
        print(f"-> Saved {veff_file}.")

    # ======================================================================
    # Stations
    # ======================================================================
    logger.info("Loading stations")
    stats_file = join(DATA_PATH, "stations.pkl")
    if not isfile(stats_file):
        stats_tab.to_pickle(stats_file)
        print(f"-> Saved {stats_file}.")

    # ======================================================================
    # Volcanoes
    # ======================================================================
    logger.info("Loading volcanoes")
    volcs_file = join(DATA_PATH, "volcanoes.pkl")
    if not isfile(volcs_file):
        volcs_tab.to_pickle(volcs_file)
        print(f"-> Saved {volcs_file}.")

    with open(join(DATA_PATH, 'processed_files.txt'), 'w') as f:
        for file_i in already_processed:
            f.write(file_i+'\n')
    logger.info(f"File {join(DATA_PATH, 'processed_files.txt')} saved.")
