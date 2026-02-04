"""
This module is used to handle Vratios
"""

from datetime import datetime, timedelta, timezone

import numpy as np

from src.infrasoundlib.vratio import Vratio
from src.logger import logger

from .util import bilinear_interpolation

# Values from Le Pichon et al., Performance of infrasound network (2012)
LST_FREQ = [0.1, 0.2, 0.4, 0.8, 1.6, 3.2]
LST_ALPHA = [-0.28, -0.33, -0.39, -0.47, -0.59, -0.69]
LST_VRATIO = [0.85, 0.88, 0.91, 0.94, 0.97, 1, 1.03, 1.06, 1.09, 1.12, 1.15, 1.18]
LST_BETA = [
    [-1.0, -1.0, -1.0, -1.0, -1.0, -0.9, -0.9, -0.85, -0.85, -0.85, -0.8, -0.8],
    [-1.2, -1.2, -1.15, -1.1, -1.05, -0.95, -0.9, -0.85, -0.85, -0.85, -0.8, -0.8],
    [-1.35, -1.4, -1.4, -1.35, -1.15, -1.0, -0.9, -0.9, -0.9, -0.85, -0.85, -0.85],
    [-1.7, -1.6, -1.6, -1.55, -1.3, -1.05, -0.9, -0.9, -0.9, -0.9, -0.95, -0.9],
    [-1.95, -1.95, -1.85, -1.75, -1.4, -1.15, -0.95, -0.95, -1.0, -1.0, -1.05, -1.05],
    [-2.4, -2.3, -2.1, -1.85, -1.45, -1.2, -1.0, -1.0, -1.05, -1.1, -1.15, -1.2],
]
LST_SIGMA = [79, 55, 43, 36, 27, 20]

DELTA = 180  # Width of shadow zone (km)


def calculate_att_coeff(freq: float, vratio: float, dist: float) -> float:
    """
    Calculate attenuation coefficient.

    Parameters
    ----------
    freq : float
        Frequency.
    vratio : float
        Vratio.
    dist : float
        Distance.

    Returns
    -------
    float
        Attenuation coefficient.
    """
    if dist:
        alpha = np.interp(freq, LST_FREQ, LST_ALPHA)
        if (
            LST_FREQ[0] < freq < LST_FREQ[-1]
            and LST_VRATIO[0] < vratio < LST_VRATIO[-1]
        ):
            # bilinear_interpolation
            xmin, xplus = find_closest_indexes(LST_FREQ, freq)
            ymin, yplus = find_closest_indexes(LST_VRATIO, vratio)
            beta = bilinear_interpolation(
                freq,
                vratio,
                [
                    (LST_FREQ[xmin], LST_VRATIO[ymin], LST_BETA[xmin][ymin]),
                    (LST_FREQ[xplus], LST_VRATIO[ymin], LST_BETA[xplus][ymin]),
                    (LST_FREQ[xmin], LST_VRATIO[yplus], LST_BETA[xmin][yplus]),
                    (LST_FREQ[xplus], LST_VRATIO[yplus], LST_BETA[xplus][yplus]),
                ],
            )
        else:
            i = find_closest(LST_FREQ, freq)
            beta = np.interp(vratio, LST_VRATIO, LST_BETA[i])
        sigma = np.interp(freq, LST_FREQ, LST_SIGMA)

        near_field_contribution = 10 ** ((alpha * dist) / 20) / dist
        far_field_contribution = np.power(dist, beta) / (
            1 + 10 ** ((DELTA - dist) / sigma)
        )
        return near_field_contribution + far_field_contribution
    logger.warning("Attenuation coefficient cannot be calculated ")
    return 0


def find_closest(arr: list[float], target: float) -> int:
    """
    Find index of the closest value of target in array.

    Parameters
    ----------
    arr : :obj:`list` of float
        Sorted array of float.
    target : float
        Target value.

    Returns
    -------
    int
        Index of closest value of target in array.

    """

    idx = np.asarray(arr).searchsorted(target)
    idx = np.clip(idx, 1, len(arr) - 1)
    left = arr[idx - 1]
    right = arr[idx]
    idx -= target - left < right - target
    return int(idx)


def find_closest_indexes(arr: list[float], target: float) -> tuple[int, int]:
    """
    Find indexes of the closest values of target in array.

    Parameters
    ----------
    arr : :obj:`list` of float
        Sorted array of float.
    target : float
        Target value.

    Returns
    -------
    (int, int)
        Indexes of closest values of target in array.

    """

    idx = np.asarray(arr).searchsorted(target)
    idx = np.clip(idx, 1, len(arr) - 1)
    return int(idx - 1), int(idx)


def interpolate_vratio(
    station_vratio: list[Vratio] | None,
    ref_time: datetime,
    vratio_time_interval: int = 21600,
) -> float:
    """
    Interpolate vratio value.

    Parameters
    ----------
    station_vratio : :obj:`list` of VratioDB
        List of VratioDB object queried from the database.
    ref_time : datetime
        Reference time at which to get vratio.
    azimuth : float
        Azimuth at which to get vratio.
    vratio_time_interval : int
        Number of seconds before and after ref_time to look for vratio.

    Returns
    -------
    float
        Vratio value at given time and azimuth.
    """
    vratio_available = False
    vratio = 1.0
    if station_vratio:
        i_dt = [
            i
            for i, v in enumerate(station_vratio)
            if v.dt >= ref_time - timedelta(seconds=vratio_time_interval)
            and v.dt <= ref_time + timedelta(seconds=vratio_time_interval)
        ]
        if len(i_dt) >= 2:
            lst_dt = [station_vratio[i].dt.timestamp() for i in i_dt]
            lst_vratio = [station_vratio[i].value for i in i_dt]
            vratio = np.interp(ref_time.timestamp(), lst_dt, lst_vratio)
            vratio_available = True
        elif len(i_dt) == 1:
            vratio = station_vratio[i_dt[0]].value
            vratio_available = True

    if not vratio_available:
        logger.warning(
            "ECMWF data not available for the concerned period. Vratio = 1 used instead"
        )
    return vratio


def convert_time_mat2py(time: float) -> datetime:
    """
    Convert matlab time to datetime.

    Parameters
    ----------
    time : float
        Matlab time.

    Returns
    -------
    datetime
        Python time.

    Notes
    -----
    See: https://stackoverflow.com/questions/13965740/converting-matlabs-datenum-format-to-python.

    """
    date = datetime.fromordinal(int(time)).replace(
        tzinfo=timezone.utc
    )  # Number of full days
    date += timedelta(days=time % 1)  # Decimal number of days
    # Substract 366 days  (difference in calendar, see link above)
    date -= timedelta(days=366)
    return date


def convert_time_bgrprod2py(time_p: str) -> datetime:
    """
    Convert BGR products time to datetime

    Parameters
    ----------
    time_p: datetime string [yyyymmddTHHMMSS (ISO 8601)

    Returns
    -------
    datetime
        Python time.
    """
    date = datetime.strptime(time_p, '%Y%m%dT%H%M%S')
    return date
