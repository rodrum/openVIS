"""
Script to get the filtered detections by amplitude, following the procedure of the VIS.

RDSNL // 2025-08-17 // CF
"""

from os.path import join
import pandas as pd
from datetime import datetime, timezone, timedelta
from obspy.geodetics.base import gps2dist_azimuth
import numpy as np
from tqdm import tqdm


# Get veff-ratio data
sta_name = 'I02AR'
data_path = '/home/rodrigo/Desktop/openVIS/cordon_caulle_2011/compiled_data/'
veff_sta = pd.read_pickle(join(data_path, 'veff_ratios.pkl'))
i02_lat, i02_lon = -54.58, -67.31
cc_lat, cc_lon = -40.59, -72.117
d_i02, _, baz_i02 = gps2dist_azimuth(cc_lat, cc_lon, i02_lat, i02_lon)

def get_vratio(sta_baz, t_start, veff_sta):
    """
    Assuming vratio is a table with one record per time of day (usually every
    6 hours) and angle (0 to 359 in steps of 1 degree).
    Hence dimension Datetimes x 360, whatever Datatimes is.
    In 365 days, four times a day, this would be [365x4 x 360].
    The headers are: 'dt', 'name', 'azimuth' and 'value' (see 'load_data')
    """
    if type(veff_sta) == int:  # a terrible way of doing this, FIXME
        return 1
    vratios = veff_sta[veff_sta['dt'] <= t_start+timedelta(hours=6)]
    vratios = vratios[vratios['dt'] >= t_start]
    vratios = vratios[vratios['azimuth'] == int(round(sta_baz))%360]
    return vratios

# Get detections
detections = pd.read_pickle(join(data_path, 'detections.pkl'))

# =====================
# Calculate attenuation
# =====================

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


def bilinear_interpolation(
    x: float, y: float, points: list[tuple[float, float, float]]
) -> float:
    """
    Interpolate (x,y) from values associated with four points.

    The four points are a list of four triplets:  (x, y, value).
    The four points can be in any order.  They should form a rectangle.


    Parameters
    ----------
    x : float
        x value of the point to interpolate.
    y : float
        y value of the point to interpolate.
    points : :obj:`list`
        Four list of triplets that form a rectangle.

    Returns
    -------
    float
        Interpolated value.


    Example
    -------
    bilinear_interpolation(12, 5.5,
    ...                   [(10, 4, 100),
    ...                   (20, 4, 200),
    ...                   (10, 6, 150),
    ...                   (20, 6, 300)])
    165.0

    Notes
    -----
        See formula at:  http://en.wikipedia.org/wiki/Bilinear_interpolation

    """

    points = sorted(points)  # order points by x, then by y
    (x1, y1, q11), (_x1, y2, q12), (x2, _y1, q21), (_x2, _y2, q22) = points

    if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
        raise ValueError("points do not form a rectangle")
    if not x1 <= x <= x2 or not y1 <= y <= y2:
        raise ValueError("(x, y) not within the rectangle")

    return (
        q11 * (x2 - x) * (y2 - y)
        + q21 * (x - x1) * (y2 - y)
        + q12 * (x2 - x) * (y - y1)
        + q22 * (x - x1) * (y - y1)
    ) / ((x2 - x1) * (y2 - y1) + 0.0)


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

detections_i02 = detections[detections['station_id'] == 'I02AR']
t_start = datetime(2011,1,1,0,0,0).replace(tzinfo=timezone.utc)
t_end = datetime(2012,1,1,0,0,0).replace(tzinfo=timezone.utc)
veff_i02 = veff_sta[veff_sta['name'] == 'I02AR']

att_coeffs = []
dets_to_use = []
max_amp = 500.0
for i, row in tqdm(detections_i02.iterrows(), total=len(detections_i02)):
    freq = row['f_mean']
    amp = row['amp']
    t_start = row['t_start']
    vratio = get_vratio(baz_i02, t_start, veff_i02)['value'].values[0]
    att_coeff = calculate_att_coeff(freq, vratio, d_i02/1e3)
    att_coeffs.append(att_coeff)
    if amp/att_coeff < max_amp:
        dets_to_use.append(row)
breakpoint()
