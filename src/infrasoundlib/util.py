"""
This module defines mathematical functions used by the VIS system.
"""

import math
from datetime import datetime

import numpy as np


def interval_processing(
    array: list[datetime], interval: list[datetime]
) -> list[datetime]:
    """
    Given a list of interval (interval are from one value to the other ex: [1,2,4,5] means one
    interval from 1 to 2 and one interval from 4 to 5) and one interval (ex:[4, 3]) return a list
    of interval without any overlap.

    The format of the list of interval was chosen in order to be able to use a binary search.
    The interval must be in sorted in descending order (more practical to insert in a list).

    Parameters
    ----------
    array : :obj:`list` of datetime
        List of intervals.
    interval : :obj:`list` of datetime
        Interval to insert in the list.

    Returns
    -------
    :obj:`list` of datetime
        New list of intervals without any overlap.

    """

    # TODO: Check type here (maybe unit testing would be a good idea)
    indices = np.searchsorted(array, interval, side="right")
    for index, value in zip(indices, interval):
        if index % 2 == 0:
            array.insert(index, value)
    if indices[1] % 2 == 1:
        ctr = 0
    else:
        ctr = 1
    del array[indices[1] + ctr : indices[0] + ctr]
    return array


def dist_az(lat1: float, lon1: float, lat2: float, lon2: float) -> dict[str, float]:
    """
    Compute the great circle arc distance in km and degrees,
    the azimuth in degrees and the back-azimuth in degrees between two positions.

    Parameters
    ----------
    lat1 : float
        Latitude of first point (+N, -S) in degrees.
    lon1 : float
        Longitude of first point (+E, -W) in degrees.
    lat2 : float
        Latitude of second point.
    lon2 : float
        Longitude of second point.

    Returns
    -------
    (float, float, float, float)
        Returns the great circle arc distance in km,
        great circle arc distance in degrees,
        azimuth from pt. 1 to pt. 2 in degrees, and
        azimuth from pt. 2 to pt. 1 in degrees.

    """

    rad = math.pi / 180.0
    sph = 1.0 / 298.257

    scolat = math.pi / 2.0 - math.atan((1.0 - sph) * (1.0 - sph) * math.tan(lat1 * rad))
    ecolat = math.pi / 2.0 - math.atan((1.0 - sph) * (1.0 - sph) * math.tan(lat2 * rad))
    slon = lon1 * rad
    elon = lon2 * rad

    a = math.sin(scolat) * math.cos(slon)
    b = math.sin(scolat) * math.sin(slon)
    c = math.cos(scolat)
    d = math.sin(slon)
    e = -math.cos(slon)
    g = -c * e
    h = c * d
    k = -math.sin(scolat)

    aa = math.sin(ecolat) * math.cos(elon)
    bb = math.sin(ecolat) * math.sin(elon)
    cc = math.cos(ecolat)
    dd = math.sin(elon)
    ee = -math.cos(elon)
    gg = -cc * ee
    hh = cc * dd
    kk = -math.sin(ecolat)

    delrad = math.acos(a * aa + b * bb + c * cc)
    ddeg = delrad / rad
    dkm = ddeg * 111.135

    rhs1 = (aa - d) * (aa - d) + (bb - e) * (bb - e) + cc * cc - 2.0
    rhs2 = (aa - g) * (aa - g) + (bb - h) * (bb - h) + (cc - k) * (cc - k) - 2.0
    dbaz = math.atan2(rhs1, rhs2)

    if dbaz < 0:
        dbaz = dbaz + 2 * math.pi

    bazi = dbaz / rad

    rhs1 = (a - dd) * (a - dd) + (b - ee) * (b - ee) + c * c - 2.0
    rhs2 = (a - gg) * (a - gg) + (b - hh) * (b - hh) + (c - kk) * (c - kk) - 2.0
    dazi = math.atan2(rhs1, rhs2)
    if dazi < 0:
        dazi = dazi + 2 * math.pi

    azi = dazi / rad

    if abs(bazi - 360.0) < 0.00001:
        bazi = 0.0
    if abs(azi - 360.0) < 0.00001:
        azi = 0.0

    return {"dkm": dkm, "ddeg": ddeg, "baz": bazi, "az": azi}


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
