"""
This module defines the Region class.
"""

import numpy as np

from src.volcano import Volcano


class Region:
    """
    A class to represent a region. It is defined with its top left and bottom
    right corners. It should at least contain one volcano. Its barycenter is
    calculated from the considered volcanoes.

    Attributes
    ----------
    volcanoes : :obj:`list` of :class:`.Volcano`
        List of volcanoes in the region.
    name : str
        Generic name of the region.
    top_left_corner : :obj:`tuple` of float
        Coordinates of the top left corner of the region.
    bot_right_corner : :obj:`tuple` of float
        Coordinates of the bottom right corner of the region.
    lat : float
        Latitude of the center of the region.
    lon : float
        Longitude of the center of the region.

    """

    def __init__(
        self,
        top_left_corner: tuple[float, float] | None,
        bot_right_corner: tuple[float, float] | None,
        volcanoes: list[Volcano],
        idx: int,
    ):
        """
        Constructor of :class:`.Region`.

        Attributes
        ----------
        top_left_corner : :obj:`tuple` of float
            Coordinates of the top left corner of the region.
        bot_right_corner : :obj:`tuple` of float
            Coordinates of the bottom right corner of the region.
        volcanoes : :obj:`list` of :class:`.Volcano`
            List of volcanoes in the region.

        """

        self.top_left_corner = top_left_corner
        self.bot_right_corner = bot_right_corner
        self.volcanoes = volcanoes
        if len(self.volcanoes) > 1:  # a list of Volcano instances
            self.name = f"Region-{idx}"
        else:  # a sole Volcano instance
            self.name = volcanoes[0].name
        self.compute_center()

    def compute_center(self) -> None:
        """
        Compute the geographic midpoint of the volcanoes located in the region.

        """

        x_cart = []
        y_cart = []
        z_cart = []
        for vol in self.volcanoes:
            lat = vol.lat * np.pi / 180
            lon = vol.lon * np.pi / 180
            x_cart.append(np.cos(lat) * np.cos(lon))
            y_cart.append(np.cos(lat) * np.sin(lon))
            z_cart.append(np.sin(lat))
        x_final = np.sum(x_cart) / len(x_cart)
        y_final = np.sum(y_cart) / len(y_cart)
        z_final = np.sum(z_cart) / len(z_cart)
        lon_center = np.arctan2(y_final, x_final)
        hyp = np.sqrt(np.power(x_final, 2) + np.power(y_final, 2))
        lat_center = np.arctan2(z_final, hyp)
        self.lon = np.round(lon_center * 180 / np.pi, 2)
        self.lat = np.round(lat_center * 180 / np.pi, 2)
