"""
This module defines the ProgressBar class.

It is used to display the progress of the analysis, and the estimated remaining time.
"""

import sys
import time


class ProgressBar:
    """
    A class to plot and handle a progress bar that shows the progress of the analysis.

    Attributes
    ----------
    start : int
        Number at which the progress bar starts.
    end : int
        Number at which the progress bar stops.
    start_time : datetime
        Time at which the progress bar started.
    now_time : datetime
        Time now, when the progress bar is displayed.

    Notes
    -----
    External progress bar module obtained from
    StackOverflow <http://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console>

    """

    # pylint: disable=too-many-instance-attributes

    DEFAULT_BAR_LENGTH = float(65)

    def __init__(self, end: int, start: int = 0):
        """
        Constructor of :class:`.ProgressBar`.

        Parameters
        ----------
        end : int
            Number at which the progress bar ends.
        start : int, optional
            Number at which the progress bar starts. Defaults to 0.

        """

        if start == end:
            self.end = end + 1
        else:
            self.end = end
        self.start = start
        self._bar_length = ProgressBar.DEFAULT_BAR_LENGTH
        self._ratio = 0.0
        self._level_chars = 0
        self.start_time = time.time()
        self.now_time = time.time()

        self.set_level(self.start)
        self._plotted = False

    def set_level(self, level: int) -> None:
        """
        Set the level of the progress bar.

        Parameters
        ----------
        level : int
            Level of the progress bar.

        """

        _level = level
        if level < self.start:
            _level = self.start
        if level > self.end:
            _level = self.end

        self._ratio = float(_level - self.start) / float(self.end - self.start)
        self._level_chars = int(self._ratio * self._bar_length)

    def plot_progress(self) -> None:
        """
        Plot the progress bar.

        """

        sys.stdout.write(
            "\r  %3i%% [%s%s] -- estimated remaining time: %8.2f seconds"
            % (
                int(self._ratio * 100.0),
                "=" * int(self._level_chars),
                " " * int(self._bar_length - self._level_chars),
                (self.now_time - self.start_time) * (1 / self._ratio - 1),
            )
        )
        sys.stdout.flush()
        self._plotted = True

    def set_and_plot(self, level: int) -> None:
        """
        Set the level of the progress bar and plot it.

        Parameters
        ----------
        level : int
            Level of the progress bar.

        """

        old_chars = self._level_chars
        self.set_level(level)
        self.now_time = time.time()
        if (not self._plotted) or (old_chars != self._level_chars):
            self.plot_progress()

    def __del__(self) -> None:
        sys.stdout.write("\n")
