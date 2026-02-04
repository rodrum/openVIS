"""
This module defines and create the logger of the VIS.

Two classes are created by inheritance from classes in the logging module. 
One to add colors to the output. The other to prevent repetitive message to be 
logged.
"""

import logging
import sys


class ColorFormatter(logging.Formatter):
    """
    This class inherit from the `logging.Formatter` class.
    It allows log to be outputed with a color
    that depends on the level of the message.

    Notes
    -----
    See https://stackoverflow.com/questions/384076/how-can-i-color-python-logging-output

    """

    green = "\x1b[32;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    format_c = "%(asctime)s %(levelname)-8s [%(filename)-14s:%(lineno)-3d] %(message)s"
    datefmt = "%Y-%m-%d %H:%M:%S"

    FORMATS = {
        logging.INFO: green + format_c + reset,
        logging.WARNING: yellow + format_c + reset,
        logging.ERROR: red + format_c + reset,
        logging.CRITICAL: bold_red + format_c + reset,
    }

    def format(self, record: logging.LogRecord) -> str:
        """
        Get record and format it with colors.

        Parameters
        ----------
        record
            Record to format.

        Returns
        -------
        record
            Formatted record.

        """

        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)


class NonRepetitiveLogger(logging.Logger):
    """
    This class inherit from the logging.Logger class. It prevents same message to be outputed
    multiple times. Instead, another line with the number of times the message was displayed is
    outputed.

    Attributes
    ----------
    last_message : str
        Last message logged.
    counter : int
        Count the number of times the same message was logged.
    last_level : int
        Keep the level of the last message.

    """

    def __init__(self, name: str, level: int = logging.NOTSET):
        """
        Constructor of :class:`.NonRepetitiveLogger`.

        Parameters
        ----------
        name : str
            Name of the logger.
        level : int
            Level of the logger.

        """

        super().__init__(name=name, level=level)
        self.last_message = ""
        self.counter = 0
        self.last_level = logging.NOTSET


    def _log(  # type: ignore
        self, level: int, msg, args, exc_info=None, extra=None, stack_info=False
    ) -> None:
        sys.stdout.write("\r")
        if msg == self.last_message and self.last_level == level:
            self.counter += 1
            return
        if self.counter == 0:
            super()._log(level, msg, args, exc_info, extra, stack_info, stacklevel=2)
            self.last_message = msg
            self.last_level = level
            self.counter = 0
            return
        if self.counter == 1:
            rep_msg = " --- The last message repeated 1 time\033[K"
            super()._log(
                self.last_level,
                rep_msg,
                args,
                exc_info,
                extra,
                stack_info,
                stacklevel=2,
            )
        else:
            rep_msg = f" --- The last message repeated {self.counter} times\033[K"
            super()._log(
                self.last_level,
                rep_msg,
                args,
                exc_info,
                extra,
                stack_info,
                stacklevel=2,
            )
        super()._log(level, msg, args, exc_info, extra, stack_info, stacklevel=2)
        self.last_message = msg
        self.last_level = level
        self.counter = 0


def init_logger() -> NonRepetitiveLogger:
    """
    Function to initialize the VIS logger.

    """

    log = NonRepetitiveLogger(__name__)
    file_h = logging.FileHandler(filename="VIS.log", mode="w")

    console = logging.StreamHandler()

    console_formatter = ColorFormatter()
    format_str = (
        "%(asctime)s,%(msecs)d %(levelname)-8s [%(pathname)s:%(lineno)d in function"
        "%(funcName)s] %(message)s"
    )
    file_formatter = logging.Formatter(format_str, datefmt="%Y-%m-%d %H:%M:%S")
    console.setFormatter(console_formatter)
    file_h.setFormatter(file_formatter)

    log.addHandler(file_h)
    log.addHandler(console)
    log.setLevel(logging.INFO)
    return log


logger = init_logger()
