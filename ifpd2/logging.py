"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import logging


class Loggable(object):
    """Shows logger instance to children classes."""

    defaultfmt = "%(asctime)s %(levelname)s:\t%(message)s"
    datefmt = "%d/%m/%Y %H:%M:%S"

    def __init__(
        self,
        logger=logging.getLogger(),
        formatter=logging.Formatter(
            "%(asctime)s %(levelname)s:\t%(message)s", datefmt="%d/%m/%Y %H:%M:%S"
        ),
    ):
        super(Loggable, self).__init__()
        self._logger = logger
        self._formatter = formatter

    @property
    def log(self):
        return self._logger

    def addFileHandler(self, path, mode="w+", level=logging.DEBUG):
        fileHandler = logging.FileHandler(path, mode)
        fileHandler.setFormatter(self._formatter)
        fileHandler.setLevel(level)
        self.log.addHandler(fileHandler)
