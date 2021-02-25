"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import argparse
from ifpd2.const import __version__
import multiprocessing as mp
import sys


def add_version_option(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument(
        "--version", action="version", version=f"{sys.argv[0]} {__version__}"
    )
    return parser


def check_threads(threads: int) -> int:
    if threads > mp.cpu_count():
        return mp.cpu_count()
    elif threads <= 0:
        return 1
    return threads
