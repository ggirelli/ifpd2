"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import argparse
from ifpd2.const import __version__
import joblib  # type: ignore
import logging
from rich.logging import RichHandler  # type: ignore
from rich.console import Console  # type: ignore
import os
import sys
from typing import Optional


def add_version_option(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument(
        "--version", action="version", version=f"{sys.argv[0]} {__version__}"
    )
    return parser


def add_log_file_handler(path: str, logger_name: Optional[str] = None) -> None:
    assert not os.path.isdir(path)
    log_dir = os.path.dirname(path)
    assert os.path.isdir(log_dir) or "" == log_dir
    fh = RichHandler(console=Console(file=open(path, mode="w+")), markup=True)
    fh.setLevel(logging.INFO)
    logging.getLogger(logger_name).addHandler(fh)
    logging.info(f"[green]Log to[/]: '{path}'")


def check_threads(threads: int) -> int:
    if threads > joblib.cpu_count():
        return joblib.cpu_count()
    elif threads <= 0:
        return 1
    return threads
