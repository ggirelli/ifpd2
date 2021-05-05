"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from ifpd2.scripts import arguments, ifpd2
from ifpd2.scripts import db, db_check, db_dump, db_info
from ifpd2.scripts import db_make, db_merge, db_reindex
from ifpd2.scripts import query, query2

import logging
from rich.logging import RichHandler  # type: ignore

logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    handlers=[RichHandler(markup=True, rich_tracebacks=True)],
)

__all__ = [
    "arguments",
    "ifpd2",
    "db",
    "db_check",
    "db_dump",
    "db_info",
    "db_make",
    "db_merge",
    "db_reindex",
    "query",
    "query2",
]
