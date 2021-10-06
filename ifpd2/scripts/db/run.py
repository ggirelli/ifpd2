"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import click  # type: ignore
from ifpd2.const import CONTEXT_SETTINGS
from ifpd2.scripts.db import check, dump, info, make


@click.group(
    name="db",
    context_settings=CONTEXT_SETTINGS,
    help="Tools to manage ifpd2 databases.",
)
def main() -> None:
    pass


main.add_command(check.main)
main.add_command(dump.main)
main.add_command(info.main)
main.add_command(make.main)
