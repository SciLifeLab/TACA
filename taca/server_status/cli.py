import json
import logging

import click

from taca.server_status import (
    cronjobs as cj,  # to avoid similar names with command, otherwise exception
)
from taca.server_status import server_status as status
from taca.utils.config import CONFIG


@click.group(name="server_status")
def server_status():
    """Monitor server status"""


# server status subcommands
@server_status.command()
@click.option("--no_update", is_flag=True, help="Do not update the statusdb")
def nases(no_update):
    """Checks the available space on all the nases"""
    if not CONFIG.get("server_status", ""):
        logging.warning("Configuration missing required entries: server_status")
    disk_space = status.get_nases_disk_space()
    if not no_update:
        status.update_status_db(disk_space, server_type="nas")
    else:
        print(json.dumps(disk_space, indent=4))


@server_status.command()
def cronjobs():
    """Monitors cronjobs and updates statusdb"""
    cj.update_cronjob_db()


@server_status.command()
def monitor_promethion():
    """Checks the status of PromethION and if ngi-nas is mounted"""
    if not CONFIG.get("promethion_status", ""):
        logging.warning("Configuration missing required entries: server_status")
    promethion_status = status.check_promethion_status()
    if promethion_status:
        logging.info("No issues encountered with the PromethION")
    else:
        logging.warning(
            "An issue with the PromethION was encountered. Operator has been notified by email."
        )
