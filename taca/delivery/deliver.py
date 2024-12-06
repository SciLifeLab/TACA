"""Delivery methods for TACA."""

import logging
import os
import sys

from taca.delivery.delivery_classes import get_staging_object, get_upload_object
from taca.utils.config import CONFIG

logger = logging.getLogger(__name__)


def locate_data_dirs(project):
    """Locate all data locations of a given project."""
    data_locations = CONFIG.get("delivery").get("data_locations")
    project_locations = []
    logger.info(f"Locating {project} directories.")
    for location in data_locations:
        if os.path.isdir(os.path.join(location, project)):
            project_locations.append(os.path.join(location, project))
            logger.info(
                f"Located {project} directory in {os.path.basename(location)}, including it in staging."
            )
        else:
            logger.info(
                f"No directory found for {project} in {os.path.basename(location)}."
            )
    return project_locations


def stage(project, flowcells, samples):
    """Determine data type and stage accordingly."""
    project_data_dirs = locate_data_dirs(project)

    if not project_data_dirs:
        logger.warning(f"Could not find any data to stage for {project}. Exiting.")
        sys.exit()

    for project_dir in project_data_dirs:
        staging_object = get_staging_object(project, project_dir, flowcells, samples)
        staging_object.stage_data()
    # future todo: update statusdb with status "staged" (project, FC or sample level? Maybe new delivery DB?)


def upload_to_dds(
    project,
    stage_dir,
    pi_email=None,
    add_user=None,
    project_description=None,
    ignore_orderportal_members=False,
):
    """Upload staged data to DDS"""
    upload_object = get_upload_object(
        project,
        stage_dir,
        pi_email,
        add_user,
        project_description,
        ignore_orderportal_members,
    )
    dds_project_id = upload_object.create_dds_project()
    delivery_status = upload_object.upload_data(dds_project_id)
    if delivery_status:
        logger.info(
            f"Successfully uploaded {stage_dir} to DDS project {dds_project_id}"
        )
        # Future todo: Update statusdb with status "uploaded" and DDS project ID
    else:
        logger.error(
            f"Something went wrong when uploading data to {dds_project_id} for project {project}."
        )


def release_dds_project(project, dds_id):
    """Release DDS project to user"""
    # Release DDS project
    # Update statusdb with status "delivered"
    pass
