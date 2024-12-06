"""Delivery classes for TACA."""

import glob
import json
import logging
import os
import re
import subprocess
from datetime import datetime

import requests

from taca.utils.config import CONFIG
from taca.utils.filesystem import do_symlink
from taca.utils.statusdb import StatusdbSession

logger = logging.getLogger(__name__)


def get_staging_object(project, project_dir, flowcells, samples):
    """Identify type of data and instantiate appropriate staging object."""
    if "ONT_TAR" in project_dir:
        return StageTar(project, project_dir, flowcells)
    elif "DATA" in project_dir:
        return StageData(project, project_dir, flowcells, samples)
    elif "ANALYSIS" in project_dir:
        return StageAnalysis(project, project_dir, flowcells, samples)


def get_upload_object(
    project,
    stage_dir,
    pi_email=None,
    add_user=None,
    project_description=None,
    ignore_orderportal_members=False,
):
    """Instantiate upload object."""
    # Future todo: determine data type and instantiate appropriate object
    return UploadNanopore(
        project,
        stage_dir,
        pi_email,
        add_user,
        project_description,
        ignore_orderportal_members,
    )


class Stage:
    """Defines a generic staging object."""

    def __init__(self, project, project_dir, flowcells=None, samples=None):
        self.project_id = project
        self.project_dir = project_dir
        self.flowcells = flowcells
        self.samples = samples
        self.staging_path = CONFIG.get("delivery").get("staging_path")

    def make_staging_path(self):
        """Create a staging dir."""
        timestamp = datetime.now().strftime("%Y%m%d%H%M")
        self.project_staging_path = os.path.join(
            self.staging_path, self.project_id, timestamp
        )
        try:
            os.makedirs(self.project_staging_path)
        except OSError as e:  # future todo: cleanup failed staging and handle rerunning staging with new data
            logger.error(
                f"An error occurred while setting up the staging directory {self.project_staging_path}. Aborting."
            )
            raise (e)


class StageTar(Stage):
    """Defines an object for staging files in ONT_TAR."""

    def __init__(self, project, project_dir, flowcells=None):
        super().__init__(project, project_dir, flowcells)

    def stage_data(self):
        """Stage tarball and md5sum-file to DELIVERY/PID/TIMESTAMP."""
        self.make_staging_path()
        if self.flowcells:
            fc_tarballs = []
            for fc in self.flowcells:
                fc_tarballs.append(os.path.join(self.project_dir, fc + ".tar"))
        else:
            fc_tarballs = glob.glob(os.path.join(self.project_dir, "*.tar"))
        # future todo: look in DB for FCs that have already been delivered/are failed or aborted and exclude them unless listed in -f option
        for tarball_path in fc_tarballs:
            try:
                tarball_filename = os.path.basename(tarball_path)
                do_symlink(
                    tarball_path,
                    os.path.join(self.project_staging_path, tarball_filename),
                )
                do_symlink(
                    tarball_path + ".md5",
                    os.path.join(self.project_staging_path, tarball_filename + ".md5"),
                )
            except FileExistsError as e:
                logger.error(
                    f"An error occurred while symlinking {tarball_path} to {self.project_staging_path}."
                )
                raise (e)
        # future todo: symlink reports
        logger.info(
            f"Successfully staged {self.project_id} to {self.project_staging_path}"
        )


class StageData(Stage):
    """Defines an object for staging files in DATA."""

    def __init__(self, project, project_dir, flowcells=None, samples=None):
        super().__init__(project, project_dir, flowcells, samples)

    def stage_data(self):
        """Stage relevant files in DATA."""
        # Make a list of samples that should be delivered
        # Generate md5sums
        # Symlink all data to DELIVERY/PID/TIMESTAMP
        # Make xml files
        pass


class StageAnalysis(Stage):
    """Defines an object for staging files in ANALYSIS."""

    def __init__(self, project, project_dir, flowcells=None, samples=None):
        super().__init__(project, project_dir, flowcells, samples)

    def stage_data(self):
        """Stage relevant files in ANALYSIS."""
        # For BP, Make a list of samples that should be delivered
        # Generate md5sums
        # Symlink all data and reports to DELIVERY/PID/DDSID
        pass


class Upload:
    """Defines a generic upload object."""

    def __init__(
        self,
        project,
        stage_dir,
        pi_email=None,
        add_user=None,
        project_description=None,
        ignore_orderportal_members=False,
    ):
        self.project_id = project
        self.stage_dir = stage_dir
        self.config_statusdb = CONFIG.get("statusdb", None)
        self.orderportal = CONFIG.get("order_portal", None)
        self.status_db_connection = StatusdbSession(self.config_statusdb)
        self.order_details = self.get_order_details()
        self.pi_email = self.get_pi_email(pi_email)
        self.other_member_details = self.get_other_member_details(
            add_user, ignore_orderportal_members
        )
        self.project_description = self.get_project_description(project_description)

    def get_order_details(self):
        """Fetch order details from order portal"""
        projects_db = self.status_db_connection.connection["projects"]
        view = projects_db.view("order_portal/ProjectID_to_PortalID")
        rows = view[self.project_id].rows
        if len(rows) < 1:
            raise AssertionError(f"Project {self.project_id} not found in StatusDB")
        if len(rows) > 1:
            raise AssertionError(
                f"Project {self.project_id} has more than one entry in StatusDB orderportal_db"
            )
        portal_id = rows[0].value
        # Get project info from order portal API
        get_project_url = "{}/v1/order/{}".format(
            self.orderportal.get("orderportal_api_url"), portal_id
        )
        headers = {
            "X-OrderPortal-API-key": self.orderportal.get("orderportal_api_token")
        }
        response = requests.get(get_project_url, headers=headers)
        if response.status_code != 200:
            raise AssertionError(
                "Status code returned when trying to get "
                "project info from the order portal: "
                f"{portal_id} was not 200. Response was: {response.content}"
            )
        return json.loads(response.content)

    def get_pi_email(self, given_pi_email):
        """Determine the PI email."""
        if given_pi_email:
            logger.warning(
                f"PI email for project {self.project_id} specified by user: {given_pi_email}"
            )
            return given_pi_email
        else:
            found_pi_email = self.order_details["fields"]["project_pi_email"]
            logger.info(
                f"PI email for project {self.project_id} found: {found_pi_email}"
            )
            return found_pi_email

    def get_other_member_details(
        self, other_member_emails=[], ignore_orderportal_members=False
    ):
        """Set other contact details if available. This is not mandatory so
        the method will not raise error if it could not find any contact
        """
        other_member_details = []
        if not ignore_orderportal_members:
            logger.info("Fetching additional members from order portal.")
            try:
                owner_email = self.order_details.get("owner", {}).get("email")
                if (
                    owner_email
                    and owner_email != self.pi_email
                    and owner_email not in other_member_details
                ):
                    other_member_details.append(owner_email)
                bioinfo_email = self.order_details.get("fields", {}).get(
                    "project_bx_email"
                )
                if (
                    bioinfo_email
                    and bioinfo_email != self.pi_email
                    and bioinfo_email not in other_member_details
                ):
                    other_member_details.append(bioinfo_email)
                lab_email = self.order_details.get("fields", {}).get(
                    "project_lab_email"
                )
                if (
                    lab_email
                    and lab_email != self.pi_email
                    and lab_email not in other_member_details
                ):
                    other_member_details.append(lab_email)
            except (AssertionError, ValueError):
                pass  # nothing to worry, just move on
        if other_member_emails:
            logger.info(
                "Other appropriate contacts were found, they will be added "
                f"to DDS delivery project: {', '.join(other_member_emails)}"
            )
            for email in other_member_emails:
                if email not in other_member_details:
                    other_member_details.append(email)
        return other_member_details

    def get_project_description(self, given_desc=None):
        """Set project description, either given or from order portal"""
        if given_desc:
            logger.warning(
                f"Project description for project {self.project_id} specified by user: {given_desc}"
            )
            return given_desc
        else:
            project_name = self.order_details.get("fields", {}).get("project_ngi_name")
            created_desc = f"{project_name} ({datetime.now().strftime('%Y-%m-%d')})"
            logger.info(
                f"Project description for project {self.project_id}: {created_desc}"
            )
            return created_desc

    def create_dds_project(self):
        """Create a DDS delivery project and return the ID."""
        dds_command = [
            "dds",
            "--no-prompt",
            "project",
            "create",
            "--title",
            self.project_id,
            "--description",
            self.project_description,
            "--principal-investigator",
            self.pi_email,
            "--owner",
            self.pi_email,
        ]
        if self.other_member_details:
            for member in self.other_member_details:
                dds_command.append(f"--researcher {member}")
        dds_project_id = ""
        try:
            output = ""
            for line in self._execute(dds_command):
                output += line
                print(line, end="")
        except subprocess.CalledProcessError as e:
            logger.exception(
                "An error occurred while setting up the DDS delivery project."
            )
            raise e
        project_pattern = re.compile("ngisthlm\d{5}")
        found_project = re.search(project_pattern, output)
        if found_project:
            dds_project_id = found_project.group()
            return dds_project_id
        else:
            raise AssertionError(f"DDS project NOT set up for {self.project_id}")

    def upload_data(self, delivery_id):
        """Upload staged data with DDS"""
        log_dir = os.path.join(
            os.path.dirname(CONFIG.get("log").get("file")), "DDS_logs"
        )
        project_log_dir = os.path.join(log_dir, self.project_id)
        cmd = [
            "dds",
            "--no-prompt",
            "data",
            "put",
            "--mount-dir",
            project_log_dir,
            "--project",
            delivery_id,
            "--source",
            self.stage_dir,
        ]
        try:
            output = ""
            for line in self._execute(cmd):
                output += line
                print(line, end="")
        except subprocess.CalledProcessError as e:
            logger.exception(
                f"DDS upload failed while uploading {self.stage_dir} to {delivery_id}"
            )
            raise e
        if "Upload completed!" in output:
            delivery_status = "uploaded"
        else:
            delivery_status = None
        return delivery_status

    def _execute(self, cmd):
        """Helper function to both capture and print subprocess output.
        Adapted from https://stackoverflow.com/a/4417735
        """
        popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
        yield from iter(popen.stdout.readline, "")
        popen.stdout.close()
        return_code = popen.wait()
        if return_code:
            raise subprocess.CalledProcessError(return_code, cmd)


class UploadNanopore(Upload):
    """Defines an object for uploading Nanopore data."""

    def __init__(
        self,
        project,
        stage_dir,
        pi_email=None,
        add_user=None,
        project_description=None,
        ignore_orderportal_members=False,
    ):
        super().__init__(
            project,
            stage_dir,
            pi_email,
            add_user,
            project_description,
            ignore_orderportal_members,
        )


class UploadIllumina(Upload):
    """Defines an object for uploading Illumina data."""

    def __init__(self, project, stage_dir):
        super().__init__(project, stage_dir)


class UploadElement(Upload):
    """Defines an object for uploading Element data."""

    def __init__(self, project, stage_dir):
        super().__init__(project, stage_dir)
