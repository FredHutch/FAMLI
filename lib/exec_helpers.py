#!/usr/bin/python
"""Functions to help with execution of system commands."""

import os
import sys
import json
import shutil
import logging
import traceback
import subprocess


def run_cmds(commands, retry=0, catchExcept=False):
    """Run commands and write out the log, combining STDOUT & STDERR."""
    logging.info("Commands:")
    logging.info(' '.join(commands))
    p = subprocess.Popen(commands,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    stdout, stderr = p.communicate()
    exitcode = p.wait()
    if stdout:
        logging.info("Standard output of subprocess:")
        for line in stdout.split('\n'):
            logging.info(line)
    if stderr:
        logging.info("Standard error of subprocess:")
        for line in stderr.split('\n'):
            logging.info(line)

    # Check the exit code
    if exitcode != 0 and retry > 0:
        msg = "Exit code {}, retrying {} more times".format(exitcode, retry)
        logging.info(msg)
        run_cmds(commands, retry=retry - 1)
    elif exitcode != 0 and catchExcept:
        msg = "Exit code was {}, but we will continue anyway"
        logging.info(msg.format(exitcode))
    else:
        assert exitcode == 0, "Exit code {}".format(exitcode)


def get_reference_database(ref_db, temp_folder):
    """Get a reference database folder."""

    # Get files from AWS S3
    if ref_db.startswith('s3://'):
        logging.info("Getting reference database from S3: " + ref_db)

        # Save the database to the local temp folder
        local_fp = os.path.join(
            temp_folder,
            ref_db.rstrip("/").split('/')[-1] + "/"
        )

        assert os.path.exists(local_fp) is False

        logging.info("Saving database to " + local_fp)
        run_cmds([
            'aws',
            's3',
            'sync',
            '--quiet',
            '--sse',
            'AES256',
            ref_db,
            local_fp
        ])

        # Look for the Paladin database within the folder
        for f in os.listdir(local_fp):
            if f.endswith(".pro"):
                return os.path.join(local_fp, f.replace(".pro", ""))

        return local_fp

    else:
        # Treat the input as a local path
        logging.info("Getting reference database from local path: " + ref_db)
        # If the path is a folder, look for the Paladin database within it
        if os.path.isdir(ref_db):
            for f in os.listdir(ref_db):
                if f.endswith(".pro"):
                    return os.path.join(ref_db, f.replace(".pro", ""))
        else:
            # If not, it should be the prefix for the database
            msg = "Please provide Paladin database prefix"
            assert os.path.exists(ref_db + ".pro"), msg

            return ref_db


def return_results(out, read_prefix, output_folder, temp_folder):
    """Write out the final results as a JSON object."""
    # Make a temporary file
    temp_fp = os.path.join(temp_folder, read_prefix + '.json')
    with open(temp_fp, 'wt') as fo:
        json.dump(out, fo)

    # Compress the output
    run_cmds(['gzip', temp_fp])
    temp_fp = temp_fp + '.gz'

    if output_folder.startswith('s3://'):
        # Copy to S3
        run_cmds([
            'aws',
            's3',
            'cp',
            '--quiet',
            '--sse',
            'AES256',
            temp_fp,
            output_folder])
        os.unlink(temp_fp)
    else:
        # Copy to local folder
        run_cmds(['mv', temp_fp, output_folder])


def exit_and_clean_up(temp_folder):
    """Log the error messages and delete the temporary folder."""
    # Capture the traceback
    logging.info("There was an unexpected failure")
    exc_type, exc_value, exc_traceback = sys.exc_info()
    for line in traceback.format_tb(exc_traceback):
        logging.info(line)

    # Delete any files that were created for this sample
    logging.info("Removing temporary folder: " + temp_folder)
    shutil.rmtree(temp_folder)

    # Exit
    logging.info("Exit type: {}".format(exc_type))
    logging.info("Exit code: {}".format(exc_value))
    sys.exit(exc_value)
