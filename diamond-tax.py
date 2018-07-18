#!/usr/bin/env python3
"""Wrapper script to run FAMLI on one or more FASTQ files."""

import os
import sys
import uuid
import time
import gzip
import json
import shutil
import logging
import argparse
import subprocess
from famli.exec_helpers import exit_and_clean_up


def get_file_from_url(
    input_str,
    temp_folder
):
    """Get a file from a URL -- return the downloaded filepath."""
    logging.info("Getting file from {}".format(input_str))

    filename = input_str.split('/')[-1]
    local_path = os.path.join(temp_folder, filename)

    logging.info("Filename: " + filename)
    logging.info("Local path: " + local_path)

    if not input_str.startswith(('s3://', 'ftp://')):
        logging.info("Treating as local path")
        msg = "Input file does not exist ({})".format(input_str)
        assert os.path.exists(input_str), msg

        logging.info("Making a symlink to temporary folder")
        os.symlink(input_str, local_path)

    # Get files from AWS S3
    elif input_str.startswith('s3://'):
        logging.info("Getting reads from S3")
        run_cmds([
            'aws', 's3', 'cp', '--quiet', '--sse',
            'AES256', input_str, temp_folder
        ])

    # Get files from an FTP server
    elif input_str.startswith('ftp://'):
        logging.info("Getting reads from FTP")
        run_cmds(['wget', '-P', temp_folder, input_str])

    else:
        raise Exception("Did not recognize prefix for input: " + input_str)

    assert os.path.exists(local_path)
    return local_path


def run_cmds(commands, retry=0, catchExcept=False, stdout=None):
    """Run commands and write out the log, combining STDOUT & STDERR."""
    logging.info("Commands:")
    logging.info(' '.join(commands))
    if stdout is None:
        p = subprocess.Popen(commands,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
        stdout, stderr = p.communicate()
    else:
        with open(stdout, "wt") as fo:
            p = subprocess.Popen(commands,
                                 stderr=subprocess.PIPE,
                                 stdout=fo)
            stdout, stderr = p.communicate()
        stdout = False
    exitcode = p.wait()
    if stdout:
        logging.info("Standard output of subprocess:")
        for line in stdout.decode("utf-8").split('\n'):
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


def diamond_tax(
    input_url=None,
    sample_name=None,
    ref_db=None,
    output_folder=None,
    blocks=5,
    threads=16,
    top_pct=1.0,
    temp_folder="/scratch"
):

    # Make sure that there are no commas or whitespaces in the input
    assert ' ' not in input_url, input_url
    assert ',' not in input_url, input_url
    assert '+' not in input_url, input_url

    # Make a temporary folder for all files to be placed in
    assert os.path.exists(temp_folder)
    temp_folder = os.path.join(temp_folder, str(uuid.uuid4())[:8]) + "/"
    assert os.path.exists(temp_folder) is False
    os.mkdir(temp_folder)

    # Set up logging
    log_fp = os.path.join(temp_folder, "log.txt")
    logFormatter = logging.Formatter(
        '%(asctime)s %(levelname)-8s [DIAMOND-tax] %(message)s'
    )
    rootLogger = logging.getLogger()
    rootLogger.setLevel(logging.INFO)

    # Write to file
    fileHandler = logging.FileHandler(log_fp)
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)
    # Also write to STDOUT
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)

    # Check to see if DIAMOND is available
    logging.info("Checking for a working copy of DIAMOND")
    run_cmds(["diamond", "--version"])

    # Get the reference database
    try:
        db_fp = get_file_from_url(
            ref_db,
            temp_folder
        )
    except:
        exit_and_clean_up(temp_folder)

    logging.info("Reference database: " + db_fp)

    # Align the input data and calculate the LCA per-query

    # Keep track of the time elapsed to process this sample
    start_time = time.time()

    logging.info("Processing input argument: " + input_url)

    # Capture each command in a try statement
    # Get the input reads
    try:
        input_fp = get_file_from_url(input_url, temp_folder)
    except:
        exit_and_clean_up(temp_folder)

    # Run the alignment
    align_fp = os.path.join(temp_folder, sample_name + ".diamond.tax.gz")
    try:
        run_cmds([
            "diamond", "blastp", "--db", db_fp, "--query", input_fp,
            "--out", align_fp, "--outfmt", "102", "--top", str(top_pct),
            "-b", str(blocks), "--threads", str(threads), "--compress", "1"
        ])
    except:
        exit_and_clean_up(temp_folder)

    # Read in the logs
    logging.info("Reading in the logs")
    logs = open(log_fp, 'rt').readlines()

    output = {
        "input_path": input_url,
        "sample": sample_name,
        "output_folder": output_folder,
        "logs": logs,
        "ref_db": db_fp,
        "ref_db_url": ref_db,
        "time_elapsed": time.time() - start_time,
        "params": {
            "blocks": blocks,
            "threads": threads,
            "top_pct": top_pct
        }
    }

    # Make sure the output folder ends with a trailing slash
    if output_folder[-1] != "/":
        output_folder += "/"

    # Make a temporary file
    temp_fp = os.path.join(temp_folder, sample_name + '.json')
    with open(temp_fp, 'wt') as fo:
        json.dump(output, fo)

    # Compress the output
    run_cmds(['gzip', temp_fp])
    temp_fp = temp_fp + '.gz'

    for file_to_upload in [temp_fp, align_fp]:

        if output_folder.startswith('s3://'):
            # Copy to S3
            run_cmds([
                'aws',
                's3',
                'cp',
                '--quiet',
                '--sse',
                'AES256',
                file_to_upload,
                output_folder])
            os.unlink(temp_fp)
        else:
            # Copy to local folder
            run_cmds(['mv', file_to_upload, output_folder])

    # Delete any files that were created for this sample
    logging.info("Removing temporary folder: " + temp_folder)
    shutil.rmtree(temp_folder)

    # Stop logging
    logging.info("Done")
    logging.shutdown()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""
        Align a set of protein sequences against a database and return the LCA per-query.
    """)

    parser.add_argument("--input",
                        type=str,
                        dest="input_url",
                        required=True,
                        help="""Location for input file(s).
                                (Supported: sra://, s3://, or ftp://).""")
    parser.add_argument("--sample-name",
                        type=str,
                        required=True,
                        help="""Name of sample, sets output filename.""")
    parser.add_argument("--ref-db",
                        type=str,
                        required=True,
                        help="""Path to DIAMOND reference database.
                                (Supported: s3://, ftp://, or local path).
                                """)
    parser.add_argument("--output-folder",
                        type=str,
                        required=True,
                        help="""Folder to place results.
                                (Supported: s3://, or local path).""")
    parser.add_argument("--blocks",
                        type=int,
                        default=5,
                        help="""Number of blocks used when aligning.
                                Value relates to the amount of memory used.
                                Roughly 6Gb RAM used by DIAMOND per block.
                                """)
    parser.add_argument("--threads",
                        type=int,
                        default=16,
                        help="Number of threads to use aligning.")
    parser.add_argument("--top-pct",
                        type=float,
                        default=1.0,
                        help="Keep all hits with X percent of the top hit.")
    parser.add_argument("--temp-folder",
                        type=str,
                        default='/share',
                        help="Folder used for temporary files.")

    args = parser.parse_args(sys.argv[1:])

    diamond_tax(**args.__dict__)
