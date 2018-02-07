#!/usr/bin/python
"""Functions to help working with FASTQ files."""

import os
import gzip
import uuid
import logging
import subprocess
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.FastaIO import SimpleFastaParser
from lib.exec_helpers import run_cmds


def get_reads_from_url(
    input_str,
    temp_folder,
    random_string=str(uuid.uuid4())[:8]
):
    """Get a set of reads from a URL -- return the downloaded filepath."""
    logging.info("Getting reads from {}".format(input_str))

    filename = input_str.split('/')[-1]
    local_path = os.path.join(temp_folder, filename)

    logging.info("Filename: " + filename)
    logging.info("Local path: " + local_path)

    if not input_str.startswith(('s3://', 'sra://', 'ftp://')):
        logging.info("Treating as local path")
        msg = "Input file does not exist ({})".format(input_str)
        assert os.path.exists(input_str), msg
        logging.info("Copying to temporary folder, cleaning up headers")

        # Add a random string to the filename
        local_path = local_path.split('/')
        if local_path[-1].endswith(".gz"):
            local_path[-1] = local_path[-1].replace(".gz", "")
        local_path[-1] = "{}-{}".format(random_string, local_path[-1])
        local_path = '/'.join(local_path)

        # Make the FASTQ headers unique
        clean_fastq_headers(input_str, local_path)

        return local_path

    # Get files from AWS S3
    if input_str.startswith('s3://'):
        logging.info("Getting reads from S3")
        run_cmds([
            'aws', 's3', 'cp', '--quiet', '--sse',
            'AES256', input_str, temp_folder
            ])

    # Get files from an FTP server
    elif input_str.startswith('ftp://'):
        logging.info("Getting reads from FTP")
        run_cmds(['wget', '-P', temp_folder, input_str])

    # Get files from SRA
    elif input_str.startswith('sra://'):
        accession = filename
        logging.info("Getting reads from SRA: " + accession)
        local_path = get_sra(accession, temp_folder)

    else:
        raise Exception("Did not recognize prefix for input: " + input_str)

    # Add a random string to the filename
    new_path = local_path.split('/')
    if new_path[-1].endswith(".gz"):
        new_path[-1] = new_path[-1].replace(".gz", "")
    new_path[-1] = "{}-{}".format(random_string, new_path[-1])
    new_path = '/'.join(new_path)
    logging.info(
        "Copying {} to {}, cleaning up FASTQ headers".format(
            local_path, new_path
            )
        )
    clean_fastq_headers(local_path, new_path)
    logging.info("Deleting old file: {}".format(local_path))
    os.unlink(local_path)
    return new_path


def get_sra(accession, temp_folder):
    """Get the FASTQ for an SRA accession."""
    local_path = os.path.join(temp_folder, accession + ".fastq")

    logging.info("Downloading {} from SRA".format(accession))
    run_cmds([
        "fastq-dump",
        "--split-files",
        "--outdir",
        temp_folder, accession
    ])

    # Combine any multiple files that were found
    logging.info("Concatenating output files")
    with open(local_path + ".temp", "wt") as fo:
        cmd = "cat {}/{}*fastq".format(temp_folder, accession)
        cat = subprocess.Popen(cmd, shell=True, stdout=fo)
        cat.wait()

    # Clean up the FASTQ headers for the downloaded file
    if os.path.exists(local_path + ".temp"):
        run_cmds(["mv", local_path + ".temp", local_path])

    # Check to see if the file was downloaded
    msg = "File could not be downloaded from SRA: {}".format(accession)
    assert os.path.exists(local_path), msg

    # Return the path to the file
    logging.info("Done fetching " + accession)
    return local_path


def count_fasta_reads(fp):
    n = 0
    if fp.endswith(".gz"):
        with gzip.open(fp, "rt") as f:
            for record in SimpleFastaParser(f):
                n += 1
    else:
        with open(fp, "rt") as f:
            for record in SimpleFastaParser(f):
                n += 1

    return n


def count_fastq_reads(fp):
    n = 0
    if fp.endswith(".gz"):
        with gzip.open(fp, "rt") as f:
            for record in FastqGeneralIterator(f):
                n += 1
    else:
        with open(fp, "rt") as f:
            for record in FastqGeneralIterator(f):
                n += 1

    # If no reads were found, try counting it as a FASTA
    if n == 0:
        logging.info("No FASTQ reads found, trying to read as FASTA")
        n = count_fasta_reads(fp)

    return n


def clean_fastq_headers(fp_in, fp_out):
    """Read in a FASTQ file and write out a copy with unique headers."""

    # Constraints
    # 1. Headers start with '@'
    # 2. Headers are stripped to the first whitespace
    # 3. Headers are unique
    # 4. Sequence lines are not empty
    # 5. Spacer lines match the header line
    # 6. Quality lines are not empty

    if fp_in.endswith(".gz"):
        f_in = gzip.open(fp_in, "rt")
    else:
        f_in = open(fp_in, "rt")

    if fp_out.endswith(".gz"):
        f_out = gzip.open(fp_out, "wt")
    else:
        f_out = open(fp_out, "wt")

    # Keep track of the line number
    for ix, line in enumerate(f_in):
        # Get the line position 0-3
        mod = ix % 4

        if mod == 0:
            # Skip lines that are blank (at the end of the file)
            if len(line) == 1:
                continue
            # 1. Headers start with '@'
            assert line[0] == '@', "Header lacks '@' ({})".format(line)

            # 2. Strip to the first whitespace
            line = line.rstrip("\n").split(" ")[0].split("\t")[0]

            # 3. Add a unique line number and the newline
            line = "{}-r{}\n".format(line, 1 + (ix / 4))

            # Save the header to use for the spacer line
            header = line[1:]

        elif mod == 1:
            # 4. Sequence lines are not empty
            assert len(line) > 1

        elif mod == 2:
            # 5. Spacer lines start with '+' and match the header
            assert line[0] == "+"
            line = "+" + header

        elif mod == 3:
            # 6. Quality lines are not empty
            assert len(line) > 1

        # Write out the line
        f_out.write(line)

    # Close the input and output file handles
    f_in.close()
    f_out.close()
