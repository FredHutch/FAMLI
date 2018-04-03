#!/usr/bin/python
"""Functions to help working with FASTQ files."""

import re
import os
import gzip
import uuid
import shutil
import logging
import subprocess
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.FastaIO import SimpleFastaParser
from famli.exec_helpers import run_cmds


def combine_fastqs(fps_in, fp_out):
    """Combine multiple FASTQs into a single FASTQ."""
    assert len(fps_in) > 0

    if len(fps_in) == 1:
        assert os.path.exists(fps_in[0])
        logging.info("Making a symlink: {} -> {}".format(fps_in[0], fp_out))
        os.symlink(fps_in[0], fp_out)
    else:
        logging.info("Combining {:,} FASTQ files".format(len(fps_in)))
        logging.info("Writing all inputs to {}".format(fp_out))
        with open(fp_out, "wt") as fo:
            for fp_ix, f in enumerate(fps_in):
                logging.info("Adding {} to {}".format(f, fp_out))
                with open(f, "rt") as fi:
                    for line_ix, line in enumerate(fi):
                        # Add a file index to the header
                        # In case there are two files with the same headers
                        mod = line_ix % 4
                        if mod == 0 or mod == 2:
                            line = line.rstrip("\n")
                            fo.write("{}-{}\n".format(line, fp_ix))
                        else:
                            fo.write(line)


def set_up_sra_cache_folder(temp_folder):
    """Set up the fastq-dump cache folder within the temp folder."""
    logging.info("Setting up fastq-dump cache within {}".format(temp_folder))
    for path in [
        "/root/ncbi",
        "/root/ncbi/public"
    ]:
        if os.path.exists(path) is False:
            os.mkdir(path)

    if os.path.exists("/root/ncbi/public/sra"):
        shutil.rmtree("/root/ncbi/public/sra")

    # Now make a folder within the temp folder
    temp_cache = os.path.join(temp_folder, "sra")
    assert os.path.exists(temp_cache) is False
    os.mkdir(temp_cache)

    # Symlink it to /root/ncbi/public/sra/
    run_cmds(["ln", "-s", "-f", temp_cache, "/root/ncbi/public/sra"])

    assert os.path.exists("/root/ncbi/public/sra")


def get_reads_from_url(
    input_str,
    temp_folder,
    random_string=str(uuid.uuid4())[:8],
    min_qual=None
):
    """Get a set of reads from a URL -- return the downloaded filepath."""
    # Fetch reads into $temp_folder/fetched_reads/
    fetched_reads_folder = os.path.join(temp_folder, "fetched_reads")

    # Reads with cleaned headers go into $temp_folder/cleaned_reads/
    cleaned_reads_folder = os.path.join(temp_folder, "cleaned_reads")

    # Quality trimmed reads go into $temp_folder/trimmed_reads/
    trimmed_reads_folder = os.path.join(temp_folder, "trimmed_reads")

    for folder in [
        fetched_reads_folder, cleaned_reads_folder, trimmed_reads_folder
    ]:
        if not os.path.exists(folder):
            logging.info("Making new folder {}".format(folder))
            os.mkdir(folder)

    logging.info("Getting reads from {}".format(input_str))

    filename = input_str.split('/')[-1]
    local_path = os.path.join(fetched_reads_folder, filename)

    logging.info("Filename: " + filename)
    logging.info("Local path: " + local_path)

    if not input_str.startswith(('s3://', 'sra://', 'ftp://')):
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
            'AES256', input_str, fetched_reads_folder
            ])

    # Get files from an FTP server
    elif input_str.startswith('ftp://'):
        logging.info("Getting reads from FTP")
        run_cmds(['wget', '-P', fetched_reads_folder, input_str])

    # Get files from SRA
    elif input_str.startswith('sra://'):
        accession = filename
        logging.info("Getting reads from SRA: " + accession)
        local_path = get_sra(accession, fetched_reads_folder)

    else:
        raise Exception("Did not recognize prefix for input: " + input_str)

    # Clean up the FASTQ headers
    logging.info("Cleaning up FASTQ headers")
    cleaned_path = clean_fastq_headers(
        local_path,
        cleaned_reads_folder
    )
    logging.info("Made new cleaned FASTQ file: {}".format(cleaned_path))
    logging.info("Deleting old file: {}".format(local_path))
    os.unlink(local_path)

    if min_qual is None:
        return cleaned_path
    else:
        # Quality trim the FASTQ
        logging.info("Quality trimming the FASTQ (Q{})".format(min_qual))
        trimmed_path = quality_trim(
            cleaned_path,
            trimmed_reads_folder,
            min_qual
        )
        logging.info("Made new quality trimmed FASTQ: {}".format(trimmed_path))
        logging.info("Deleting old file: {}".format(cleaned_path))
        os.unlink(cleaned_path)
        return trimmed_path


def quality_trim(fp_in, folder_out, min_qual, min_len=30):
    """Trim a FASTQ to a minimum quality score."""
    assert os.path.exists(fp_in)
    assert fp_in.endswith(".gz") is False
    assert os.path.exists(folder_out)
    assert isinstance(min_qual, int)

    fp_out = os.path.join(folder_out, fp_in.split("/")[-1])

    run_cmds([
        "fastq_quality_trimmer",
        "-Q", "33",
        "-t", str(min_qual),
        "-i", fp_in,
        "-o", fp_out,
        "-l", str(min_len),
        "-v"
    ])

    assert os.path.exists(fp_out)

    return fp_out


def get_sra(accession, temp_folder):
    """Get the FASTQ for an SRA accession."""
    logging.info("Downloading {} from SRA".format(accession))

    local_path = os.path.join(temp_folder, accession + ".fastq")
    logging.info("Local path: {}".format(local_path))

    # Download via fastq-dump
    logging.info("Downloading via fastq-dump")
    run_cmds([
        "prefetch", accession
    ])
    run_cmds([
        "fastq-dump",
        "--split-files",
        "--outdir",
        temp_folder, accession
    ])

    # Make sure that some files were created
    msg = "File could not be downloaded from SRA: {}".format(accession)
    assert any([
        fp.startswith(accession) and fp.endswith("fastq")
        for fp in os.listdir(temp_folder)
    ]), msg

    # Combine any multiple files that were found
    logging.info("Concatenating output files")
    with open(local_path + ".temp", "wt") as fo:
        cmd = "cat {}/{}*fastq".format(temp_folder, accession)
        cat = subprocess.Popen(cmd, shell=True, stdout=fo)
        cat.wait()

    # Remove the temp files
    for fp in os.listdir(temp_folder):
        if fp.startswith(accession) and fp.endswith("fastq"):
            fp = os.path.join(temp_folder, fp)
            logging.info("Removing {}".format(fp))
            os.unlink(fp)

    # Remove the cache file, if any
    cache_fp = "/root/ncbi/public/sra/{}.sra".format(accession)
    if os.path.exists(cache_fp):
        logging.info("Removing {}".format(cache_fp))
        os.unlink(cache_fp)

    # Clean up the FASTQ headers for the downloaded file
    run_cmds(["mv", local_path + ".temp", local_path])

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


def clean_fastq_headers(fp_in, folder_out):
    """Read in a FASTQ file and write out a copy with unique headers."""

    # Constraints
    # 1. Headers start with '@'
    # 2. Headers are stripped to the first whitespace
    # 3. Headers are unique
    # 4. Sequence lines are not empty
    # 5. Spacer lines match the header line
    # 6. Quality lines are not empty

    # Make a new file in the output folder
    fp_out = os.path.join(folder_out, fp_in.split("/")[-1])
    # Don't gzip the output
    if fp_out.endswith(".gz"):
        fp_out = fp_out[:-3]

    if fp_in.endswith(".gz"):
        f_in = gzip.open(fp_in, "rt")
    else:
        f_in = open(fp_in, "rt")

    f_out = open(fp_out, "wt")

    # Compile a regex to mask non ATCG
    atcg = re.compile('[^ATCG\n]')

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

            # Replace any non-ATCG with N
            line = atcg.sub('N', line)

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

    # Return the path to the file that was written
    return fp_out
