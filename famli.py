#!/usr/bin/python
"""Wrapper script to run FAMLI on one or more FASTQ files."""

import os
import uuid
import time
import shutil
import logging
import argparse
from famli.exec_helpers import align_reads
from famli.exec_helpers import return_results
from famli.exec_helpers import exit_and_clean_up
from famli.exec_helpers import get_reference_database
from famli.fastq_helpers import get_reads_from_url
from famli.fastq_helpers import count_fastq_reads
from famli.fastq_helpers import combine_fastqs
from famli.famli_helpers import parse_alignment


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Align a set of reads against a reference database with DIAMOND,
    process those alignments with FAMLI, and save the results.
    """)

    parser.add_argument("--input",
                        type=str,
                        help="""Location for input file(s). Combine multiple files with +.
                                (Supported: sra://, s3://, or ftp://).""")
    parser.add_argument("--sample-name",
                        type=str,
                        help="""Name of sample, sets output filename.""")
    parser.add_argument("--ref-db",
                        type=str,
                        help="""Folder containing reference database.
                                (Supported: s3://, ftp://, or local path).""")
    parser.add_argument("--output-folder",
                        type=str,
                        help="""Folder to place results.
                                (Supported: s3://, or local path).""")
    parser.add_argument("--min-score",
                        type=float,
                        default=20,
                        help="Minimum alignment score to report.")
    parser.add_argument("--blocks",
                        type=int,
                        default=5,
                        help="""Number of blocks used when aligning.
                              Value relates to the amount of memory used.""")
    parser.add_argument("--query-gencode",
                        type=int,
                        default=11,
                        help="Genetic code used to translate nucleotides.")
    parser.add_argument("--threads",
                        type=int,
                        default=16,
                        help="Number of threads to use aligning.")
    parser.add_argument("--min-qual",
                        type=int,
                        default=None,
                        help="If specified, trim reads to a minimum Q score.")
    parser.add_argument("--temp-folder",
                        type=str,
                        default='/share',
                        help="Folder used for temporary files.")

    args = parser.parse_args()

    # Make sure that there are no commas or whitespaces in the input
    input_str = args.input
    assert ' ' not in input_str, input_str
    assert ',' not in input_str, input_str

    # Make a temporary folder for all files to be placed in
    temp_folder = os.path.join(args.temp_folder, str(uuid.uuid4())[:8])
    assert os.path.exists(temp_folder) is False
    os.mkdir(temp_folder)

    # Set up logging
    log_fp = os.path.join(temp_folder, "log.txt")
    logFormatter = logging.Formatter(
        '%(asctime)s %(levelname)-8s [FAMLI] %(message)s'
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

    # Get the reference database
    try:
        db_fp = get_reference_database(
            args.ref_db,
            temp_folder
        )
    except:
        exit_and_clean_up(temp_folder)

    logging.info("Reference database: " + db_fp)

    # Align the input data and calculate the overall abundance

    # Keep track of the time elapsed to process this sample
    start_time = time.time()

    logging.info("Processing input argument: " + input_str)

    # Multiple input reads may be separated with a '+'
    input_str = input_str.split("+")
    # Make sure that they are all unique arguments
    assert len(input_str) == len(set(input_str)), "Duplicate arguments"
    # Make sure that the filenames are also all unique
    assert len(input_str) == len(set([
        s.split('/')[-1] for s in input_str
    ])), "Duplicate filenames"

    # Capture each command in a try statement
    # Get the input reads
    read_fps = []
    for s in input_str:
        logging.info("Fetching {}".format(s))
        try:
            read_fps.append(get_reads_from_url(
                s, temp_folder, min_qual=args.min_qual))
        except:
            exit_and_clean_up(temp_folder)

    # Combine the files into a single FASTQ
    read_fp = os.path.join(temp_folder, "input.fastq")
    combine_fastqs(read_fps, read_fp)

    # Run the alignment
    try:
        align_fp = align_reads(
            read_fp,               # FASTQ file path
            db_fp,                 # Local path to DB
            temp_folder,    # Folder for results
            query_gencode=args.query_gencode,
            threads=args.threads,
            min_score=args.min_score,
            blocks=args.blocks,
        )
    except:
        exit_and_clean_up(temp_folder)

    # Process the alignments, reassigning multi-mapped reads
    try:
        with open(align_fp, "rt") as align_handle:
            aligned_reads, abund = parse_alignment(
                align_handle,
            )
    except:
        exit_and_clean_up(temp_folder)

    # Calculate the number of deduplicated reads
    deduplicated_reads = sum([d["nreads"] for d in abund])

    # Name the output file based on the input file
    # Ultimately adding ".json.gz" to the input file name
    if args.sample_name is not None:
        output_prefix = args.sample_name
    else:
        output_prefix = input_str[0].split("/")[-1]
    logging.info("Using sample name {} for output prefix".format(
        output_prefix))

    # Count the total number of reads
    logging.info("Counting the total number of reads")
    n_reads = count_fastq_reads(read_fp)
    logging.info("Reads in input file: {}".format(n_reads))

    # Read in the logs
    logging.info("Reading in the logs")
    logs = open(log_fp, 'rt').readlines()

    # Wrap up all of the results into a single JSON
    # and write it to the output folder
    output = {
        "input_path": "+".join(input_str),
        "input": output_prefix,
        "output_folder": args.output_folder,
        "logs": logs,
        "ref_db": db_fp,
        "ref_db_url": args.ref_db,
        "results": abund,
        "total_reads": n_reads,
        "aligned_reads": aligned_reads,
        "deduplicated_reads": deduplicated_reads,
        "time_elapsed": time.time() - start_time
    }
    return_results(
        output, output_prefix, args.output_folder, temp_folder
    )

    # Delete any files that were created for this sample
    logging.info("Removing temporary folder: " + temp_folder)
    shutil.rmtree(temp_folder)

    # Stop logging
    logging.info("Done")
    logging.shutdown()
