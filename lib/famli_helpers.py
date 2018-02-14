#!/usr/bin/python

import os
import gzip
import json
import logging
import argparse
from .exec_helpers import run_cmds

from collections import defaultdict
import numpy as np


def align_reads(read_fp,               # FASTQ file path
                db_fp,                 # Local path to DB
                temp_folder,           # Folder for results
                query_gencode=11,      # Genetic code
                threads=1,             # Threads
                min_score=20,          # Minimum alignment score
                blocks=4):             # Memory block size

    """Align a set of reads with DIAMOND."""

    align_fp = "{}.sam".format(read_fp)
    logging.info("Input reads: {}".format(read_fp))
    logging.info("Reference database: {}".format(db_fp))
    logging.info("Genetic code: {}".format(query_gencode))
    logging.info("Threads: {}".format(threads))
    logging.info("Output: {}".format(align_fp))

    run_cmds([
            "diamond",
            "blastx",
            "--query", read_fp,             # Input FASTQ
            "--out", align_fp,              # Alignment file
            "--threads", str(threads),      # Threads
            "--db", db_fp,                  # Reference database
            "--outfmt", "6",                # Output format
            "qseqid", "sseqid",
            "pident", "length",
            "mismatch", "gapopen",
            "qstart", "qend",
            "sstart", "send",
            "evalue", "bitscore",
            "qlen", "slen",
            "--min-score", str(min_score),  # Minimum alignment score
            "--query-cover", "50",          # Minimum query coverage
            "--id", "80",                   # Minimum alignment identity
            "--max-target-seqs", "0",       # Report all alignments
            "--block-size", str(blocks),    # Memory block size
            "--query-gencode",              # Genetic code
            str(query_gencode),
            "--unal", "0",                  # Don't report unaligned reads
            ])

    return align_fp


class BLAST6Parser:
    """Object to help with parsing alignments in BLAST 6 format."""

    def __init__(self):
        # Keep track of the length of all subjects
        self.subject_len = {}

        # Set of all possible query_ids
        self.query_set = set()

    def parse(self,
              align_handle,
              QSEQID_i=0,
              SSEQID_i=1,
              SSTART_i=8,
              SEND_i=9,
              BITSCORE_i=11,
              SLEN_i=13):
        """Parse a file, while keeping track of the subject lengths."""
        for i, line in enumerate(align_handle):
            if i % 1000000 == 0 and i > 0:
                logging.info("{} lines of alignment parsed".format(i))
            line_list = line.strip().split()
            # Get the query and subject
            query = line_list[QSEQID_i]
            subject = line_list[SSEQID_i]

            # Save information for the query and subject
            self.query_set.add(query)
            self.subject_len[subject] = int(line_list[SLEN_i])

            # Yield a tuple with query, subject, sstart, send, bitscore
            yield (
                query,
                subject,
                int(line_list[SSTART_i]),
                int(line_list[SEND_i]),
                float(line_list[BITSCORE_i])
            )


def parse_alignment(align_handle,
                    QSEQID_i=0,
                    SSEQID_i=1,
                    QSTART_i=6,
                    QEND_i=7,
                    SSTART_i=8,
                    SEND_i=9,
                    BITSCORE_i=11,
                    SLEN_i=13,
                    SD_MEAN_CUTOFF=0.33,
                    STRIM_5=18,
                    STRIM_3=18,
                    ITERATIONS_MAX=1000):
    """
    Parse an alignment in BLAST6 format and determine which subjects are likely to be present. This is the core of FAMLI.
    BLAST 6 columns by default (in order): qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
                                              0     1       2       3       4       5       6   7   8       9   10      11      12   13  
    """

    # 1. Look at the how well each subject could be covered and to get a sense of the subjects and queries present

    # Fill our data structures by parsing our alignment.
    logging.info("Starting to parse the alignment for the first time.")
    parser = BLAST6Parser()
    alignments = [
        a
        for i, a in enumerate(parser.parse(
            align_handle,
            QSEQID_i=QSEQID_i,
            SSEQID_i=SSEQID_i,
            SSTART_i=SSTART_i,
            SEND_i=SEND_i,
            BITSCORE_i=BITSCORE_i,
            SLEN_i=SLEN_i))
    ]

    # After parsing, the subject lengths are in parser.subject_len
    # and the query set is in parser.query_set

    logging.info("Done parsing {} subjects and {} queries".format(
        len(parser.subject_len),
        len(parser.query_set)
    ))

    # Now we create coverage-o-grams for each subject,
    # considering how it COULD be covered
    subject_could_coverage = {}
    # Iterate over the alignments
    for query, subject, sstart, send, bitscore in alignments:
        # For new subjects, make an empty vector
        if subject not in subject_could_coverage:
            subject_could_coverage[subject] = np.zeros(
                shape=parser.subject_len[subject],
                dtype=int
                )
        # Increment the coverage vector
        subject_could_coverage[subject][sstart:send] += 1

    # Now go through those coverage vectors and reject those for which the
    # coverage is uneven enough to be implausible
    subjects_below_cutoff = set()
    # Iterate over all subjects
    for subject in parser.subject_len:
        # Trim the 5' and 3' ends of the subject.
        # Due to overhangs, they tend to be less covered.
        subj_cov_trimmed = subject_could_coverage[subject][STRIM_5:-STRIM_3]
        # Using the trimmed coverage-o-gram, if the STD / mean of coverage is
        # LESS than a threshold (empirically set to 0.33), we keep it
        subject_std = np.std(subj_cov_trimmed)
        subject_mean = np.mean(subj_cov_trimmed)
        if (subject_std / subject_mean) <= SD_MEAN_CUTOFF:
            # Add to the set of subjects passing the filter
            subjects_below_cutoff.add(subject)

    logging.info("Kept {:,} of {:,} total subjects after filtering for evenness of possible coverage".format(len(subjects_below_cutoff),len(parser.subject_len)))

    # Reduce the alignments to just those subjects passing the filter
    alignments = [a for a in alignments if a[1] in subjects_below_cutoff]

    logging.info(
        "Number of alignments passing the first evenness filter: {:,}".format(
            len(alignments)
        ))

    # We will eventually report output for every subject, including:
    # id: name of the subject
    # nreads: number of reads
    # coverage: proportion of the subject that is covered
    # depth: average depth across the subject
    # length: length of the subject
    output = []

    # 2. Find the groups of subjects that share any queries
    for subjects, queries, group_alignments in group_cooccurring_subjects(
        alignments
    ):

        # 3. Use the combination of the bitscores (alignment quality) and
        # subject-read-depth to iterative filter low-likely alignments of
        # queries against subjects

        # Indicies for quick lookup
        subject_i = {n: i for i, n in enumerate(subjects)}
        query_i = {n: i for i, n in enumerate(queries)}

        # Subject Length vector to use for later normalization of the subject read depths
        subject_lengths = np.zeros(shape=(len(subject_i)), dtype=int)
        for subject in subject_i:
            subject_lengths[subject_i[subject]] = parser.subject_len[subject]

        # Numpy matrix to store the bitscores
        bitscore_mat = np.zeros(
            shape=(len(subject_i), len(query_i)),
            dtype='float64'
        )

        logging.info("Fill in the bitscore matrix")
        for query, subject, sstart, send, bitscore in group_alignments:
            bitscore_mat[subject_i[subject], query_i[query]] = bitscore

        logging.info("Completed filling the bitscore matrix")

        # Iterative alignment filtering
        # Idea here is to take the alignemnt scores (bitscores here) plus the length-normalized subject coverage to calcualate a likelihood [0-1] that a given query came from a given subject and then prune the lowest likelihood alignments.
        # We repeat this process iteratively (reweighting after pruning) to until we either converge OR reach our maximum iterations. 
        iterations_n = 0

        # Initialize some variables
        prior_align_norm_min_max = -1

        while iterations_n < ITERATIONS_MAX:
            iterations_n += 1
            # Weight each alignment score by the total mass of the subject
            align_mat_w = (
                bitscore_mat.T * (np.sum(bitscore_mat, axis=1)/subject_lengths)
            ).T

            # And then normalize so each query sums to 1.0
            align_mat_w = (
                align_mat_w.astype('float') / np.sum(align_mat_w, axis=0)
            )

            # Generate a per-query vector of maximum normalized [0-1] alignment scores 
            # Then take the minimum of that vector.
            align_norm_min_max = np.min(align_mat_w.max(axis=0))

            # Zero out any alignments below our minumum maximum alignment score
            bitscore_mat[align_mat_w < align_norm_min_max] = 0.0

            logging.info("Iteration {}. Min {:.4f} Mean {:.3f}".format(iterations_n,align_mat_w[align_mat_w > 0.0].min(), align_mat_w[align_mat_w > 0.0].mean()))
            if prior_align_norm_min_max == align_norm_min_max:
                print("Interations complete")
                break
            # Implicit else
            prior_align_norm_min_max = align_norm_min_max
            # and proceed to next iteration

        # Now let us get the subjects with some remaining alignments.
        subjects_with_iteratively_aligned_reads = {[p[0] for p in subject_i.items() if p[1] == i][0]: i for i in np.argwhere(align_mat_w.sum(axis=1)).flatten()}

        # 3. Regenerate coverage-o-grams for the subjects that still have iteratively mapped queries. 
        # Create coverage-O-grams for the subjects with iteratively assigned reads
        # Use our same SD / mean coverage metric to screen now with the reassigned reads.
        logging.info("Starting final coverage evenness screening of {} subjects with filtered alignments.".format(len(subjects_with_iteratively_aligned_reads)))
        subject_final_passed = set()

        for subject in subjects_with_iteratively_aligned_reads:
            # np zero vector of ints to be our coverage-o-gram
            subject_coverage = np.zeros(shape=(parser.subject_len[subject],), dtype=int)

            # For subject s, which query_ids still have a non-zero normalized alignment score. 
            subject_queries = {[p[0] for p in query_i.items() if p[1] == i][0] for i in np.argwhere(align_mat_w[subjects_with_iteratively_aligned_reads[subject]]).flatten()}
            for q, s, sstart, send, bitscore in group_alignments:
                # Fill in the coverage-o-gram with these queries
                if s == subject:
                    subject_coverage[sstart:send] += 1

            # Our testing point. With the trimmed coverage-o-gram are we still below our evenness threshold?
            subject_std = np.std(subject_coverage[STRIM_5:-STRIM_3])
            subject_mean = np.mean(subject_coverage[STRIM_5:-STRIM_3])
            if subject_std / subject_mean < SD_MEAN_CUTOFF:
                output.append({
                    "id": subject,
                    "length": parser.subject_len[subject],
                    "depth": np.mean(subject_coverage),
                    "coverage": np.mean(subject_coverage > 0),
                    "nreads": len(subject_queries),
                })

        logging.info("Done filtering subjects. {} subjects are felt to likely be present".format(len(subject_final_passed)))

    return output


def group_cooccurring_subjects(alignments):
    """Get the sets of subjects that share any queries."""

    # Get the dict of the subjects that each query aligns to
    query_subject_dict = defaultdict(set)
    for query, subject, sstart, send, bitscore in alignments:
        query_subject_dict[query].add(subject)

    # Now assemble the groups of subjects that share any queries
    subject_groups = []
    for query_ix, g in enumerate(query_subject_dict.values()):
        # See which pre-existing groups this query matches
        matching_groups = [
            ix
            for ix, q in enumerate(subject_groups)
            if g & q
        ]
        if len(matching_groups) == 0:
            # No existing groups match
            # Add to the end of the list
            subject_groups.append(g)
        elif len(matching_groups) == 1:
            # A single existing group matches
            # Add to the existing group
            subject_groups[matching_groups[0]] |= g
        else:
            # Multiple groups match
            # Make a new larger group
            for ix in matching_groups:
                g |= subject_groups[ix]
            # Remove the existing matches
            subject_groups = [
                q
                for ix, q in enumerate(subject_groups)
                if ix not in matching_groups
            ]
            # Add to the end of the list
            subject_groups.append(g)
        if query_ix % 1e6 == 0 and query_ix > 0:
            msg = "Processing {:,} queries to make {:,} co-occuring subject groups"
            logging.info(
                msg.format(query_ix, len(subject_groups)))

    # Return the results
    msg = "Processed {:,} queries to make {:,} co-occuring subject groups"
    logging.info(
        msg.format(query_ix, len(subject_groups)))
    for g in subject_groups:
        # Get the alignments matching this group
        group_alignments = [a for a in alignments if a[1] in g]
        group_queries = set([
            query
            for query, subject, sstart, send, bitscore in group_alignments
        ])
        yield g, group_queries, group_alignments


def parse_alignments_by_query(align_fp):
    """Function to parse the SAM alignments, grouping them by query ID."""
    # Column names
    header = [
        "query", "subject", "pct_iden",
        "alen", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send",
        "evalue", "bitscore", "qlen", "slen",
    ]
    # Keep all of the alignments for a given read here
    buff = []
    # Keep track of the query ID for the previous alignment
    last_query = None
    with open(align_fp, "rt") as f:
        for line in f:
            line = line.rstrip("\n").split("\t")
            if len(line) != len(header):
                continue

            # Get the query ID
            query = line[0]

            # If this is a new query ID
            if query != last_query:
                # If there is data for the last query, yield it
                if len(buff) > 0:
                    yield buff
                    # Reset the buffer
                    buff = []
                # Reset the last query sequence
                last_query = query

            # Add this alignment's information to the buffer
            buff.append(dict(zip(header, line)))

    # If there is data for the last query, yield it
    if len(buff) > 0:
        yield buff


def parse_alignment_score(tags):
    """Parse the alignment score from a set of SAM tags."""
    for t in tags:
        if t[:2] == "AS":
            return float(t[5:])


def parse_mismatch_tag(tags):
    """Parse the alignment score from a set of SAM tags."""
    for t in tags:
        if t[:2] == "NM":
            return int(t[5:])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Parse a DIAMOND output file and save the deduplicated coverage data.
    """)
    parser.add_argument("--input",
                        type=str,
                        help="""DIAMOND output file in tabular format.""")
    parser.add_argument("--output",
                        type=str,
                        help="""Output file in JSON format.""")

    args = parser.parse_args()

    assert os.path.exists(args.input)

    # Set up logging
    logFormatter = logging.Formatter(
        '%(asctime)s %(levelname)-8s [FAMLI Parser] %(message)s'
    )
    rootLogger = logging.getLogger()
    rootLogger.setLevel(logging.INFO)

    # Write to STDOUT
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)

    if args.input.endswith(".gz"):
        with gzip.open(args.input, "rt") as f:
            total_reads, deduplicated_reads, output = parse_alignment(f)
    else:
        with open(args.input, "rt") as f:
            total_reads, deduplicated_reads, output = parse_alignment(f)

    with open(args.output, "wt") as fo:
        json.dump(output, fo)
