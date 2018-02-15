#!/usr/bin/python

import os
import gzip
import json
import logging
import argparse
import numpy as np
from multiprocessing import Pool


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
                logging.info("{:,} lines of alignment parsed".format(i))
            line_list = line.strip().split()
            # Get the query and subject
            query = line_list[QSEQID_i]
            subject = line_list[SSEQID_i]

            # Save information for the query and subject
            self.query_set.add(query)
            self.subject_len[subject] = int(line_list[SLEN_i])

            # Yield a tuple with query, subject, sstart, send, bitscore
            sstart = int(line_list[SSTART_i])
            send = int(line_list[SEND_i])

            assert send > sstart

            yield (
                query,
                subject,
                sstart - 1,  # Convert to 0-based indexing
                send,
                float(line_list[BITSCORE_i])
            )


def coverage_filter(args):
    """Test whether a given subject passes the coverage filter."""
    subject, subject_len, subject_alignments, STRIM_5, STRIM_3, SD_MEAN_CUTOFF = args
    subject_could_coverage = np.zeros(
        shape=subject_len,
        dtype=int
    )
    # Add the alignments to the coverage
    for query, subject, sstart, send, bitscore in subject_alignments:
        # Increment the coverage vector
        subject_could_coverage[sstart:send] += 1

    # Reject this subject if coverage is uneven enough to be implausible

    # Trim the 5' and 3' ends of the subject.
    # Due to overhangs, they tend to be less covered.
    # Only trim subjects that have >= 30aa after trimming
    if subject_len - (STRIM_5 + STRIM_3) >= 30:
        subject_could_coverage = subject_could_coverage[STRIM_5:-STRIM_3]

    # Using the trimmed coverage-o-gram, if the STD / mean of coverage is
    # LESS than a threshold (empirically set to 0.33), we keep it
    # THIS IS THE MOST COMPUTATIONALLY COSTLY PART OF THIS LOOP
    subject_std = np.std(subject_could_coverage)
    subject_mean = np.mean(subject_could_coverage)
    passes_filter = subject_mean > 0 and (subject_std / subject_mean) <= SD_MEAN_CUTOFF
    return subject, passes_filter


def parse_alignment(align_handle,
                    QSEQID_i=0,
                    SSEQID_i=1,
                    QSTART_i=6,
                    QEND_i=7,
                    SSTART_i=8,
                    SEND_i=9,
                    BITSCORE_i=11,
                    SLEN_i=13,
                    SD_MEAN_CUTOFF=1.0,
                    STRIM_5=18,
                    STRIM_3=18,
                    ITERATIONS_MAX=1000,
                    threads=4):
    """
    Parse an alignment in BLAST6 format and determine which subjects are likely to be present. This is the core of FAMLI.
    BLAST 6 columns by default (in order): qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
                                              0     1       2       3       4       5       6   7   8       9   10      11      12   13  
    """

    # ######################
    # 1. Parse the alignment
    # ######################

    # Fill our data structures by parsing our alignment.
    logging.info("Starting to parse the alignment.")
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

    logging.info("Done parsing {:,} alignments".format(len(alignments)))
    logging.info("Subjects: {:,}".format(
        len(parser.subject_len),
    ))
    logging.info("Queries: {:,}".format(
        len(parser.query_set)
    ))

    # Sort the alignments by subject
    logging.info("Sorting alignments by subject")
    alignments = sorted(alignments, key=lambda a: a[1])
    logging.info("Done sorting alignments")

    # Collect the index range for each subject
    logging.info("Indexing alignments by subject")
    subject_index_ranges = {}
    last_subject = None
    last_start_ix = None
    # Iterate through the alignments
    for ix, a in enumerate(alignments):
        # Check to see if this is a new subject
        if a[1] != last_subject:
            # Add the subject to the dict of index ranges
            if last_subject is not None:
                # Make sure the subject isn't being added twice
                # (would indicate the list wasn't sorted completely)
                assert last_subject not in subject_index_ranges
                # The value is the start and end index positions
                subject_index_ranges[last_subject] = (last_start_ix, ix)
            # Set the start index for the new subject
            last_start_ix = ix
            # Set the name for the new subject
            last_subject = a[1]
    # Add the final subject
    subject_index_ranges[last_subject] = (last_start_ix, ix + 1)
    # Make sure that we have ranges for all of the subjects
    assert len(subject_index_ranges) == len(parser.subject_len)
    logging.info("Done indexing by subject")

    # #############################
    # 2. Calculate coverage metrics
    # #############################

    # Now we create coverage-o-grams for each subject,
    # considering how it COULD be covered
    # Keep a list of the subjects that will be KEPT
    subjects_below_cutoff = set()
    logging.info("Calculating coverage metrics for {:,} subjects".format(
        len(subject_index_ranges)))
    # Iterate over the subjects
    # THIS IS PRIME FOR MULTITHREADING

    p = Pool(threads)
    for subject, passes_filter in p.map(coverage_filter, [
        [
            subject,
            parser.subject_len[subject],
            alignments[
                align_ix[0]: align_ix[1]
            ],
            STRIM_5,
            STRIM_3,
            SD_MEAN_CUTOFF,
        ]
        for subject, align_ix in subject_index_ranges.items()
    ]):
        if passes_filter:
            # Add to the set of subjects passing the filter
            subjects_below_cutoff.add(subject)
        else:
            # Remove from the list of index ranges
            del subject_index_ranges[subject]

    logging.info("Kept {:,} of {:,} total subjects after filtering for evenness of possible coverage".format(len(subjects_below_cutoff),len(parser.subject_len)))

    # We will eventually report output for every subject, including:
    # id: name of the subject
    # nreads: number of reads
    # coverage: proportion of the subject that is covered
    # depth: average depth across the subject
    # length: length of the subject
    output = []

    # 2. Find the groups of subjects that share any queries
    for subjects, queries, group_alignments in group_cooccurring_subjects(
        alignments, subject_index_ranges
    ):

        # TODO: Only perform more filtering if there is > 1 reference

        logging.info("")
        logging.info("Assigning {:,} reads among {:,} references".format(
            len(queries), len(subjects)
            ))

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

            logging.info("Iteration {:,}. Min {:.4f} Mean {:.3f}".format(iterations_n,align_mat_w[align_mat_w > 0.0].min(), align_mat_w[align_mat_w > 0.0].mean()))
            if prior_align_norm_min_max == align_norm_min_max:
                logging.info("Interations complete")
                break
            # Implicit else
            prior_align_norm_min_max = align_norm_min_max
            # and proceed to next iteration

        # Now let us get the subjects with some remaining alignments.
        subjects_with_iteratively_aligned_reads = {[p[0] for p in subject_i.items() if p[1] == i][0]: i for i in np.argwhere(align_mat_w.sum(axis=1)).flatten()}

        # 3. Regenerate coverage-o-grams for the subjects that still have iteratively mapped queries. 
        # Create coverage-O-grams for the subjects with iteratively assigned reads
        # Use our same SD / mean coverage metric to screen now with the reassigned reads.
        logging.info("Starting coverage evenness screening for {:,} subjects with filtered alignments.".format(len(subjects_with_iteratively_aligned_reads)))
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
                subject_final_passed.add(subject)

        logging.info("Done filtering subjects. {:,} subjects are felt to likely be present".format(len(subject_final_passed)))

    logging.info("Number of input queries: {:,}".format(len(parser.query_set)))
    logging.info("Number of input subjects: {:,}".format(len(parser.subject_len)))
    logging.info("Number of input alignments: {:,}".format(len(alignments)))
    logging.info("Number of deduplicated queries: {:,}".format(
        sum([d["nreads"] for d in output])
    ))
    logging.info("Number of final subjects: {:,}".format(
        len(output)
    ))

    return len(parser.query_set), output


def group_cooccurring_subjects(alignments, subject_index_ranges):
    """Get the sets of subjects that share any queries."""

    # Get the dict of the queries that each subject aligns to
    subject_queries = {
        subject: set([
            a[0]
            for a in alignments[index_range[0]: index_range[1]]
        ])
        for subject, index_range in subject_index_ranges.items()
    }

    # Make a list of all the subjects
    all_subjects = [s for s in subject_queries]

    # Make a list for the group number for those subjects
    subject_groups = []

    # Iterate through the list of subjects
    for new_ix, new_subject in enumerate(all_subjects):
        if new_ix % 100 == 0 and new_ix > 0:
            logging.info("Grouping {:,} subjects by queries".format(new_ix))
        # Check to see which of the previous subjects have overlap
        for old_ix, old_subject in enumerate(all_subjects[:new_ix]):
            # There is overlap
            if subject_queries[old_subject] & subject_queries[new_subject]:
                # Assign that previous subject to the new group
                subject_groups[old_ix] = new_ix
        # Add the new subject to this new group
        subject_groups.append(new_ix)

    # Return the results
    msg = "Processed {:,} subjects to make {:,} co-occuring subject groups"
    logging.info(
        msg.format(len(subject_groups), len(set(subject_groups))))

    # Yield each of the subject groups
    for group in list(set(subject_groups)):
        # Get the subjects in this group
        subjects = [
            subject
            for subject, g in zip(all_subjects, subject_groups)
            if group == g
        ]
        # Get the alignments matching this group
        group_alignments = [
            a
            for s in subjects
            for a in alignments[
                subject_index_ranges[s][0]:
                subject_index_ranges[s][1]
            ]
        ]
        group_queries = set([])
        for s in subjects:
            group_queries |= subject_queries[s]

        yield subjects, group_queries, group_alignments


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
    parser.add_argument("--threads",
                        type=int,
                        help="""Number of processors available.""")

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
            aligned_reads, output = parse_alignment(f, threads=args.threads)
    else:
        with open(args.input, "rt") as f:
            aligned_reads, output = parse_alignment(f, threads=args.threads)

    with open(args.output, "wt") as fo:
        json.dump(output, fo)
