#!/usr/bin/python

import os
import gzip
import json
import time
import logging
import argparse
import numpy as np
from multiprocessing import Pool
from collections import defaultdict


class BLAST6Parser:
    """Object to help with parsing alignments in BLAST 6 format."""

    def __init__(self):
        # Keep track of the length of all subjects
        self.subject_len = {}
        # Keep track of the number of queries
        self.unique_queries = set([])

    def parse(self,
              align_handle,
              QSEQID_i=0,
              SSEQID_i=1,
              SSTART_i=8,
              SEND_i=9,
              EVALUE_i=10,
              BITSCORE_i=11,
              SLEN_i=13):
        """Parse a file, while keeping track of the subject lengths."""
        logging.info("Reading in alignments")

        for i, line in enumerate(align_handle):
            if i % 1000000 == 0 and i > 0:
                logging.info("{:,} lines of alignment parsed".format(i))
            line_list = line.rstrip("\n").split("\t")
            # Get the query and subject
            query = line_list[QSEQID_i]
            subject = line_list[SSEQID_i]
            bitscore = float(line_list[BITSCORE_i])
            sstart = int(line_list[SSTART_i])
            send = int(line_list[SEND_i])

            # Save information for the subject length
            if subject not in self.subject_len:
                self.subject_len[subject] = int(line_list[SLEN_i])

            # Add the query to the set of unique queries
            self.unique_queries.add(query)

            # Yield a tuple with query, subject, sstart, send, bitscore
            yield (
                query,
                subject,
                sstart - 1,  # Convert to 0-based indexing
                send,
                bitscore
            )

        if len(self.unique_queries) == 0:
            logging.info("Zero alignments were found")
        else:
            logging.info("Done reading in {:,} alignments".format(i + 1))
            logging.info("Number of unique subjects: {:,}".format(
                len(self.subject_len)))
            logging.info("Number of unique queries: {:,}".format(
                len(self.unique_queries)))

    def yield_alignments(
        self,
        align_handle,
        batchsize=None,
        QSEQID_i=0,
        SSEQID_i=1,
        SSTART_i=8,
        SEND_i=9,
        EVALUE_i=10,
        BITSCORE_i=11,
        SLEN_i=13
    ):
        """Yield batches of alignments, never to splitting up queries."""

        # If the batchsize function has been disabled, return the full list
        if batchsize is None:
            logging.info("Batchsize feature disabled, reading all alignments")
            yield [
                a
                for a in self.parse(
                    align_handle,
                    QSEQID_i=QSEQID_i,
                    SSEQID_i=SSEQID_i,
                    SSTART_i=SSTART_i,
                    SEND_i=SEND_i,
                    EVALUE_i=EVALUE_i,
                    BITSCORE_i=BITSCORE_i,
                    SLEN_i=SLEN_i
                )
            ]
            return

        logging.info("Reading in batches of {:,} alignments".format(batchsize))

        assert batchsize > 0, "Batches must be larger than 0"

        # Allocate a list that is 10% larger than the batch size
        alignments = [None] * int(batchsize * 1.1)
        # Keep track of how many alignments have been added to the list
        n_alignments_filled = 0
        # Keep track of the number of unique queries in this batch
        n_unique_queries = 0
        # Keep track of which batch of alignments we are on
        batch_num = 1
        # Keep track of the last unique query ID
        last_query = None

        # Add to the alignments, stopping when we get to the batch size limit
        for a in self.parse(
            align_handle,
            QSEQID_i=QSEQID_i,
            SSEQID_i=SSEQID_i,
            SSTART_i=SSTART_i,
            SEND_i=SEND_i,
            EVALUE_i=EVALUE_i,
            BITSCORE_i=BITSCORE_i,
            SLEN_i=SLEN_i
        ):

            # When we reach the batch size and new query, yield the alignments
            if a[0] != last_query and n_alignments_filled >= batchsize:
                logging.info("Processing batch {:,}".format(batch_num))
                logging.info("Alignments: {:,}".format(n_alignments_filled))
                logging.info("Unique queries: {:,}".format(n_unique_queries))

                if n_alignments_filled < int(batchsize * 1.1):
                    yield alignments[:n_alignments_filled]
                else:
                    yield alignments
                # Reset the alignment data
                alignments = [None] * int(batchsize * 1.1)
                n_alignments_filled = 0
                n_unique_queries = 0
                batch_num += 1

            # Add to the alignment list
            if n_alignments_filled < int(batchsize * 1.1):
                alignments[n_alignments_filled] = a
            else:
                alignments.append(a)

            # Increment the counter of filled alignments
            n_alignments_filled += 1
            # Increment the counter of unique queries
            if a[0] != last_query:
                n_unique_queries += 1
            # Record which query we saw in this alignment
            last_query = a[0]

        # Yield the final block of alignments
        logging.info("Processing batch {:,}".format(batch_num))
        logging.info("Alignments: {:,}".format(n_alignments_filled))
        logging.info("Unique queries {:,}".format(n_unique_queries))
        if n_alignments_filled < int(batchsize * 1.1):
            yield alignments[:n_alignments_filled]
        else:
            yield alignments


def recalc_subject_weight_worker(args):
    """Calculate the new weight of a subject."""
    query_weights, subject_len, subject_name = args
    return subject_name, sum(query_weights) / subject_len


class FAMLI_Reassignment:
    """Object to help with parsing alignments in BLAST 6 format."""

    def __init__(self, alignments, subject_len, pool=None):
        # Save the length of all subjects
        self.subject_len = subject_len

        # Keep track of the bitscores for each query/subject
        self.bitscores = defaultdict(lambda: defaultdict(float))

        # Keep track of the number of uniquely aligned queries
        self.n_unique = 0

        # Pool of workers
        if pool is None:
            self.pool = Pool(4)
        else:
            self.pool = pool

        logging.info("Adding alignments to read re-assignment model")
        for query, subject, sstart, send, bitscore in alignments:
            # Save information for the query and subject
            self.bitscores[query][subject] = bitscore

    def init_subject_weight(self):
        """Initialize the subject weights, equal and normalized to their length."""
        logging.info("Initializing subject weights")
        self.subject_weight = {
            subject: 1.0 / length
            for subject, length in self.subject_len.items()
        }

        # Also initialize the alignment probabilities as the bitscores
        logging.info("Initializing alignment probabilities")
        self.aln_prob = defaultdict(dict)
        self.aln_prob_T = defaultdict(dict)

        # Keep track of the queries and subjects that will need to be updated
        self.subjects_to_update = set([])
        self.queries_to_update = set(self.bitscores.keys())

        # Keep track of the queryies that have multiple alignments
        self.multimapped_queries = set(self.bitscores.keys())

        for query, bitscores in self.bitscores.items():
            if len(bitscores) == 1:
                self.n_unique += 1
            bitscore_sum = sum(bitscores.values())
            for subject, bitscore in bitscores.items():
                self.subjects_to_update.add(subject)
                v = bitscore / bitscore_sum
                self.aln_prob[query][subject] = v
                self.aln_prob_T[subject][query] = v

    def recalc_subject_weight(self):
        """Recalculate the subject weights."""
        self.queries_to_update = set([])

        self.queries_to_update = set([
            query
            for subject in self.subjects_to_update
            for query in self.aln_prob_T[subject].keys()
        ])

        for subject, new_weight in self.pool.imap(
            recalc_subject_weight_worker,
            [
                [
                    list(self.aln_prob_T[subject].values()),
                    self.subject_len[subject],
                    subject
                ]
                for subject in self.subjects_to_update
            ]
        ):
            self.subject_weight[subject] = new_weight

    def recalc_aln_prob(self):
        """Recalculate the alignment probabilities."""
        self.subjects_to_update = set([])

        # Iterate over every query
        for query in self.queries_to_update:
            aln_prob = self.aln_prob[query]
            new_probs = [
                prob * self.subject_weight[subject]
                for subject, prob in aln_prob.items()
            ]
            new_probs_sum = sum(new_probs)
            for ix, subject in enumerate(aln_prob):
                self.subjects_to_update.add(subject)
                v = new_probs[ix] / new_probs_sum
                self.aln_prob[query][subject] = v
                self.aln_prob_T[subject][query] = v

    def trim_least_likely(self, scale=0.9):
        """Remove the least likely alignments."""
        logging.info("Removing alignments below {} of the max".format(scale))
        n_trimmed = 0
        newly_unique_queries = set([])
        for query in self.multimapped_queries:
            aln_prob = self.aln_prob[query]
            # Skip queries with no remaining alignments
            if len(aln_prob) == 0:
                continue
            # Skip queries with only a single possible subject
            if len(aln_prob) == 1:
                newly_unique_queries.add(query)
                self.n_unique += 1
                continue
            # Figure out our maximum score for this query
            max_likely = max(list(aln_prob.values()))

            # Trim anyone BELOW the maximum possible value for this query.
            to_remove = [
                subject for subject, prob in aln_prob.items()
                if prob < scale * max_likely
            ]

            n_trimmed += len(to_remove)

            # Remove the subjects
            for subject in to_remove:
                del self.aln_prob[query][subject]
                del self.aln_prob_T[subject][query]

            if len(self.aln_prob[query]) == 1:
                newly_unique_queries.add(query)
                self.n_unique += 1

        self.multimapped_queries = self.multimapped_queries - newly_unique_queries

        logging.info("Removed {:,} unlikely alignments".format(n_trimmed))
        logging.info("Number of uniquely aligned queries: {:,}".format(
            self.n_unique))
        return n_trimmed


def filter_subjects_by_coverage(args):
    """Check whether the subject passes the coverage filter."""

    subject, cov, SD_MEAN_CUTOFF, STRIM_5, STRIM_3 = args

    # Trim the ends
    if cov.shape[0] >= STRIM_5 + STRIM_3 + 10:
        cov = cov[STRIM_5: -STRIM_3]

    passes_filter = cov.mean() > 0 and cov.std() / cov.mean() <= SD_MEAN_CUTOFF
    return subject, passes_filter


def calc_cov_by_subject(alignments, subject_len):
    """Index a set of sorted alignments by subject."""
    # The output will be a dict with the start:stop index for each subject
    assert len(alignments) > 0
    index = {}

    coverages = {
        subject: np.zeros(length, dtype=int)
        for subject, length in subject_len.items()
    }

    last_subject = None
    last_start_ix = None
    for ix, a in enumerate(alignments):
        # query, subject, sstart, send, bitstore = a
        if a[1] != last_subject:
            if last_subject is not None:
                index[last_subject] = (last_start_ix, ix)
            last_subject = a[1]
            last_start_ix = ix

        # Add to the cov_proc
        coverages[a[1]][a[2]:a[3]] += 1

    index[last_subject] = (last_start_ix, ix)

    return coverages, index


def parse_alignment(align_handle,
                    batchsize=None,   # Number of reads to process at a time
                    QSEQID_i=0,
                    SSEQID_i=1,
                    SSTART_i=8,
                    SEND_i=9,
                    BITSCORE_i=11,
                    SLEN_i=13,
                    SD_MEAN_CUTOFF=1.0,
                    STRIM_5=18,
                    STRIM_3=18,
                    threads=4,
                    MAX_ITERATIONS=1000):
    """
    Parse an alignment in BLAST6 format and determine which subjects are likely to be present. This is the core of FAMLI.
    BLAST 6 columns by default (in order): qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
                                              0     1       2       3       4       5       6   7   8       9   10      11      12   13  
    """
    # Initialize the multithreading pool before reading in the alignments
    pool = Pool(threads)

    # Initialize the alignment parser
    parser = BLAST6Parser()

    # Keep a list of the final, filtered alignments
    final_alignments = []

    # Iterate over blocks of alignments
    for alignments in parser.yield_alignments(
        align_handle,
        batchsize=batchsize,
        QSEQID_i=QSEQID_i,
        SSEQID_i=SSEQID_i,
        SSTART_i=SSTART_i,
        SEND_i=SEND_i,
        BITSCORE_i=BITSCORE_i,
        SLEN_i=SLEN_i
    ):
        if len(alignments) == 0:
            logging.info("Zero alignments to process in this batch")
            continue

        # Count the total number of reads that were aligned
        logging.info("Counting unique queries")
        aligned_reads = len(set([a[0] for a in alignments]))

        logging.info("Number of unique queries: {:,}".format(aligned_reads))
        logging.info("Number of alignments: {:,}".format(len(alignments)))

        # Sort alignments by subject
        logging.info("Sorting alignments by subject")
        alignments.sort(key=lambda a: a[1])

        logging.info("Calculating coverage by subject")
        subject_coverages, subject_index = calc_cov_by_subject(
            alignments, parser.subject_len)

        # STEP 1. FILTER SUBJECTS BY "COULD" COVERAGE EVENNESS
        logging.info("FILTER 1: Even coverage of all alignments")

        filter_1 = pool.map(filter_subjects_by_coverage, [
            [
                subject,
                coverage,
                SD_MEAN_CUTOFF,
                STRIM_5,
                STRIM_3
            ]
            for subject, coverage in subject_coverages.items()
        ])

        # Reformat as a dict
        filter_1 = dict(filter_1)

        logging.info("Subjects passing FILTER 1: {:,} / {:,}".format(
            sum(filter_1.values()), len(filter_1)
        ))

        # Subset the alignments to only those subjects passing the filter
        alignments = [
            a
            for subject, start_stop in subject_index.items()
            for a in alignments[start_stop[0]: start_stop[1]]
            if filter_1[subject]
        ]

        logging.info("Queries passing FILTER 1: {:,} / {:,}".format(
            len(set([a[0] for a in alignments])), aligned_reads
        ))

        # STEP 2: Reassign multi-mapped reads to a single subject
        logging.info("FILTER 2: Reassign queries to a single subject")

        # Add the alignments to a model to optimally re-assign reads
        model = FAMLI_Reassignment(alignments, parser.subject_len, pool=pool)

        # Initialize the subject weights
        model.init_subject_weight()

        ix = 0
        while ix <= MAX_ITERATIONS:
            ix += 1
            logging.info("Iteration: {:,}".format(ix))
            # Recalculate the subject weight, given the naive alignment probabliities
            model.recalc_subject_weight()
            # Recalculate the alignment probabilities, given the subject weights
            model.recalc_aln_prob()

            # Trim the least likely alignment for each read
            n_trimmed = model.trim_least_likely()

            if n_trimmed == 0:
                break

        logging.info("Subsetting the alignment to only the reassigned queries")

        # Subset the alignment to only the reassigned queries
        alignments = [
            a
            for a in alignments
            if a[1] in model.aln_prob[a[0]] and
            len(model.aln_prob[a[0]]) == 1
        ]

        logging.info("Finished reassigning reads ({:,} remaining)".format(
            len(alignments)))

        # Add to the batch of final alignments
        final_alignments.extend(alignments)
        logging.info("A total of {:,} filtered alignments collected".format(
            len(final_alignments)))
        del alignments

    if len(final_alignments) == 0:
        logging.info("The entire sample contains zero alignments")
        # Return some empty data
        return 0, []

    # Change the name, for brevity
    alignments = final_alignments
    queries_after_reassignment = len(alignments)

    # Cull the subjects that were removed during FILTER 2
    all_subjects = set([a[1] for a in alignments])
    parser.subject_len = {
        subject: length
        for subject, length in parser.subject_len.items()
        if subject in all_subjects
    }

    # STEP 3: Filter subjects by even coverage of reassigned queries
    logging.info("FILTER 3: Filtering subjects by even sequencing coverage")

    subject_coverages, subject_index = calc_cov_by_subject(
        alignments, parser.subject_len)

    filter_3 = pool.map(filter_subjects_by_coverage, [
        [
            subject,
            coverage,
            SD_MEAN_CUTOFF,
            STRIM_5,
            STRIM_3
        ]
        for subject, coverage in subject_coverages.items()
    ])
    # Reformat as a dict
    filter_3 = dict(filter_3)

    logging.info("Subjects passing FILTER 3: {:,} / {:,}".format(
        sum(filter_3.values()), len(filter_3)
    ))

    # Subset the alignments to only those subjects passing the filter
    alignments = [
        a
        for subject, start_stop in subject_index.items()
        for a in alignments[start_stop[0]: start_stop[1]]
        if filter_3[subject]
    ]

    logging.info("Queries passing FILTER 3: {:,} / {:,}".format(
        len(set([a[0] for a in alignments])), queries_after_reassignment
    ))

    # Make the output by calculating coverage per subject
    output = []

    # Make a dict of the alignment ranges
    logging.info("Collecting final alignments")
    alignment_ranges = defaultdict(list)
    for query, subject, sstart, send, bitscore in alignments:
        alignment_ranges[subject].append((sstart, send))

    logging.info("Calculating final stats")

    for subject, aln_ranges in alignment_ranges.items():
        # Make a coverage map
        cov = np.zeros(model.subject_len[subject], dtype=int)

        # Add to the coverage
        for sstart, send in aln_ranges:
            cov[sstart:send] += 1

        output.append({
            "id": subject,
            "nreads": len(aln_ranges),
            "coverage": (cov > 0).mean(),
            "depth": cov.mean(),
            "std": cov.std(),
            "length": model.subject_len[subject],
        })

    logging.info("Results: assigned {:,} queries to {:,} subjects".format(
        sum([d["nreads"] for d in output]),
        len(output),
    ))

    return aligned_reads, output
