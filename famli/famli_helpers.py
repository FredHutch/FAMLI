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

    def __init__(
        self,
        QSEQID_i=0,
        SSEQID_i=1,
        SSTART_i=8,
        SEND_i=9,
        EVALUE_i=10,
        BITSCORE_i=11,
        SLEN_i=13
    ):
        # Keep track of the length of all subjects
        self.subject_len = {}
        # Keep track of the number of queries
        self.unique_queries = set([])

        # Column index values
        self.QSEQID_i = QSEQID_i
        self.SSEQID_i = SSEQID_i
        self.SSTART_i = SSTART_i
        self.SEND_i = SEND_i
        self.EVALUE_i = EVALUE_i
        self.BITSCORE_i = BITSCORE_i
        self.SLEN_i = SLEN_i

    def parse(self, align_handle):
        """Parse a file, while keeping track of the subject lengths."""
        logging.info("Reading in alignments")

        for i, line in enumerate(align_handle):
            if i % 1000000 == 0 and i > 0:
                logging.info("{:,} lines of alignment parsed".format(i))
            line_list = line.rstrip("\n").split("\t")
            # Get the query and subject
            query = line_list[self.QSEQID_i]
            subject = line_list[self.SSEQID_i]
            bitscore = float(line_list[self.BITSCORE_i])
            sstart = int(line_list[self.SSTART_i])
            send = int(line_list[self.SEND_i])

            # Save information for the subject length
            if subject not in self.subject_len:
                self.subject_len[subject] = float(line_list[self.SLEN_i])

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

    def yield_queries(self, align_handle):
        """Yield all of the alignments for a single query."""
        logging.info("Reading in alignments, grouping by query")

        # Keep track of the number of unique queries
        self.n_unique_queries = 0
        # Keep track of the last unique query ID
        last_query = None
        # Keep a buffer of alignments for each query
        alignments = []

        # Add to the alignments, stopping when we get to the batch size limit
        for a in self.parse(align_handle):

            # When we reach the batch size and new query, yield the alignments
            if a[0] != last_query and len(alignments) > 0:
                yield alignments
                # Reset the alignment data
                alignments = []
                self.n_unique_queries += 1

            # Add to the alignment list
            alignments.append(a)

            # Record which query we saw in this alignment
            last_query = a[0]

        # Yield the final block of alignments
        self.n_unique_queries += 1
        logging.info("Total queries read in: {:,}".format(self.n_unique_queries))
        yield alignments


class FAMLI_Reassignment:
    """Reassign alignments to a single reference using FAMLI."""

    def __init__(
        self,
        threads=4,
        SD_MEAN_CUTOFF=1.0,
        QSEQID_i=0,
        SSEQID_i=1,
        SSTART_i=8,
        SEND_i=9,
        EVALUE_i=10,
        BITSCORE_i=11,
        SLEN_i=13,
    ):

        # Threshold to prune unlikely alignments
        self.SD_MEAN_CUTOFF = SD_MEAN_CUTOFF

        # Keep track of the weight for every subject
        self.subject_weight = defaultdict(float)

        # Keep track of the coverage vectors per subject
        self.coverage = {}

        # Number of total and uniquely aligned queries
        self.n_queries = 0
        self.n_unique = 0

        # Keep track of the number of reads per subject
        self.subject_n_reads = defaultdict(int)

        # Column index values
        self.QSEQID_i = QSEQID_i
        self.SSEQID_i = SSEQID_i
        self.SSTART_i = SSTART_i
        self.SEND_i = SEND_i
        self.EVALUE_i = EVALUE_i
        self.BITSCORE_i = BITSCORE_i
        self.SLEN_i = SLEN_i

    def filter_subjects_evenness(self, align_handle):
        """Return the set of subjects with coverage evenness below the threshold."""
        # Track the coverage for each subject
        coverage = {}

        parser = BLAST6Parser(
            QSEQID_i = self.QSEQID_i,
            SSEQID_i = self.SSEQID_i,
            SSTART_i = self.SSTART_i,
            SEND_i = self.SEND_i,
            EVALUE_i = self.EVALUE_i,
            BITSCORE_i = self.BITSCORE_i,
            SLEN_i = self.SLEN_i
        )
        for query_alignments in parser.yield_queries(align_handle):
            # Weight the coverage by the number of alignments
            alignment_weight = 1. / len(query_alignments)
            for _, subject, sstart, send, _ in query_alignments:
                if subject not in coverage:
                    coverage[subject] = np.zeros(
                        int(parser.subject_len[subject]),
                        dtype=np.float
                    )
                coverage[subject][sstart: send] += alignment_weight

        # Return the set of subjects that are below the evenness cutoff
        to_keep = set([])
        for subject, cov in coverage.items():
            # Trim the ends
            if cov.shape[0] >= 46:
                cov = cov[18: -18]
            if cov.max() > 0 and cov.std() / cov.mean() <= self.SD_MEAN_CUTOFF:
                to_keep.add(subject)

        return to_keep

    def calc_init_subject_weights(self, align_handle):
        """Calculate the initial set of subject weights."""
        # Set up the parser
        self.alignment_parser = BLAST6Parser(
            QSEQID_i = self.QSEQID_i,
            SSEQID_i = self.SSEQID_i,
            SSTART_i = self.SSTART_i,
            SEND_i = self.SEND_i,
            EVALUE_i = self.EVALUE_i,
            BITSCORE_i = self.BITSCORE_i,
            SLEN_i = self.SLEN_i
        )
        for query_alignments in self.alignment_parser.yield_queries(align_handle):
            # Filter to the subjects passing FILTER 1
            query_alignments = [
                a for a in query_alignments
                if a[1] in self.filter_1_subjects
            ]

            # Skip if no subjects pass filter 1
            if len(query_alignments) == 0:
                continue

            # Adjust the bitscores to add up to 1
            tot_bitscore = np.sum([a[4] for a in query_alignments])

            # Add those bitscores to the running totals
            for _, subject, _, _, bitscore in query_alignments:
                # Weights sum to 1
                bitscore = bitscore / tot_bitscore
                
                # Subject weight
                self.subject_weight[subject] += bitscore

    def parse(self, align_handle):
        """Parse a set of reads, optimizing as we go."""

        # Get the subjects which pass the first evenness filter
        self.filter_1_subjects = self.filter_subjects_evenness(align_handle)
        logging.info("Subjects passing FILTER 1: {:,}".format(
            len(self.filter_1_subjects)
        ))

        # Get the initial weights for those subjects passing filter 1
        align_handle.seek(0)
        self.calc_init_subject_weights(align_handle)

        # Now go through and greedily assign reads to subjects
        self.alignment_parser = BLAST6Parser(
            QSEQID_i = self.QSEQID_i,
            SSEQID_i = self.SSEQID_i,
            SSTART_i = self.SSTART_i,
            SEND_i = self.SEND_i,
            EVALUE_i = self.EVALUE_i,
            BITSCORE_i = self.BITSCORE_i,
            SLEN_i = self.SLEN_i
        )
        align_handle.seek(0)
        for query_alignments in self.alignment_parser.yield_queries(align_handle):
            # Increment the number of total queries
            self.n_queries += 1

            # Calculate the likelihood for each subject
            # Filter to the subjects passing FILTER 1
            # Likelihood = bitscore * subject weight / subject length
            likelihood = {
                subject: bitscore * self.subject_weight[subject] / self.alignment_parser.subject_len[subject]
                for _, subject, _, _, bitscore in query_alignments
                if subject in self.filter_1_subjects
            }

            # Skip if no subjects pass filter 1
            if len(likelihood) == 0:
                continue

            # Maximum likelihood value
            max_likelihood = max(likelihood.values())

            # Sum of bitscores
            tot_bitscore = sum([
                bitscore
                for _, subject, _, _, bitscore in query_alignments
                if subject in self.filter_1_subjects
            ])

            # Skip if more than one subject has the top hit
            if sum([v == max_likelihood for v in likelihood.values()]) != 1:
                continue

            # Assign the query to the most likely subject
            for _, subject, sstart, send, bitscore in query_alignments:
                # Skip subject not passing filter 1
                if subject not in self.filter_1_subjects:
                    continue

                # Weights sum to 1
                bitscore = bitscore / tot_bitscore
                
                # This is the TOP HIT
                if likelihood[subject] == max_likelihood:
                    # Increment the counter of unique reads
                    self.n_unique += 1

                    # Adjust the subject weight
                    self.subject_weight[subject] += (1 - bitscore)

                    # Add to the coverage vectors
                    if subject not in self.coverage:
                        self.coverage[subject] = np.zeros(
                            int(self.alignment_parser.subject_len[subject]),
                            dtype=int
                        )

                    # Increment the coverage vector
                    self.coverage[subject][sstart: send] += 1

                    # Increment the read counter
                    self.subject_n_reads[subject] += 1

                else:
                    # This is a suboptimal hit
                    self.subject_weight[subject] -= bitscore

            if self.n_queries % 100000 == 0:
                logging.info("Total queries read in: {:,}".format(self.n_queries))
                logging.info("Uniquely aligned queries: {:,}".format(self.n_unique))

        logging.info("Done reassigning reads")

    def summary(self):
        """Calculate coverage metrics for every reference."""
        
        # Output the coverage metrics
        output = [{
                "id": subject,
                "nreads": n,
                "depth": self.coverage[subject].mean(),
                "std": self.coverage[subject].std(),
                "length": self.alignment_parser.subject_len[subject]
            }
            for subject, n in self.subject_n_reads.items()
        ]

        assert sum([a["nreads"] for a in output]) == self.n_unique

        logging.info("Number of queries passing FILTER 2: {:,}".format(
            sum([a["nreads"] for a in output])
        ))
        logging.info("Number of subjects passing FILTER 2: {:,}".format(
            len(output)
        ))

        # FILTER 3: Evenness
        output = [
            a for a in output
            if a["std"] / a["depth"] <= self.SD_MEAN_CUTOFF
        ]

        logging.info("Number of queries passing FILTER 3: {:,}".format(
            sum([a["nreads"] for a in output])
        ))
        logging.info("Number of subjects passing FILTER 3: {:,}".format(
            len(output)
        ))

        return self.n_queries, output
