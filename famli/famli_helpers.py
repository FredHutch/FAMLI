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
        burn_in=1000000,
        SD_MEAN_CUTOFF=1.0,
        QSEQID_i=0,
        SSEQID_i=1,
        SSTART_i=8,
        SEND_i=9,
        EVALUE_i=10,
        BITSCORE_i=11,
        SLEN_i=13,
    ):

        # Number of reads to process before optimizing
        self.burn_in = burn_in

        # Threshold to prune unlikely alignments
        self.SD_MEAN_CUTOFF = SD_MEAN_CUTOFF

        # Pool of workers
        # self.pool = Pool(threads)

        # Keep track of the weight for every subject
        self.subject_weight = defaultdict(float)

        # Keep track of the bitscore for each query ~ subject
        self.query_bitscore = defaultdict(dict)

        # Keep track of the alignments for every query (sstart, send)
        self.alignments = defaultdict(dict)

        # Number of total and uniquely aligned queries
        self.n_queries = 0
        self.n_unique = 0

        # Keep a flag indicating whether any suboptimal alignments have been pruned
        self.any_pruned = False
        
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
            if cov.std() / cov.mean() <= self.SD_MEAN_CUTOFF:
                to_keep.add(subject)

        return to_keep

    def parse(self, align_handle):
        """Parse a set of reads, optimizing as we go."""

        # Get the subjects which pass the first evenness filter
        filter_1_subjects = self.filter_subjects_evenness(align_handle)
        logging.info("Subjects passing FILTER 1: {:,}".format(
            len(filter_1_subjects)
        ))

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
            # Filter to the subjects passing FILTER 1
            query_alignments = [
                a for a in query_alignments
                if a[1] in filter_1_subjects
            ]

            # Skip if no subjects pass filter 1
            if len(query_alignments) == 0:
                continue

            # Adjust the bitscores to add up to 1
            tot_bitscore = np.sum([a[4] for a in query_alignments])

            # Add those bitscores to the running totals
            for query_name, subject, sstart, send, bitscore in query_alignments:
                # Weights sum to 1
                bitscore = bitscore / tot_bitscore
                
                # Subject weight
                self.subject_weight[subject] += bitscore
                # Query ~ subject bitscore
                self.query_bitscore[query_name][subject] = bitscore
                # Add the alignment information
                self.alignments[query_name][subject] = (sstart, send)

            # Add to the counter for total queries
            self.n_queries += 1

            # Keep track of uniquely aligned queries
            if len(query_alignments) == 1:
                self.n_unique += 1

            # If we're past the burn-in, try to optimize this query
            if self.n_queries > self.burn_in:
                self.optimize_one_query(query_name)
            # If we've just now reached the end of the burn-in, optimize everything
            elif self.n_queries == self.burn_in:
                logging.info("Finished burn-in of {:,} queries".format(self.burn_in))
                self.optimize_all()
                logging.info("Total queries read in: {:,}".format(self.n_queries))
                logging.info("Uniquely aligned queries: {:,}".format(self.n_unique))

            # Otherwise, just keep on adding queries

        # Once all queries have been added, keep optimizing until no longer possible
        n_rounds = 1
        logging.info("Number of unique queries: {:,}".format(self.n_unique))
        while n_rounds <= 1000:
            self.any_pruned = False
            logging.info("Round {:,}: Reassigning multi-mapped queries".format(n_rounds))
            self.optimize_all()
            logging.info("Number of unique queries: {:,}".format(self.n_unique))
            n_rounds += 1
            if self.any_pruned is False:
                break
        logging.info("Done reassigning reads")

    def optimize_one_query(self, query_name):
        """Remove suboptimal alignments for a single query."""

        # Only optimize if there is more than one possible subject
        if len(self.alignments[query_name]) == 1:
            return

        # Get the likelihoods for this query
        # Likelihood = bitscore * evenness score * total reference weight / reference length
        likelihoods = {
            subject: bitscore * self.subject_weight[subject] / self.alignment_parser.subject_len[subject]
            for subject, bitscore in self.query_bitscore[query_name].items()
        }
        # Calculate the maximum likelihood
        max_likelihood = np.max(likelihoods.values())

        # Get the bitscores for those subjects with likelihoods above the threshold
        new_bitscores = {
            subject: bitscore
            for subject, bitscore in self.query_bitscore[query_name].items()
            if likelihoods[subject] >= max_likelihood
        }

        assert len(new_bitscores) > 0

        # If none of the subjects are removed, don't do anything
        if len(new_bitscores) == len(self.query_bitscore[query_name]):
            return
        else:
            # Flag indicating that a suboptimal alignment was pruned
            self.any_pruned = True

            # Adjust the new bitscores to sum to 1
            new_bitscore_sum = np.sum(new_bitscores.values())
            new_bitscores = {
                subject: bitscore / new_bitscore_sum
                for subject, bitscore in new_bitscores.items()
            }

            # Newly unique query
            if len(new_bitscores) == 1:
                self.n_unique += 1

            # Now adjust the subject weights
            for subject, old_bitscore in self.query_bitscore[query_name].items():
                # If this subject has been eliminated
                if subject not in new_bitscores:
                    # Lower the total score by the amount for this query
                    new_score = self.subject_weight[subject] - old_bitscore

                    # Remove the alignment information
                    sstart, send = self.alignments[query_name][subject]

                    # Remove the coverage information
                    del self.alignments[query_name][subject]
                else:
                    # Calculate the new score based on the change in this query's contribution
                    new_score = self.subject_weight[subject] + new_bitscores[subject] - old_bitscore

                # Assign the new score
                self.subject_weight[subject] = new_score

            # And update the per-query bitscores
            self.query_bitscore[query_name] = new_bitscores

    def optimize_all(self):
        """Remove suboptimal alignments for all queries."""
        logging.info("Optimizing all {:,} queries".format(self.n_queries))

        for query_name in self.query_bitscore.keys():
            self.optimize_one_query(query_name)

    def summary(self):
        """Calculate coverage metrics for every reference."""

        # Calculate the coverage using uniquely assigned queries
        coverage = {}
        nreads = defaultdict(int)

        for query_name, bitscores in self.query_bitscore.items():
            if len(bitscores) != 1:
                continue
            subject = list(bitscores.keys())[0]
            sstart, send = self.alignments[query_name][subject]

            if subject not in coverage:
                slen = self.alignment_parser.subject_len[subject]
                coverage[subject] = np.zeros(int(slen), dtype=int)
            
            coverage[subject][sstart: send] += 1
            nreads[subject] += 1
        
        # Output the coverage metrics
        output = [{
                "id": subject,
                "nreads": n,
                "depth": coverage[subject].mean(),
                "std": coverage[subject].std(),
                "length": self.alignment_parser.subject_len[subject]
            }
            for subject, n in nreads.items()
        ]

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

        return self.alignment_parser.n_unique_queries, output
