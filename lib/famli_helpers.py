#!/usr/bin/python

import json
import logging
import pandas as pd
import numpy as np
from math import ceil
from scipy.stats import binom
from Bio.Data import CodonTable
from .exec_helpers import run_cmds
from itertools import permutations
from collections import defaultdict


class ErrorModelFAMLI:
    """Error model for the likelihood of non-synonymous nucleotide errors."""

    def __init__(self, error_rate, genetic_code=11):
        """For the genetic code, calculate the proportion of non-syn errors."""

        # The likelihood that a given read is truly from a given reference,
        # given only the information from a single assignment, is defined
        # by the number of mismatches and the sequence length (using the
        # binomial)

        # The likelihood that a given read is truly from a given reference,
        # given all of the information we have access to about the reference
        # that the read _could_ align to is defined as the likelihood due to
        # error, divided by the sum of the likelihood that it is truly from any
        # other reference due to random sequencing error

        # Keep track of the likelihood that any single read is truly from a
        # given reference with the dict of dicts, psingle

        self.psingle = defaultdict(lambda: defaultdict(float))

        # Keep track of the names of each input
        self.query_names = {}

        # For this genetic code, determine the proportion of nucleotide
        # substitutions that result in a change in the encoded amino acid
        table = CodonTable.unambiguous_dna_by_id[genetic_code]

        all_codons = list(set(list(permutations('ATCGATCGATCG', 3))))
        assert len(all_codons) == 64

        # The total number of substitutions
        n_subs = 0
        # The number of non-synonymous substitutions
        n_nonsyn_subs = 0
        for codon1 in all_codons:
            codon1 = ''.join(codon1)
            if codon1 in table.stop_codons:
                aa1 = "X"
            else:
                aa1 = table.forward_table[codon1]

            for codon2 in all_codons:
                codon2 = ''.join(codon2)
                if codon2 in table.stop_codons:
                    aa2 = "X"
                else:
                    aa2 = table.forward_table[codon2]
                # Check to see if there is a single substitution
                if sum([c1 == c2 for c1, c2 in zip(codon1, codon2)]) == 2:
                    n_subs += 1
                    if aa1 != aa2:
                        n_nonsyn_subs += 1

        nonsyn_rate = n_nonsyn_subs / float(n_subs)
        msg = "With genetic code {}, the non-synonymous proportion is {}"
        logging.info(msg.format(genetic_code, nonsyn_rate))

        # Given the nucleotide error rate, calculate the
        # effective amino acid error rate
        self.aa_error_rate = error_rate * 3. * nonsyn_rate
        logging.info("Nucleotide error rate: {}".format(error_rate))
        logging.info("Amino acid error rate: {}".format(self.aa_error_rate))

        # Keep a cache with the proability of AT LEAST {N} subsitutions
        # in a sequence of {LEN} amino acids in length
        self.likelihood_cache = {}

    def edit_dist_prob(self, n_mismatch, seq_len):
        """Calculate the probability of mismatches due to random error."""
        assert n_mismatch < seq_len
        # Add the likelihood values to the cache, if not already present
        if seq_len not in self.likelihood_cache:
            # Calculate the probability of exactly N substitutions
            pmf = [
                binom.pmf(N, seq_len, self.aa_error_rate)
                for N in range(seq_len)
            ]
            # Calculate the likelihood of AT LEAST N subsitutions
            sum_pmf = [
                sum(pmf[ix:])
                for ix in range(len(pmf))
            ]
            self.likelihood_cache[seq_len] = sum_pmf

        return self.likelihood_cache[seq_len][n_mismatch]

    def add_prob_to_alignments(self, alignments, query_ix):
        """Add the relative likelihood metric to a set of alignments."""
        # This is to be run for all alignments BEFORE calculating
        # the maximum likelihood for any individual alignment
        # Because: this function adds to the total probability
        # mass for each reference, which is used to calculate the
        # maximum likelihood in the function below

        for a in alignments:
            a["likelihood"] = self.edit_dist_prob(a["n_mismatch"], a["len"])

            # Set the query name
            self.query_names[query_ix] = a["query"]

        # Adjust the likelihood values to sum to 1
        likelihood_sum = np.sum([a["likelihood"] for a in alignments])
        for a in alignments:
            a["likelihood"] = a["likelihood"] / likelihood_sum

            # Add to the probability mass for this reference
            self.psingle[query_ix][a["ref"]] = a["likelihood"]

    def calc_max_likelihood(self, alignments):
        """Determine which alignment has the maximum likelihood."""
        # This is to be run AFTER all of the individual alignments
        # have had their individual proabilities calculated

        # Calculate the weighted likelihood for each
        weights = [
            a["likelihood"] * self.ref_prob_mass[a["ref"]]
            for a in alignments
        ]

        max_weight = max(weights)

        return [
            a for ix, a in enumerate(alignments)
            if weights[ix] == max_weight
        ]

    def optimize_assignments(self, cutoff=0.05):
        """Iteratively optimize the assignments of reads."""

        keep_iterating = True
        n_removed = 0
        iterations = 0
        while keep_iterating:
            iterations += 1

            print("\n\nITERATION {:,} ({:,} removed)".format(iterations, n_removed))
            keep_iterating = False

            # Calculate the estimated abundance of each reference
            ref_abund = defaultdict(float)
            for query, query_ref_abund in self.psingle.items():
                for ref, abund in query_ref_abund.items():
                    ref_abund[ref] += abund

            print(pd.Series(ref_abund).sort_values().tail())

            # Calculate the total likelihood for all queries
            # as the single probability * reference abundance
            # ptotal = defaultdict(lambda: defaultdict(float))

            # for query, query_single_probs in self.psingle.items():
            #     for ref, prob in query_single_probs.items():
            #         if prob > 0:
            #             ptotal[query][ref] = prob * ref_abund[ref]

            #     # Divide by the total
            #     t = sum(ptotal[query].values())
            #     for ref, prob in ptotal[query].items():
            #         ptotal[query][ref] = prob / t

            # print(pd.DataFrame({
            #     query: pd.Series(ptotal[query]).sort_values().describe()
            #     for query in ptotal.keys()
            # }).T.describe())

            # # Start with the lowest likelyhood and iterate upwards
            # all_probs = sorted(list(set([
            #     prob
            #     for query_probs in ptotal.values()
            #     for prob in query_probs.values()
            # ])))

            # # Check to see if any of the query - references meet the criteria
            # # and if so, remove it from self.psingle
            # lowest_removed = None
            # highest_removed = None

            # assigned = defaultdict(int)
            # for query, query_total_probs in ptotal.items():
            #     for ref, p in query_total_probs.items():
            #         if round(p, 2) == 1.00:
            #             assigned[ref] += 1

            # assigned = pd.Series(assigned)
            # assigned.sort_values(ascending=False, inplace=True)
            # print(assigned.head())

            # Go through the queries and remove the aggregate set of alignments
            # that fall below the cutoff

            for query, ref_prob in self.psingle.items():

                # Calculate the reference probability based on the total sample
                ptotal = pd.Series({
                    ref: prob * ref_abund[ref]
                    for ref, prob in ref_prob.items()
                }).sort_values()

                # Divide by the sum
                ptotal = ptotal / ptotal.sum()

                if query == 2:
                    print(self.query_names[query])
                    print(ptotal)

                # Calculate the aggregate probabilty
                
                for ref, agg_prop in ptotal.cumsum()[:ceil(ptotal.shape[0] / float(5))].items():
                    # If the aggregate probability is below the cutoff
                    if agg_prop <= cutoff:
                        # Disregard the likelihood of this reference
                        del self.psingle[query][ref]
                        n_removed += 1
                        keep_iterating = True
                        t = sum(self.psingle[query].values())
                        self.psingle[query] = {k: v / t for k, v in self.psingle[query].items()}

            # for least_likely in np.arange(0, cutoff, cutoff / 100.):
            #     # If the least likely is above the cutoff, skip
            #     if least_likely == 0:
            #         continue

            #     for query, query_total_probs in ptotal.items():
            #         assert round(sum(query_total_probs.values()), 2) == 1.00, sum(query_total_probs.values())

            #         # If everything is the least likely, skip
            #         if all([
            #             p <= least_likely
            #             for p in query_total_probs.values()
            #         ]):
            #             continue

            #         x = False
            #         # Remove the least likely references
            #         for ref, p in query_total_probs.items():
            #             if p <= least_likely:
            #                 # print("Removing {} @ {}".format(ref, p))
            #                 keep_iterating = True
            #                 self.psingle[query][ref] = 0
            #                 n_removed += 1
            #                 x = True
            #         if x:
            #             self.psingle[query] = {
            #                 k: v / sum(self.psingle[query].values())
            #                 for k, v in self.psingle[query].items()
            #             }

            #     # If any references were removed, recalculate ptotal
            #     if keep_iterating:
            #         # assert lowest_removed is not None
            #         # assert highest_removed >= lowest_removed
            #         # if highest_removed > lowest_removed * 2:
            #         #     break
            #         break

        # Calculate the estimated abundance of each reference
        ref_abund = defaultdict(float)
        for query, query_ref_abund in self.psingle.items():
            for ref, abund in query_ref_abund.items():
                ref_abund[ref] += abund

        # Loop through the reads
        final_counts = defaultdict(int)
        assignments = {}
        for query, query_ref_abund in self.psingle.items():
            # Calculate the total probability
            # single probability * reference abundance
            ptotal = {
                ref: ref_single_prob * ref_abund[ref]
                for ref, ref_single_prob in query_ref_abund.items()
                if ref_single_prob > 0
            }
            # Divide by the total
            t = sum(ptotal.values())

            ptotal = {
                k: v / t
                for k, v in ptotal.items()
            }

            # print(query, len(query_ref_abund), len(ptotal))
            if len(ptotal) == 1:
                final_counts[list(ptotal.keys())[0]] += 1

                assignments[query] = list(ptotal.keys())[0]
            else:
                print(pd.Series(ptotal).sort_values(ascending=False))

        final_counts = pd.Series(final_counts)

        for query, assignment in assignments.items():
            if final_counts[assignment] <= 5:
                print(query, assignment)

        print("\n\nIterations: {:,}".format(iterations))
        print("\n\nInput: {:,}".format(sum([len(v.values()) for v in self.psingle.values()])))
        print("\n\nRemoved: {:,}".format(n_removed))
        print("\n\nOutput: {:,}".format(sum([sum([x > 0 for x in v.values()]) for v in self.psingle.values()])))
        print("\n\nDuplicated: {:,}".format(sum([sum([x > 0 for x in v.values()]) > 1 for v in self.psingle.values()])))
        print("\n\nDeduplicated: {:,}".format(sum([sum([x > 0 for x in v.values()]) == 1 for v in self.psingle.values()])))
        print("\n\nUnique references: {}".format(final_counts.shape[0]))
        print(final_counts.sort_values(ascending=False))

        # print("\n".join(list(set([
        #     ref
        #     for query, query_probs in self.psingle.items()
        #     for ref, ref_prob in query_probs.items()
        #     if ref_prob > 0
        # ]))))


def align_reads(read_fp,               # FASTQ file path
                db_fp,                 # Local path to DB
                temp_folder,           # Folder for results
                query_gencode=11,
                threads=1):
    """Align a set of reads with Paladin."""

    align_fp = "{}.sam".format(read_fp)
    logging.info("Input reads: {}".format(read_fp))
    logging.info("Reference database: {}".format(db_fp))
    logging.info("Genetic code: {}".format(query_gencode))
    logging.info("Threads: {}".format(threads))
    logging.info("Output: {}".format(align_fp))

    run_cmds([
            "paladin",
            "align",
            "-a",
            "-t",
            str(threads),
            "-z",
            str(query_gencode),
            db_fp,
            read_fp
        ], stdout=align_fp)

    return align_fp


def parse_alignment(align_fp, error_rate=0.001, genetic_code=11):
    """Parse an alignment in SAM format."""

    # Initialize the error model
    error_model = ErrorModelFAMLI(error_rate, genetic_code=genetic_code)

    # Read over the file once to get all of the alignments
    # Group the alignments by query sequence
    alignments = [a for a in parse_alignments_by_query(align_fp)]

    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    # ONLY USE A FEW ALIGNMENTS
    # alignments = alignments[:200]
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    # Keep track of the total set of reference with any alignments
    all_refs = set([
        a["ref"]
        for query in alignments
        for a in query
    ])

    # Length of each reference sequence
    ref_length = {}

    # Read just the header to get the length of each reference
    with open(align_fp, "rt") as f:
        for line in f:
            # Marker for reference length lines
            if line[:3] == "@SQ":
                # Split up the line
                line = line.rstrip("\n").split("\t")
                # Reference name is in the second field
                assert len(line) == 3, line
                ref = line[1][3:]
                # Check if this reference has any alignments
                if ref in all_refs:
                    length = int(line[2][3:])
                    # Add to the dict
                    ref_length[ref] = length

    # Add every alignment into the error model
    for query_ix, a in enumerate(alignments):
        error_model.add_prob_to_alignments(a, query_ix)

    # Optimize the read assignments, iteratively calculating the
    # total reference abunance and removing potential alignments
    # that fall below the specified level of likelihood
    error_model.optimize_assignments(cutoff=0.5)

    return None, None, None

    # Keep track of the number of reads and coverage across each reference
    # Keep track of the RAW stats, as well as the deduplicated alignments

    # Number of reads
    total_reads = 0
    deduplicated_reads = 0
    raw_reads = defaultdict(int)
    final_reads = defaultdict(int)
    # Coverage across each reference
    raw_cov = {}
    final_cov = {}

    # Now go over the alignments again and reassign each read
    # to the most likely reference
    for query in alignments:
        total_reads += 1
        # Add to the raw abundances
        for a in query:
            raw_reads[a["ref"]] += 1

            # Add to the coverage metrics
            if a["ref"] not in raw_cov:
                # Initialize the coverage array with zeros
                raw_cov[a["ref"]] = np.zeros(ref_length[a["ref"]])
                final_cov[a["ref"]] = np.zeros(ref_length[a["ref"]])

            # Add to the raw coverage metric
            start = a["pos"] - 1
            end = a["pos"] - 1 + a["len"]
            raw_cov[a["ref"]][start:end] += 1

        # Get the most likely alignment for this query
        ml_alignments = error_model.calc_max_likelihood(query)

        # Skip any reads that cannot be deduplicated
        assert len(ml_alignments) > 0
        if len(ml_alignments) == 1:
            deduplicated_reads += 1

            a = ml_alignments[0]

            # Increment the read count
            final_reads[a["ref"]] += 1

            # Add to the final coverage metric
            start = a["pos"] - 1
            end = a["pos"] - 1 + a["len"]
            final_cov[a["ref"]][start:end] += 1

    # Now calculate summary statistics for each reference
    output = [{
        "id": ref,
        "length": ref_length[ref],
        "raw_reads": raw_read_count,
        "raw_coverage": (raw_cov[ref] > 0).mean(),
        "raw_depth": raw_cov[ref].mean(),
        "final_reads": final_reads[ref],
        "final_coverage": (final_cov[ref] > 0).mean(),
        "final_depth": final_cov[ref].mean(),
    }
        for ref, raw_read_count in raw_reads.items()
    ]

    return total_reads, deduplicated_reads, output


def parse_alignments_by_query(align_fp):
    """Function to parse the SAM alignments, grouping them by query ID."""
    # Keep all of the alignments for a given read here
    buff = []
    # Keep track of the query ID for the previous alignment
    last_query = None
    with open(align_fp, "rt") as f:
        for line in f:
            if line[0] == "@":
                continue
            line = line.rstrip("\n").split("\t")
            if len(line) < 10:
                continue
            ref = line[2]
            # Skip unaligned reads
            if ref == "*":
                continue
            # Get the query ID
            query = line[0]
            # Get the start position
            start_pos = int(line[3])

            # Get the alignment score
            n_mismatch = parse_mismatch_tag(line[11:])

            # If this is a new query ID
            if query != last_query:
                # If there is data for the last query, yield it
                if len(buff) > 0:
                    yield buff
                    # Reset the buffer
                    buff = []
                # The query sequence is only printed in the first row
                query_len = len(line[9])

                # Reset the last query sequence
                last_query = query

            # Add this alignment's information to the buffer
            buff.append({
                "query": query,
                "ref": ref,
                "pos": start_pos,
                "len": query_len,
                "n_mismatch": n_mismatch
            })

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
