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

        # The total number of substitutions at each of the three positions
        n_subs = [0, 0, 0]
        # The number of non-synonymous substitutions
        n_nonsyn_subs = [0, 0, 0]
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
                mismatches = [c1 != c2 for c1, c2 in zip(codon1, codon2)]
                if sum(mismatches) == 1:
                    # Increment the counter of all nucleotide subs
                    n_subs[mismatches.index(True)] += 1
                    if aa1 != aa2:
                        # Increment the counter of all non-syn subs
                        n_nonsyn_subs[mismatches.index(True)] += 1

        # Calculate the non-synonymous rate at each codon position
        nonsyn_rate = [
            x / float(y)
            for x, y in zip(n_nonsyn_subs, n_subs)
        ]
        logging.info("The proportion of non-synonymous nucleotides")
        for ix, v in enumerate(nonsyn_rate):
            logging.info("Codon {}: {}".format(ix + 1, round(v, 5)))

        # Calculate the effective amino acid substitution rate
        self.aa_error_rate = 1 - np.prod([
            1 - (r * error_rate)
            for r in nonsyn_rate
        ])

        logging.info("Given a nucleotide error rate of {},".format(
            round(error_rate, 4)))
        logging.info("the amino acid error rate is {}.".format(
            round(self.aa_error_rate, 4)))

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

    def add_references(self, all_refs):
        """Set the starting abundance for all references."""
        ref_prop = 1. / len(all_refs)
        self.ref_abund = {
            ref: ref_prop
            for ref in all_refs
        }

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

            # Add the likelihood for this query and reference
            self.psingle[query_ix][a["ref"]] = a["likelihood"]

    def update_ref_abund(self):
        """Calculate the estimated abundance per reference."""
        self.ref_abund = defaultdict(float)
        for query, ref_abund in self.psingle.items():
            for ref, abund in ref_abund.items():
                self.ref_abund[ref] += abund

    def optimize_assignments(self, cutoff=0.5, max_prop=0.1):
        """Iteratively optimize the assignments of reads."""

        # Keep track of the total number of alignments that we started with
        n_input = sum([len(v.values()) for v in self.psingle.values()])

        keep_iterating = True
        n_removed = 0
        iterations = 0

        # Recalculate the reference abundance
        self.update_ref_abund()

        while keep_iterating:
            iterations += 1

            logging.info("ITERATION {:,} ({:,} removed)".format(
                iterations, n_removed))
            keep_iterating = False

            # Calculate the likelihood for every query,
            # given the alignment likelihood and the reference abundance
            for query, ref_prob in self.psingle.items():

                # Calculate the reference probability based on the total sample
                ptotal = pd.Series({
                    ref: prob * self.ref_abund[ref]
                    for ref, prob in ref_prob.items()
                }).sort_values()

                # Divide by the sum
                ptotal = ptotal / ptotal.sum()

                # Filter out alignments for which the aggregate probability is
                # below the cuttoff, and also don't remove more than 10% of the
                # alignments (rounding up to the nearest integer)

                # Maximum number of alignments to consider for removal
                max_to_consider = ceil(ptotal.shape[0] * max_prop)

                # Iterate through each reference for this query
                for ref, agg_prop in ptotal.cumsum()[:max_to_consider].items():
                    # If the aggregate probability is below the cutoff
                    if agg_prop <= cutoff:
                        # Remove the information for this alignment
                        del self.psingle[query][ref]
                        n_removed += 1
                        keep_iterating = True

            # Recalculate the reference abundance
            self.update_ref_abund()

        # Loop through the reads
        final_counts = defaultdict(int)
        self.assignments = {}
        for query, query_ref_abund in self.psingle.items():
            # If there is only a single alignment left, assign the read
            if len(query_ref_abund) == 1:
                for ref in query_ref_abund.keys():
                    final_counts[ref] += 1

                    self.assignments[query] = ref

        final_counts = pd.Series(final_counts)

        logging.info("Iterations: {:,}".format(iterations))
        logging.info("Input: {:,}".format(n_input))
        logging.info("Removed: {:,}".format(n_removed))
        logging.info("Duplicated: {:,}".format(
            len(self.psingle) - final_counts.sum()))
        logging.info("Deduplicated: {:,}".format(final_counts.sum()))
        logging.info("Unique references: {}".format(final_counts.shape[0]))
        for ref, count in final_counts.items():
            logging.info("{}: {:,}".format(ref, count))


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


def parse_alignment(align_fp, error_rate=0.001, genetic_code=11):
    """Parse an alignment in SAM format."""

    # Initialize the error model
    error_model = ErrorModelFAMLI(error_rate, genetic_code=genetic_code)

    # Read over the file once to get all of the alignments
    # This will group the alignments by query sequence
    alignments = [a for a in parse_alignments_by_query(align_fp)]

    # Keep track of the total set of reference with any alignments
    all_refs = set([
        a["ref"]
        for query in alignments
        for a in query
    ])

    # Add the references to the error model
    error_model.add_references(list(all_refs))

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

    # Add a "likelihood" to every alignment
    # This is the likelihood that a single alignment is due to
    # sequencing error
    for query_ix, a in enumerate(alignments):
        error_model.add_prob_to_alignments(a, query_ix)

    # Optimize the read assignments, iteratively calculating the
    # total reference abundance and removing potential alignments
    # that fall below the specified level of likelihood
    error_model.optimize_assignments(cutoff=0.5, max_prop=0.1)

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
    for query_ix, query in enumerate(alignments):
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

        # Get the deduplicated reference for this query
        ref = error_model.assignments.get(query_ix)

        # Skip any reads that cannot be deduplicated
        if ref is None:
            continue

        deduplicated_reads += 1

        a = [a for a in query if a["ref"] == ref]

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
