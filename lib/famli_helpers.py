#!/usr/bin/python

import json
import logging
import pandas as pd
from math import ceil
from scipy.stats import binom
from Bio.Data import CodonTable
from .exec_helpers import run_cmds
from itertools import permutations

from collections import defaultdict
import numpy as np


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

        # Keep track of the abundance of all references
        self.ref_abund = {}

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

    def add_reference(self, ref):
        """Set the starting abundance for all references."""
        self.ref_abund[ref] = 1

    def add_prob_to_alignments(self, alignments, query_ix):
        """Add the relative likelihood metric to a set of alignments."""
        # This is to be run for all alignments BEFORE calculating
        # the maximum likelihood for any individual alignment
        # Because: this function adds to the total probability
        # mass for each reference, which is used to calculate the
        # maximum likelihood in the function below

        for a in alignments:
            a["likelihood"] = self.edit_dist_prob(int(a["mismatch"]), int(a["qlen"]))

            # Set the query name
            self.query_names[query_ix] = a["query"]

            # Add the likelihood for this query and reference
            self.psingle[query_ix][a["subject"]] = a["likelihood"]

            self.add_reference(a["subject"])

    def make_psingle_matrix(self):
        """Calculate the prob of each read being from each reference."""
        # Initialize with zeros (float64)
        self.matrix = np.zeros((len(self.psingle), len(self.ref_abund)))

        # Set the indices for the rows and columns
        self.row_loc = {
            query: ix for ix, query in enumerate(self.psingle.keys())
        }
        self.row_ix = {ix: query for query, ix in self.row_loc.items()}
        self.col_loc = {
            ref: ix for ix, ref in enumerate(self.ref_abund.keys())
        }
        self.col_ix = {ix: query for query, ix in self.col_loc.items()}

        # Add the values for each alignment
        for query, ref_prob in self.psingle.items():
            for ref, prob in ref_prob.items():
                self.matrix[
                    self.row_loc[query],
                    self.col_loc[ref]
                ] += prob

    def optimize_assignments(self, cutoff=0.05):
        """Iteratively optimize the assignments of reads."""

        # Keep track of the total number of alignments that we started with
        n_input = sum([len(v.values()) for v in self.psingle.values()])

        keep_iterating = True
        iterations = 0

        # Make a matrix reformating the information in self.psingle
        # into a new matrix called self.matrix
        self.make_psingle_matrix()

        # Initialize the reference abundance (using self.matrix)
        # as the sum of the likelihood for each single read
        self.ref_abund = self.matrix.sum(axis=0)

        while keep_iterating:
            iterations += 1

            # Make a ptotal matrix, combining self.ref_abund and self.matrix
            # This represents the likelihood of each query truly being from
            # each reference, given the alignment quality and reference abund
            ptotal = self.matrix * self.ref_abund
            # Normalize each row (query) to sum to 1
            ptotal = ptotal / ptotal.sum(axis=1, keepdims=True)

            logging.info("ITERATION {:,} ({:,} alignments remaining)".format(
                iterations, (ptotal > 0).sum()))
            keep_iterating = False

            # For each query (row)
            for ix in range(ptotal.shape[0]):
                # Get the column indices to remove
                value_to_remove = ptotal[ix, ptotal[ix, :] > 0].min()
                if value_to_remove > cutoff:
                    continue
                if value_to_remove == ptotal[ix, :].max():
                    continue

                # Zero out the likelihood for any alignments with this value
                self.matrix[ix, ptotal[ix, :] == value_to_remove] = 0
                keep_iterating = True

            # Recalculate the reference abundance
            self.ref_abund = ptotal.sum(axis=0)

        # Loop through the reads
        final_counts = defaultdict(int)
        self.assignments = {}
        for query, query_ix in self.row_loc.items():
            # If there is only a single alignment left, assign the read
            potential_refs = [
                ref
                for ref, ref_ix in self.col_loc.items()
                if self.matrix[query_ix, ref_ix] > 0
            ]
            if len(potential_refs) == 1:
                ref = potential_refs[0]

                final_counts[ref] += 1

                self.assignments[query] = ref

        final_counts = pd.Series(final_counts)

        logging.info("Iterations: {:,}".format(iterations))
        logging.info("Input: {:,}".format(n_input))
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


def parse_alignment(align_fp,
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

    # Nested dictionary. First key is subject_id. Second key is query_id. Value is a tuple of (SSTART, SEND) for this subject-query pair
    subject_query_coverage = defaultdict(dict)

    # Dictionary. Key is subject_id. Value is the SLEN
    subject_len = {}

    # Set of all possible query_ids
    query_set = set()

    # Fill our data structures by parsing our alignment.
    logging.info("Starting to parse the alignment for the first time.")
    i = 0
    with open(align_fp, "rt") as f:
        for line in f:
            i += 1
            if i % 1000000 == 0:
                logging.info("{} lines of alignment parsed".format(i))
            line_list = line.strip().split()
            query = line_list[QSEQID_i]
            query_set.add(query)
            subject = line_list[SSEQID_i]
            subject_query_coverage[subject][query] = (
                int(line_list[SSTART_i]),
                int(line_list[SEND_i])
            )
            subject_len[subject] = int(line_list[SLEN_i])

    logging.info("Completed parsing for coverage {} subjects and {} queries".format(len(subject_query_coverage), len(query_set)))

    # Now we create coverage-o-grams for each subject, considering how it COULD be covered. Reject those for which the coverage is uneven enough to be implausible
    subjects_below_cutoff = set()
    for subject in subject_len:
        # Get our first layer dict lookup done
        subject_query = subject_query_coverage[subject]
        # Create a zero array of ints the length of this subject
        subject_could_coverage = np.zeros(
            shape=(subject_len[subject],),
            dtype=int
        )
        # For each query, add 1 for the relevant positions
        for q in subject_query:
            start_stop = subject_query[q]
            subject_could_coverage[start_stop[0]:start_stop[1]] += 1
        # Trim the 5' and 3' ends of the subject. Due to overhangs, they tend to be less covered. 
        subject_could_coverage_trimmed = subject_could_coverage[STRIM_5:-STRIM_3]
        # Using the trimmed coverage-o-gram, if the STD / mean of coverage is LESS than a threshold (empirically set to 0.33), it is a subject to keep
        subject_std = np.std(subject_could_coverage_trimmed)
        subject_mean = np.mean(subject_could_coverage_trimmed)
        if (subject_std / subject_mean) <= SD_MEAN_CUTOFF:
            subjects_below_cutoff.add(subject)
    logging.info("Kept {} of {} total subjects after filtering for evenness of possible coverage".format(len(subjects_below_cutoff),len(subject_len)))

    # 2. Use the combination of the bitscores (alignment quality) and
    # subject-read-depth to iterative filter low-likely alignments of
    # queries against subjects
    # Indicies for quick lookup
    subject_i = {n: i for i, n in enumerate(subjects_below_cutoff)}
    query_i = {n: i for i, n in enumerate(query_set)}

    # Subject Length vector to use for later normalization of the subject read depths
    subject_lengths = np.zeros(shape=(len(subject_i)), dtype=int)
    for subject in subject_i:
        subject_lengths[subject_i[subject]] = subject_len[subject]

    # Free up memory for subject-query pairs already filtered
    for subject in subject_len:
        if subject in subjects_below_cutoff is False:
            del subject_query_coverage[subject]

    # Numpy matrix to store the bitscores
    bitscore_mat = np.zeros(
        shape=(len(subject_i), len(query_i)),
        dtype='float64'
    )

    logging.info("Starting to parse alignment file for the second (and last) time, to fill the bitscore matrix.")
    i = 0
    with open(align_fp, "rt") as f:
        for line in f:
            i += 1
            if i % 1000000 == 0:
                logging.info("{} lines of alignment parsed".format(i))
            line_list = line.strip().split()
            query = line_list[QSEQID_i]
            subject = line_list[SSEQID_i]
            if subject in subject_i:
                bitscore_mat[subject_i[subject], query_i[query]] = float(line_list[BITSCORE_i])

    logging.info("Completed parsing alignment to fill bitscore matrix (second and final pass)")

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
        subject_coverage = np.zeros(shape=(subject_len[subject],), dtype=int)

        # For subject s, which query_ids still have a non-zero normalized alignment score. 
        subject_queries = {[p[0] for p in query_i.items() if p[1] == i][0] for i in np.argwhere(align_mat_w[subjects_with_iteratively_aligned_reads[subject]]).flatten()}
        for sq in subject_queries:
            # Fill in the coverage-o-gram with these queries
            start_stop = subject_query_coverage[subject][sq]
            subject_coverage[start_stop[0]:start_stop[1]] += 1

        # Our testing point. With the trimmed coverage-o-gram are we still below our evenness threshold?
        subject_std = np.std(subject_coverage[STRIM_5:-STRIM_3])
        subject_mean = np.mean(subject_coverage[STRIM_5:-STRIM_3])
        if subject_std / subject_mean < SD_MEAN_CUTOFF:
            subject_final_passed.add(subject)

    logging.info("Done filtering subjects. {} subjects are felt to likely be present".format(len(subject_final_passed)))

    return subject_final_passed


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
