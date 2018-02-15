#!/usr/bin/python

import os
import gzip
import json
import logging
import argparse
from collections import defaultdict
import numpy as np


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

def load_and_filter_alignment_stage_1(align_handle,
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
                                      STRIM_3=18):
    # Loads an alignment from file. Filters subjects by evenness. Groups the alignments into subjects with shared reads. 
    # Returns the groups of subject-alignments with shared reads, queries in each group, and the subject lengths for the subjects that passed filter. 

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

    logging.info("Done parsing {:,} subjects and {:,} queries".format(
        len(parser.subject_len),
        len(parser.query_set)
    ))

    # Sort the alignment list in place, by subject
    alignments.sort(key=lambda v: v[1])
    logging.info("Completed sorting the alignments by subject")
    # Pop through the sorted alignments, taking advantage of the sorting 
    # For each subject's block, build a could-coverage-o-gram 
    # Pass on alignments that meet criteria, grouping into blocks of subjects with shared reads

    # Output variables: A list of lists. Each entry is for a group of alignments for subject that share ANY queries. 
    group_alignments = []
    # This is a list of sets, where each set contains all of the queries for a given group. 
    group_queries = []
    # This is a dict of subject lengths only for those subjects that passed the filter
    subject_length = {}

    # Temporary variable to support the state-machine-looping strategy
    subject_aln = []
    cur_subject = None
    subject_could_coverage = np.zeros(shape=(1,0))
    aln_i =0
    start_aln_num = len(alignments)

    while len(alignments) > 0:
        aln_i +=1
        if aln_i %100000 == 0:
            logging.info("Alignment {:,} of {:,} considered in stage 1".format(aln_i, start_aln_num))        
        aln_line = alignments.pop()
        # New subject block?
        if cur_subject != aln_line[1]:
            if subject_could_coverage[STRIM_5:-STRIM_3].sum() > 0:
                # Only test if we have some alignments in our trimmed region
                if (np.std(subject_could_coverage[STRIM_5:-STRIM_3]) / subject_could_coverage[STRIM_5:-STRIM_3].mean()) < SD_MEAN_CUTOFF:
                    # Did it pass? If yes, add these alignments to the next stage, grouping here.

                    # Add this subject length to the subject_length dict
                    subject_length[cur_subject] = parser.subject_len[cur_subject]

                    # Get all the queries for this subject via set comprehension
                    subject_queries = {v[0] for v in subject_aln}
                    # Screen the existing group_queries for overlaps
                    group_query_overlaps = [len(subject_queries.intersection(gq))> 0 for gq in group_queries]
                    # Handle each case with the groups
                    if np.sum(group_query_overlaps) == 0: 
                        # No overlap with existing groups. Create a new one. 
                        group_queries.append(subject_queries)
                        group_alignments.append(subject_aln)
                    elif np.sum(group_query_overlaps) ==1: 
                        # Only ONE group overlaps. Easy. Just add to this group
                        # Get the index of the group with which we have an overlap
                        overlap_group_idx  =[i for i, o in enumerate(group_query_overlaps) if o][0]
                        group_queries[overlap_group_idx].update(subject_queries)
                        group_alignments[overlap_group_idx]+=subject_aln
                    else: 
                        # Implicitly there must be multiple overlapping groups
                        # Use list comprension to get the indicies of the overlapping groups
                        overlapping_group_idxs = [i for i, o in enumerate(group_query_overlaps) if o]
                        # Combine all of these groups into new combined groups of queries and alignments, starting with this subject's
                        combined_group_aln = subject_aln
                        combined_group_queries = subject_queries
                        # To adjust for the removed groups, we need to keep track of how many groups we've already popped.
                        num_popped=0
                        for g_i in sorted(overlapping_group_idxs):
                            # Use pop to remove the old groups from the lists
                            combined_group_aln+=group_alignments.pop(g_i-num_popped)
                            combined_group_queries.update(group_queries.pop(g_i-num_popped))
                            num_popped+=1
                        # Append the combined group to the end. 
                        group_queries.append(combined_group_queries)
                        group_alignments.append(combined_group_aln)    
            
            # Initialize for this new subject
            cur_subject = aln_line[1]
            subject_aln = []
            subject_could_coverage =  np.zeros(shape=(parser.subject_len[aln_line[1]],), dtype=int)
            # End we've found a new subject
        # Every loop through
        subject_aln.append(aln_line)
        subject_could_coverage[aln_line[2]:aln_line[3]] += 1
    


    return (group_alignments, group_queries, subject_length)

def iterative_alignment_prune_and_evenness_filter(group_aln, group_q, group_subjects, subject_length, STRIM_5=18, STRIM_3=18, ITERATIONS_MAX=1000, SD_MEAN_CUTOFF = 1.0):
    group_output = []
    logging.info("Stage 2: Iterative reassignments of {:,} subjects and {:,} queries".format(len(group_subjects),len(group_q)))
    # Indicies for quick lookup
    subject_i = {n: i for i, n in enumerate(group_subjects)}
    query_i = {n: i for i, n in enumerate(group_q)}

    # Subject Length vector to use for later normalization of the subject read depths
    subject_length_vector = np.zeros(shape=(len(subject_i)), dtype=int)
    for subject in subject_i:
        subject_length_vector[subject_i[subject]] = subject_length[subject]

    # Numpy matrix to store the bitscores
    bitscore_mat = np.zeros(
        shape=(len(subject_i), len(query_i)),
        dtype='float64'
    )

    logging.info("Filling in the bitscore matrix")
    for query, subject, sstart, send, bitscore in group_aln:
        bitscore_mat[subject_i[subject], query_i[query]] = bitscore

    logging.info("Completed filling the bitscore matrix.")

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
            bitscore_mat.T * (np.sum(bitscore_mat, axis=1)/subject_length_vector)
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
            logging.info("Iterations complete")
            break
        # Implicit else
        prior_align_norm_min_max = align_norm_min_max
        # and proceed to next iteration

    # Now let us get the subjects with some remaining alignments.
    subjects_with_iteratively_aligned_reads = {[p[0] for p in subject_i.items() if p[1] == i][0]: i for i in np.argwhere(align_mat_w.sum(axis=1)).flatten()}

    # Stage 3. Regenerate coverage-o-grams for the subjects that still have iteratively mapped queries. 
    # Create coverage-O-grams for the subjects with iteratively assigned reads
    # Use our same SD / mean coverage metric to screen now with the reassigned reads.
    logging.info("Starting stage 3 coverage evenness screening of {:,} subjects with filtered alignments".format(len(subjects_with_iteratively_aligned_reads)))
    group_subject_final_passed = set()

    for subject in subjects_with_iteratively_aligned_reads:
        # np zero vector of ints to be our coverage-o-gram
        subject_coverage = np.zeros(shape=(subject_length[subject],), dtype=int)

        # For subject s, which query_ids still have a non-zero normalized alignment score. 
        subject_queries = {[p[0] for p in query_i.items() if p[1] == i][0] for i in np.argwhere(align_mat_w[subjects_with_iteratively_aligned_reads[subject]]).flatten()}
        for q, s, sstart, send, bitscore in [a for a in group_aln if a[0] in subject_queries and a[1] == subject]:
            # Fill in the coverage-o-gram with these queries
            subject_coverage[sstart:send] += 1

        # Our testing point. With the trimmed coverage-o-gram are we still below our evenness threshold?
        if (np.std(subject_coverage[STRIM_5:-STRIM_3]) / subject_coverage[STRIM_5:-STRIM_3].mean()) < SD_MEAN_CUTOFF:
            group_output.append({
                "id": subject,
                "length": subject_length[subject],
                "depth": np.mean(subject_coverage),
                "coverage": np.mean(subject_coverage > 0),
                "nreads": len(subject_queries),
            })
            group_subject_final_passed.add(subject)

    logging.info("Done filtering {:,} subjects for which {:,} passed".format(len(subject_i), len(group_subject_final_passed)))
    
    return group_output

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
                    ITERATIONS_MAX=1000):
    """
    Parse an alignment in BLAST6 format and determine which subjects are likely to be present. This is the core of FAMLI.
    BLAST 6 columns by default (in order): qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
                                              0     1       2       3       4       5       6   7   8       9   10      11      12   13  
    """

    # Stage 1: Load and filter possible alignments by eveness of the subject coverage, and group by subjects with overlapping queries

    group_alignments, group_queries, subject_length = load_and_filter_alignment_stage_1(
                                                                                    align_handle,
                                                                                    QSEQID_i,
                                                                                    SSEQID_i,
                                                                                    QSTART_i,
                                                                                    QEND_i,
                                                                                    SSTART_i,
                                                                                    SEND_i,
                                                                                    BITSCORE_i,
                                                                                    SLEN_i,
                                                                                    SD_MEAN_CUTOFF,
                                                                                    STRIM_5,
                                                                                    STRIM_3
    )

    logging.info("{:,} groups of of queries with shared queries after filtering for evenness of possible coverage".format(len(group_alignments)))


    # Stage 2. PER GROUP Use the combination of the bitscores (alignment quality) and
    # subject-read-depth to iterative filter low-likely alignments of
    # queries against subjects

    output = []
    for i, (group_aln, group_q) in enumerate(zip(group_alignments, group_queries)):
        # Set comprehension to get all the possible subjects for this group
        group_subjects = {a[1] for a in group_aln}
        logging.info("Starting stage 2 and 3 for group {:,} of {:,}, containing {:,} subjects and {:,} queries".format(i+1, len(group_alignments), len(group_subjects),len(group_q)))
        output+=iterative_alignment_prune_and_evenness_filter(group_aln, group_q, group_subjects, subject_length, STRIM_5, STRIM_3, ITERATIONS_MAX, SD_MEAN_CUTOFF)
        logging.info("Completed for group {:,}. {:,} subjects passed so far".format(i+1,len(output)))
    logging.info("Completed for all groups. {:,} subjects passed".format(len(output)))
    

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
                        required=True,
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
            output = parse_alignment(f)
    else:
        with open(args.input, "rt") as f:
            output = parse_alignment(f)

    if args.output:
        with open(args.output, "wt") as fo:
            json.dump(output, fo)
