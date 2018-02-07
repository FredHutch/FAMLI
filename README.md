# FAMLI
Functional Analysis of Metagenomes by Likelihood Inference

Authors: 

  * Samuel Minot, Ph.D.
  * Jonathan Golob, M.D., Ph.D.


### Introduction

The goal of this work is to improve the accuracy of identifying protein-coding sequences
from short-read shotgun metagenomic sequencing data. The core challenge we consider here
is that of 'multi-mapping' reads – short nucleotide sequences that are equally similar to
two different reference protein sequences. In other domains such multi-mapping reads can
be dealt with in a variety of ways. For example, in the realm of taxonomic identification
it can be appropriate to assign them to the lowest common ancestor (LCA) of both references. 
However in the case of mapping short reads to a database of protein sequences we can not
assume that there is an underlying taxonomic structure to support such aggregation. Instead,
we can use the following simplifying assumptions to address this problem in a more general
manner. 


  1. All reference sequences contain some amount of unique amino acid sequence
  2. References can only be detected when they have been sequenced across such a unique region
  3. Errors in genome sequencing are distributed randomly


NB: The third assumption is more true for Illumina sequencing chemistry, and notably not
true for PacBio or MinION sequencing, which seem to have clustered or homopolymer-related
errors


### Approach

The approach that we use here is to reassign multi-mapping reads to the most likely reference
sequence that it aligns to, while taking into account the potential for mis-mapping due to
random sequencing error. Operationally, this consists of the following steps:


  1. Map all reads in amino acid space against a reference database
  2. Use the user-supplied nucleotide error rate and the user-supplied genetic code to
  compute the expected rate of amino acid substitutions due to random nucleotide error
  3. Calculate the total likelihood mass for each reference, which is the sum of the likelihoods
  that each individual read is truly derived from that reference (given the relative alignment
  score, the amino acid substitution rate, the length of the sequence, and the binomial distribution)
  4. Assign each read to a single reference, using the total likelihood mass for that reference
  and the computed probability for the single alignment.


Here are some examples:

  * If read #1 aligns equally well to reference A and reference B, but there is _2x more_ evidence
  for reference A across the entire sample, then read #1 gets **assigned to reference A**

  * If read #1 aligns equally well to reference A and reference B, but there is _equal_ evidence
  for reference A across the entire sample, then read #1 is **not assigned to any reference** (see assumption 2 above)

  * If read #1 aligns to reference A with 0 mismatches and reference B with 1 mismatch, but there is _1,000x more_ 
  evidence for reference B across the entire sample, then read #1 gets **assigned to reference B**

  * If read #1 aligns to reference A with 0 mismatches and reference B with 1 mismatch, but there is _equal_
  evidence for reference A and reference B across the entire sample, then read #1 gets **assigned to reference A**

