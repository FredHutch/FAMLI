# FAMLI
Functional Analysis of Metagenomes by Likelihood Inference

[![Docker Repository on Quay](https://quay.io/repository/fhcrc-microbiome/famli/status "Docker Repository on Quay")](https://quay.io/repository/fhcrc-microbiome/famli)

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


  1. Map all input nucleotide reads in amino acid space against a reference database
  2. Use the user-supplied nucleotide error rate and the user-supplied genetic code to
  compute the expected rate of amino acid substitutions due to random nucleotide error
  3. Calculate the total likelihood mass for each reference, which is the sum of the likelihoods
  that each individual read is truly derived from that reference (given the relative alignment
  score, the amino acid substitution rate, the length of the sequence, and the binomial distribution)
  4. Assign each read to a single reference, using the total likelihood mass for that reference
  and the computed probability for the single alignment.
  5. Calculate summary statistics for the coverage and abundance of each reference using the
  deduplicated set of reads


Here are some examples:

  * If read #1 aligns equally well to reference A and reference B, but there is _2x more_ evidence
  for reference A across the entire sample, then read #1 gets **assigned to reference A**

  * If read #1 aligns equally well to reference A and reference B, but there is _equal_ evidence
  for reference A across the entire sample, then read #1 is **not assigned to any reference**
  (see assumption 2 above)

  * If read #1 aligns to reference A with 0 mismatches and reference B with 1 mismatch, but there
  is _1,000x more_ evidence for reference B across the entire sample, then read #1 gets
  **assigned to reference B**

  * If read #1 aligns to reference A with 0 mismatches and reference B with 1 mismatch, but there
  is _equal_ evidence for reference A and reference B across the entire sample, then read #1 gets
  **assigned to reference A**


### Math

#### Defining the effective amino acid substitution rate

*Terms:*
  * Amino acid error rate (Eaa)
  * Nucleotide error rate (Enuc)
  * Number of possible non-synonymous nucleotide errors (N-nonsyn)
  * Number of possible nucleotide errors (N-total)

**Eaa = Enuc * 3 * N-nonsyn / N-total**

#### Defining alignment likelihood

Consider a set of *i* nucleotide sequences being aligned against a set of *j* protein reference sequences.


*Edit distance (ED)*

The quality of each alignment can be quantified with an alignment score (AS) in which larger
values indicate a higher quality alignment. The edit distance (ED) for each alignment is calculated as 

ED*ij* = max(AS*i*) - AS*ij*

Where AS*i* is the set of all alignment scores for query *i* and AS*ij* is the alignment score for query
*i* against reference *j*.

*Likelihood (L)*

Let us consider two types of likelihood, the likelihood that a given query is truly from a given reference
*given only the evidence from that single query*, and the likelihood that a given query is truly from a 
given reference *given all of the evidence from all of the queries in a sample.* For the terms of this
discussion, we will describe the first (likelihood based on the evidence from a single read) as the *weight*
(W*ij*) for a given assignment, and the second (likelihood based on the evidence from the entire sample) as 
the *likelihood* (L*ij*) for a given assignment. 

To calculate the weight, we use a binomial distribution to calculate the probability that the observed number
of amino acid substitutions would occur due to random chance, given a query sequence of a given length and an 
expected amino acid error rate. Note that this is additive with the probability that a *larger* number of
substitutions would occur due to random sequencing error. That metric, based on the binomial distribution
is depicted here as B(LENGTH*i*, AS*ij*, Eaa)

The weight for each alignment is therefore calculated as:

W*ij* = B(LENGTH*i*, AS*ij*, Eaa)

Next, we calculate the total weight for every reference *j*

TOT*j* = sum(W*ij* for all *i*)

Finally, we calculate the likelihood that any individual query *i* is truly derived from a query *j*

L*ij* = W*ij* * TOT*j*

A query is assigned to a single reference sequence when the condition L*ij* == max(L*ij* for all *j*)
is satisfied for a single query *j*. 

NB: The only cases in which a query can not be assigned to a single reference is the case in which either
(a) the reference sequence does not contain any unique amino acid sequence, or (b) the sample has not been
sequenced to sufficient depth to detect any of that unique sequence.

### Implementation

**Aligner**: For alignment of nucleotide sequences against a protein database, we are currently using
(Paladin)[https://github.com/twestbrookunh/paladin], which is a similar algorithm to BWA-MEM and
attempts to align the entire read, rather than performing a local alignment (such as BLAST or DIAMOND).

**Alignment score**: The alignment score for each alignment is encoded in the BAM tag field as `AS:i:<INT>`.

**Error model**: The amino acid error model being used is the binomial distribution (`scipy.stats.binom`), which takes
into account both the length of the query sequence and the expected amino acid substitution rate.

### Wrapper Script

The entire process can be run as a single command with the script `famli.py`. That script encompasses:

  1. Downloading a reference database (if needed)
  2. Downloading the input data from SRA, AWS S3, or FTP (if needed)
  3. Aligning the input data (FASTQ) against the reference database
  4. Parsing the aligned reads
  5. Assigning the multi-mapped reads in the manner described above
  6. Computing coverage & depth metrics for each reference
  7. Writing the output as a single JSON file, either locally or to AWS S3

Example invocation of `famli.py` inside of the docker image for this repo (`famli:latest`):

```
docker run \
  -v $PWD/tests:/share \
  --rm \
  famli:latest \
    famli.py \
      --input /share/example.fastq \
      --ref-db /share/ \
      --output-folder /share/ \
      --temp-folder /share/ \
      --query-gencode 11 \
      --error-rate 0.001

```

Running the above command in the base directory for this repo should create an output file
(`example.fastq.json.gz`) in the `tests/` directory.


### Caveats

For the approach described above, we should note that there are situations in which the observed abundance
of highly abundant references will be inflated in a sample, in cases when there are more lowly abundant 
protein-coding sequences present that share a significant amount of homology to the dominant sequence. 
In other words, for two truly present references sharing a large region of exact amino acid identity, the 
multi-mapping reads from that redundant region will be entirely assigned to the dominant reference, instead
of being split between the two. That said, the less-abundant reference will still be detected in the output,
and all of the reads mapping to regions with unique amino acid sequences should still be assigned correctly.
