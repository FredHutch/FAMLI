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
assume that there is an underlying directed acyclic graph structure. Peptides
can evolve by duplication events, homologous recombination, and other means of sharing highly conserved
domains (leading to shared short reads). If one simply includes all peptides for which there is a read, we find the false positives outnumber the true positive by as much as 1000:1. 

We developed a method to iteratively assign shared reads to the most likely true peptides, bringing the precision (TP / (TP+FP)) to close to 90%. To do so, we used the following principles:


  1. In peptides that are truly positive in the sample, there should be relatively even sequence 
  coverage across the length of the peptide. 
Present:
  C:23445432
    ||||||||
  P:--------

Not present, but with a shared domain with a peptide that is present:
  C:23445432000000000000
    ||||||||
  P:--------------------

  2. We can use the total depth of coverage for a peptide (normalized to the peptide length) to 
  iteratively reassign multiply aligned sequences to the more likely peptide to be present

### Approach

  1. Align all input nucleotide reads in amino acid space against a reference database.
  2. Filter out all recruited reference sequences with highly uneven coverage (assuming all
  possible aligning sequences are truly from this peptide):
  Standard deviation / Mean of coverage depth per amino acid of the peptide > 1.0
  3. Iteratively, until no further references are pruned: 
  i) NORMALIZATION: For sequences that align to multiple possible reference sequences, weight the alignment quality
  (bitscore) by the length-normalized total alignment quality for each candidate reference peptide. 
  ii) PRUNING. Remove from the candidate reference sequences for this sequence all references with 
  weighted alignment scores less than 90% of the maximum for this sequence. 
  4. Filter out all recruited reference sequences with highly uneven coverage after pruning of references in step 3.

Here are some examples:

  * If reads align to reference A and reference B, but there is _uneven_ depth for reference A but relatively even depth across reference B, then reference A is **removed from the candidate list**

  * If read #1 aligns equally-well to reference A and reference C, but there is _2x more_ length normalized read depth for reference A as compared to reference B across the entire sample, then **reference C's alignment is removed from the list of candidates for read #1**.


### Math

#### Coverage Evenness
This is considered on a per-reference basis. On a per-amino-acid basis, alignment-depth is calculated using an integer vector. It is expected that the 5' and 3' ends of the reference will have trail offs, thus the vector is trimmed on both the 5' and 3' ends. A mean coverage depth and the standard deviation of the mean are calculated. The standard deviation is divided by the mean. Both based on the Poisson distribution and some empirical efforts on our part, we set a threshold of 1.0 for this ratio as a cutoff of uneveness; references with SD / MEAN ratio > 1.0 are filtered. 

#### Defining alignment likelihood


*Likelihood (L)*

Let us consider the likelihood that a given query *i* is truly from a given reference *j* *given all of the evidence from all of the queries in a sample.* For the terms of this discussion, we will describe this as the *likelihood* (L*ij*) for a given assignment. 

Each alignment has a quality score; for our application here we use the bitscore--an integrated consideration of the alignment length, number of mismatches, gaps, and overhangs. 

Therefore, we can use the bitscore of an alignment divided by the sum of bitscores for all the alignments for a given query sequence as a normalized weight:

W*ij* = Bitscore*i* / Sum(Bitscore*ij* for all *j*) 

Next, we calculate the total weight for every reference *j*

TOT*j* = sum(W*ij* for all *i*)

Finally, we calculate the likelihood that any individual query *i* is truly derived from a query *j*

L*ij* = W*ij* * TOT*j*

The maximum likelihood for query i is determined Lmax*i* = max(L*ij* for all *j*). 

If the L*ij* falls below the scaled maximum likelihood for query i, the alignment is removed from consideration:

For all i, 
if L*ij* < scale*Lmax*i*, Bitscore*ij* 
then L*ij* is set to zero.

By default the scale here is set to 0.9 (or 90% of the maximum likelihood for query i).

This process (recalculate W*ij*, calculate the TOT*j* for each refrence *j*, and then calculate a L*ij* using the new W*ij* and TOT*j*) is repeated iteratively until no more alignments are culled or a maximum number of iterations is reached. 


### Implementation

**Aligner**: For alignment of nucleotide sequences against a protein database, we are currently using
(DIAMOND)[https://github.com/bbuchfink/diamond]. We specifically ran DIAMOND with the following alignment options:
--query-cover 90
--min-score 20
--top 10
--id 80

**Alignment score**: We use bitscores as calculated by DIAMOND as an integrated assessment of alignment quality (considering alignment length, gaps, mismatches, and query sequence quality).


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
