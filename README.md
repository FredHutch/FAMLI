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
multiple different reference protein sequences. In other domains such multi-mapping reads can
be dealt with in a variety of ways. For example, in the realm of taxonomic identification
it can be appropriate to assign them to the lowest common ancestor (LCA) of both references. 

However in the case of mapping short reads to a database of protein sequences (or peptides) we can not
assume that there is an underlying directed acyclic graph structure (e.g. a taxonomy). Peptides
can evolve by duplication events, homologous recombination, and other means of sharing highly conserved
domains (leading to shared short reads). If one simply includes all peptides for which there is a read,
we find the false positives outnumber the true positive by as much as 1000:1. 

We developed a method to iteratively assign shared reads to the most likely true peptides, bringing the 
precision (TP / (TP+FP)) to close to 90%. To do so, we used the following principles:


  1. In peptides that are truly positive in the sample, there should be relatively even sequence 
  coverage across the length of the peptide. 
 
 Present:
 
```
  C:23445432
    ||||||||
  P:--------
```

Not present, but with a shared domain with a peptide that is present:
```
  C:23445432000000000000
    ||||||||
  P:--------------------
```

  2. We use the total depth of coverage for a peptide (normalized to the peptide length) to 
  iteratively reassign multiply aligned sequences to the more likely peptide to be present

### Approach

  1. Align all input nucleotide reads in amino acid space against a reference database of peptides.
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

  * For reference A and reference B that both have some aligning query reads, if **there is _uneven_ depth for reference A** 
  but relatively even depth across reference B, then **reference A is removed from the candidate list** while reference B 
  is kept as a candidate.

  * If **read #1 aligns equally-well to reference A and reference C**, but **there is _2x more_ read depth for reference A as 
  compared to reference C** across the entire sample, then **reference C's alignment is removed from the list of candidates 
  for read #1**.


### Math

#### Coverage Evenness
This is considered on a per-reference basis. On a per-amino-acid basis, alignment-depth is calculated using an integer vector. 
It is expected that the 5' and 3' ends of the reference will have trail offs, thus the vector is trimmed on both the 5' and 3' 
ends. A mean coverage depth and the standard deviation of the mean are calculated. The standard deviation is divided by the 
mean. Both based on the Poisson distribution and some empirical efforts on our part, we set a threshold of 1.0 for this ratio 
as a cutoff of uneveness; **references with a coverage SD / MEAN ratio > 1.0 are filtered**. 

#### Defining alignment likelihood

Let us consider the **likelihood that a given query i is truly from a given reference j** considering all of the evidence from all of the queries in a sample. For the terms of this discussion, we will describe this as the **likelihood** (L<sub>ij</sub>) for a given assignment. 

For our application here we use the **bitscore**--an integrated consideration of the alignment length, number of mismatches, gaps, and overhangs--as a way of comparing alignment quality for weighting: Bitscore<sub>ij</sub> is the quality of the alignment of query read *i* to reference *j*.

W can use the bitscore of an alignment divided by the sum of bitscores for all the alignments for a given query sequence as a **normalized weight** W<sub>ij</sub>. 

>W<sub>ij</sub> = Bitscore<sub>ij</sub> / Sum(Bitscore<sub>ij</sub> for all *j*) 

Next, we calculate the **total weight** for every reference *j*, **TOT<sub>j</sub>**

>TOT<sub>j</sub> = sum(W<sub>ij</sub> for all *i*)

Finally, we calculate the **likelihood** that any individual query *i* is truly derived from a reference *j*, **L<sub>ij</sub>**

>L<sub>ij</sub> = W<sub>ij</sub> * TOT<sub>j*

The **maximum likelihood for query i, Lmax<sub>i</sub>** is determined 
>Lmax<sub>i</sub> = max(L<sub>ij</sub> for all *j*).

If the L<sub>ij</sub> falls below the scaled maximum likelihood for query *i*, the **alignment is removed from consideration**:

>For all query *i*, 
>if L<sub>ij</sub> < scale * Lmax<sub>i</sub>, 
>then Bitscore<sub>ij</sub> is set to zero.


By default the scale here is set to 0.9 (or 90% of the maximum likelihood for query *i*).

This process (recalculate W<sub>ij</sub>, calculate the TOT<sub>j</sub> for each refrence *j*, and then calculate a 
L<sub>ij</sub> using the new W<sub>ij</sub> and TOT<sub>j</sub>) is **repeated iteratively until no more alignments 
are culled** or a maximum number of iterations is reached. 


### Implementation

**Aligner**: For alignment of nucleotide sequences against a protein database, we are currently using
DIAMOND [https://github.com/bbuchfink/diamond]. We specifically ran DIAMOND with the following alignment options:
```
--query-cover 90
--min-score 20
--top 10
--id 80
```

**Alignment score**: We use bitscores as calculated by DIAMOND as an integrated assessment of alignment quality 
(considering alignment length, gaps, mismatches, and query sequence quality).


### Wrapper Script

The entire process can be run as a single command with the script `famli`. That script encompasses:

  1. Downloading a reference database (if needed)
  2. Downloading the input data from SRA, AWS S3, or FTP (if needed)
  3. Aligning the input data (FASTQ) against the reference database
  4. Parsing the aligned reads
  5. Assigning the multi-mapped reads in the manner described above
  6. Computing coverage & depth metrics for each reference
  7. Writing the output as a single JSON file, either locally or to AWS S3

The script `famil` can run two commands: `filter` and `align`. 

#### `filter`

`filter` runs the core algorithm of FAMLI, processing a set of BLASTX-like alignments (in tabular format), 
filtering out unlikely proteins, and assigning multi-mapped reads to individual references. 

The options available when invoking `famli filter` are as follows:

```
usage: famli [-h] [--input INPUT] [--output OUTPUT] [--threads THREADS]
             [--logfile LOGFILE] [--qseqid-ix QSEQID_IX]
             [--sseqid-ix SSEQID_IX] [--qstart-ix QSTART_IX]
             [--qend-ix QEND_IX] [--sstart-ix SSTART_IX]
             [--send-ix SEND_IX] [--bitscore-ix BITSCORE_IX]
             [--slen-ix SLEN_IX] [--sd-mean-cutoff SD_MEAN_CUTOFF]
             [--strim-5 STRIM_5] [--strim-3 STRIM_3]

Filter a set of existing alignments in tabular format with FAMLI

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT         Location for input alignement file.
  --output OUTPUT       Location for output JSON file.
  --threads THREADS     Number of processors to use.
  --logfile LOGFILE     (Optional) Write log to this file.
  --qseqid-ix QSEQID_IX
                        Alignment column for query sequence ID. (0-indexed
                        column ix)
  --sseqid-ix SSEQID_IX
                        Alignment column for subject sequence ID. (0-indexed
                        column ix)
  --qstart-ix QSTART_IX
                        Alignment column for query start position. (0-indexed
                        column ix, 1-indexed start position)
  --qend-ix QEND_IX     Alignment column for query end position. (0-indexed
                        column ix, 1-indexed end position)
  --sstart-ix SSTART_IX
                        Alignment column for subject start position.
                        (0-indexed column ix, 1-indexed start position)
  --send-ix SEND_IX     Alignment column for subject end position. (0-indexed
                        column ix, 1-indexed end position)
  --bitscore-ix BITSCORE_IX
                        Alignment column for alignment bitscore. (0-indexed
                        column ix)
  --slen-ix SLEN_IX     Alignment column for subject length. (0-indexed column
                        ix)
  --sd-mean-cutoff SD_MEAN_CUTOFF
                        Threshold for filtering max SD / MEAN
  --strim-5 STRIM_5     Amount to trim from 5' end of subject
  --strim-3 STRIM_3     Amount to trim from 3' end of subject
```

#### `align`

`align` is used to process a set of nucleotide sequences in FASTQ format, aligning against a 
reference database using DIAMOND, processing the alignments, filtering out unlikely proteins,
and assigning multi-mapped reads to individual references. 

The options available when invoking `famli align` are as follows:

```
usage: famli [-h] --input INPUT --sample-name SAMPLE_NAME --ref-db REF_DB
              --output-folder OUTPUT_FOLDER [--min-score MIN_SCORE]
              [--blocks BLOCKS] [--query-gencode QUERY_GENCODE]
              [--threads THREADS] [--min-qual MIN_QUAL]
              [--temp-folder TEMP_FOLDER]

Align a set of reads with DIAMOND, filter alignments with FAMLI, and return
the results

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT         Location for input file(s). Combine multiple files
                        with +. (Supported: sra://, s3://, or ftp://).
  --sample-name SAMPLE_NAME
                        Name of sample, sets output filename.
  --ref-db REF_DB       Folder containing reference database. (Supported:
                        s3://, ftp://, or local path).
  --output-folder OUTPUT_FOLDER
                        Folder to place results. (Supported: s3://, or local
                        path).
  --min-score MIN_SCORE
                        Minimum alignment score to report.
  --blocks BLOCKS       Number of blocks used when aligning. Value relates to
                        the amount of memory used. Roughly 6Gb RAM used by
                        DIAMOND per block.
  --query-gencode QUERY_GENCODE
                        Genetic code used to translate nucleotides.
  --threads THREADS     Number of threads to use aligning.
  --min-qual MIN_QUAL   Trim reads to a minimum Q score.
  --temp-folder TEMP_FOLDER
                        Folder used for temporary files.
```

#### Installation

The software can be installed with the command `pip install famli`. If you choose to install
FAMLI in this way, you will also need to install a working copy of the DIAMOND aligner,
accessible at runtime via `diamond`. That is the only dependency not installed by `pip install famli`.
FAMLI has not been tested with Python3  -- it is only advisable to run it with Python2. 

#### Docker

Example invocation of `famli` inside of the docker image for this repo (`famli:latest`):

```
docker run \
  -v $PWD/tests:/share \
  --rm \
  famli:latest \
    famli \
    align \
      --input /share/example.fastq \
      --sample-name example \
      --ref-db /share/refdb.dmnd \
      --output-folder /share/ \
      --blocks 5 \
      --query-gencode 11 \
      --threads 16 \
      --min-qual 30 \
      --temp-folder /share/

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
