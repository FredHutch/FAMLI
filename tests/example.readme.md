Description of example datasets


1. Selected a set of 3 genes from the E. coli genome, including 30bp of flanking sequence.
The combined reference nucleotide data is saved as `example.ref.fasta`.

  * mreD (NC_000913.3:3397815-3398408)
  * ftsZ (NC_000913.3:105275-106486)
  * ftsQ (NC_000913.3:102203-103183)


2. Simulate a set of short sequences using the Illumina HiSeq2500 error profile using ART.
(`art_illumina -ss HS25 -i example.ref.fasta -l 150 -f 20 -o example`). The simulated reads
can be found in the file `example.fastq`.

3. Identify the UniRef100 accessions that perfectly match the three genes above.

  * UniRef100_P0A9A8
  * UniRef100_P25536
  * UniRef100_P07862

4. Collect the complete set of UniRef100 accessions that are found within the UniRef50 groups
that contain the three UniRef100 accessions listed above (matching the genes from the simulation).
That set of reference protein sequences can be found in `refdb.fastp`.

5. Index that set of reference proteins with Paladin:

```
docker run --rm -v $PWD:/share famli:latest paladin index -r3 refdb.fastp
```

5. Align the simulated FASTQ data against the reference UniRef100 groups with the command:

```
docker run --rm -v $PWD:/share famli:latest paladin align -a -z 11 -v 0 refdb.fastp example.fastq > example.sam
```
