#!/usr/bin/python

from famli_helpers import parse_alignment

with open("/usr/famli/tests/example.diamond.aln", "rt") as f:
    n_reads, output = parse_alignment(f)

n_dedup = sum([d["nreads"] for d in output])

assert n_reads == 360, n_reads
assert n_dedup == 359, n_dedup

# Three references survived filtering
assert len(output) == 3

print("PASSED TESTS")
