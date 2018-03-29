#!/usr/bin/python

from famli_helpers import FAMLI_Reassignment

famli = FAMLI_Reassignment()
with open("/usr/famli/tests/example.diamond.aln", "rt") as f:
    famli.parse(f)
    n_reads, output = famli.summary()

n_dedup = sum([d["nreads"] for d in output])

assert n_reads == 360, n_reads
assert n_dedup == 360, n_dedup

# Three references survived filtering
assert len(output) == 3

print("PASSED TESTS")
