#!/usr/bin/python

from famli_helpers import parse_alignment

n_reads, n_dedup, output = parse_alignment("/usr/famli/tests/example.sam")

assert n_reads == 319
assert n_dedup == 319

# Five references were aligned against
assert len(output) == 5

# Only two had deduplicated reads
assert len([a for a in output if a["final_reads"] > 0]) == 2

# The number of total and deduplicated reads adds up to the totals
n_reads = sum([a["raw_reads"] for a in output])
n_dedup = sum([a["final_reads"] for a in output])

print("PASSED TESTS")
