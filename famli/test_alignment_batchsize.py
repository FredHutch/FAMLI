#!/usr/bin/python

from famli_helpers import BLAST6Parser

parser = BLAST6Parser()

# Read the alignments as a single batch
with open("/usr/famli/tests/example.diamond.aln", "rt") as f:
    alignments = [a for a in parser.parse(f)]

assert len(alignments) == 100468, len(alignments)

# Now read the alignments in smaller batches, and make sure
# all of them make it through
for batchsize in [10, 20, 50, 100, 1e7]:
    with open("/usr/famli/tests/example.diamond.aln", "rt") as f:
        chunks = [a for a in parser.yield_alignments(f, batchsize=batchsize)]
    assert sum(map(len, chunks)) == len(alignments), chunks

print("PASSED TESTS")
