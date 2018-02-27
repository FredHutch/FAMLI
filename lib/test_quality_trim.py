#!/usr/bin/python

import os
from fastq_helpers import quality_trim

os.mkdir("/usr/famli/tests/quality_trim")

for min_qual, nchar in [[30, 120359], [20, 120439]]:

    trimmed_fp = quality_trim(
        "/usr/famli/tests/example.fastq",
        "/usr/famli/tests/quality_trim",
        min_qual=min_qual
    )

    assert trimmed_fp == "/usr/famli/tests/quality_trim/example.fastq"

    trimmed_reads = ''.join(open(trimmed_fp, "rt").readlines())
    assert len(trimmed_reads) == nchar

    os.remove(trimmed_fp)

print("PASSED TESTS")
