#!/bin/python

import gzip
import json
import subprocess

# Run the entire process on the test data

p = subprocess.Popen([
    "famli.py",
    "--input", "/usr/famli/tests/example.fastq",
    "--ref-db", "/usr/famli/tests/",
    "--output-folder", "/usr/famli/tests",
    "--temp-folder", "/usr/famli/tests"
])
stdout, stderr = p.communicate()
exitcode = p.wait()

assert exitcode == 0

output = json.load(gzip.open("/usr/famli/tests/example.fastq.json.gz"))

assert output["aligned_reads"] == 319
assert output["total_reads"] == 2270
assert output["ref_db"] == "/usr/famli/tests/refdb.fastp"
assert output["deduplicated_reads"] == 319
assert output["output_folder"] == "/usr/famli/tests", output["output_folder"]
assert output["input_path"] == "/usr/famli/tests/example.fastq"

# Five references were aligned against
assert len(output["results"]) == 5

# Only two had deduplicated reads
assert len([a for a in output["results"] if a["final_reads"] > 0]) == 2

# The number of total and deduplicated reads adds up to the totals
n_reads = sum([a["raw_reads"] for a in output["results"]])
n_dedup = sum([a["final_reads"] for a in output["results"]])

print("PASSED TESTS")
