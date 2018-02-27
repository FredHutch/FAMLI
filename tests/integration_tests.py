#!/bin/python

import gzip
import json
import subprocess

# Run the entire process on the test data

p = subprocess.Popen([
    "famli.py",
    "--input", "/usr/famli/tests/example.fastq",
    "--ref-db", "/usr/famli/tests/refdb.dmnd",
    "--output-folder", "/usr/famli/tests",
    "--temp-folder", "/usr/famli/tests"
])
stdout, stderr = p.communicate()
exitcode = p.wait()

assert exitcode == 0

output = json.load(gzip.open("/usr/famli/tests/example.fastq.json.gz"))

assert output["aligned_reads"] == 338, output["aligned_reads"]
assert output["total_reads"] == 360, output["total_reads"]
assert output["ref_db"] == "/usr/famli/tests/refdb.dmnd", output["ref_db"]
assert output["deduplicated_reads"] == 154, output["deduplicated_reads"]
assert output["output_folder"] == "/usr/famli/tests", output["output_folder"]
assert output["input_path"] == "/usr/famli/tests/example.fastq", output["input_path"]

# Two references had deduplicated reads
assert len(output["results"]) == 2, len(output["results"])

# Run the same tests with quality trimming
p = subprocess.Popen([
    "famli.py",
    "--input", "/usr/famli/tests/example.fastq",
    "--ref-db", "/usr/famli/tests/refdb.dmnd",
    "--output-folder", "/usr/famli/tests",
    "--temp-folder", "/usr/famli/tests",
    "--min-qual", "35"
])
stdout, stderr = p.communicate()
exitcode = p.wait()

assert exitcode == 0

output = json.load(gzip.open("/usr/famli/tests/example.fastq.json.gz"))

assert output["aligned_reads"] == 338, output["aligned_reads"]
assert output["total_reads"] == 360, output["total_reads"]
assert output["ref_db"] == "/usr/famli/tests/refdb.dmnd", output["ref_db"]
assert output["deduplicated_reads"] == 253, output["deduplicated_reads"]
assert output["output_folder"] == "/usr/famli/tests", output["output_folder"]
assert output["input_path"] == "/usr/famli/tests/example.fastq", output["input_path"]

# THREE references now have deduplicated reads
assert len(output["results"]) == 3, len(output["results"])

print("PASSED TESTS")
