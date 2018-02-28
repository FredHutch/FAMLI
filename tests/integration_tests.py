#!/bin/python

import gzip
import json
import subprocess

# Run the entire process on the test data

p = subprocess.Popen([
    "famli.py",
    "align",
    "--input", "/usr/famli/tests/example.fastq",
    "--sample-name", "example1",
    "--ref-db", "/usr/famli/tests/refdb.dmnd",
    "--output-folder", "/usr/famli/tests",
    "--temp-folder", "/usr/famli/tests"
])
stdout, stderr = p.communicate()
exitcode = p.wait()

assert exitcode == 0

output = json.load(gzip.open("/usr/famli/tests/example1.json.gz"))

assert output["aligned_reads"] == 338, output["aligned_reads"]
assert output["total_reads"] == 360, output["total_reads"]
assert output["ref_db"] == "/usr/famli/tests/refdb.dmnd", output["ref_db"]
assert output["deduplicated_reads"] == 154, output["deduplicated_reads"]
assert output["output_folder"] == "/usr/famli/tests", output["output_folder"]
assert output["input_path"] == "/usr/famli/tests/example.fastq", output["input_path"]

# Two references had deduplicated reads
assert len(output["results"]) == 2, len(output["results"])

# QUALITY TRIMMING

# Run the same tests with quality trimming
p = subprocess.Popen([
    "famli.py",
    "align",
    "--input", "/usr/famli/tests/example.fastq",
    "--sample-name", "example2",
    "--ref-db", "/usr/famli/tests/refdb.dmnd",
    "--output-folder", "/usr/famli/tests",
    "--temp-folder", "/usr/famli/tests",
    "--min-qual", "35"
])
stdout, stderr = p.communicate()
exitcode = p.wait()

assert exitcode == 0

output = json.load(gzip.open("/usr/famli/tests/example2.json.gz"))

assert output["aligned_reads"] == 338, output["aligned_reads"]
assert output["total_reads"] == 360, output["total_reads"]
assert output["ref_db"] == "/usr/famli/tests/refdb.dmnd", output["ref_db"]
assert output["deduplicated_reads"] == 253, output["deduplicated_reads"]
assert output["output_folder"] == "/usr/famli/tests", output["output_folder"]
assert output["input_path"] == "/usr/famli/tests/example.fastq", output["input_path"]

# THREE references now have deduplicated reads
assert len(output["results"]) == 3, len(output["results"])

# MULTIPLE INPUT FILES, COMBINED

input_fps = "/usr/famli/tests/example.fastq+/usr/famli/tests/example2.fastq"

# Run the same tests with quality trimming
p = subprocess.Popen([
    "famli.py",
    "align",
    "--input",
    input_fps,
    "--sample-name", "combined",
    "--ref-db", "/usr/famli/tests/refdb.dmnd",
    "--output-folder", "/usr/famli/tests",
    "--temp-folder", "/usr/famli/tests",
    "--min-qual", "35"
])
stdout, stderr = p.communicate()
exitcode = p.wait()

assert exitcode == 0

output = json.load(gzip.open("/usr/famli/tests/combined.json.gz"))

assert output["aligned_reads"] == 676, output["aligned_reads"]
assert output["total_reads"] == 720, output["total_reads"]
assert output["ref_db"] == "/usr/famli/tests/refdb.dmnd", output["ref_db"]
assert output["deduplicated_reads"] == 507, output["deduplicated_reads"]
assert output["output_folder"] == "/usr/famli/tests", output["output_folder"]
assert output["input_path"] == input_fps, output["input_path"]

# THREE references now have deduplicated reads
assert len(output["results"]) == 3, len(output["results"])

# Just run the 'filter' module
p = subprocess.Popen([
    "famli.py",
    "filter",
    "--input",
    "/usr/famli/tests/example.diamond.aln",
    "--output", "/usr/famli/tests/example3.json"
])
stdout, stderr = p.communicate()
exitcode = p.wait()

assert exitcode == 0

output = json.load(open("/usr/famli/tests/example3.json", "rt"))
n_dedup = sum([d["nreads"] for d in output])

assert n_dedup == 359, n_dedup

# Three references survived filtering
assert len(output) == 3

print("PASSED TESTS")
