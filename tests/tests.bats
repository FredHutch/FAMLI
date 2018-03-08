#!/usr/bin/env bats

@test "SRA Toolkit v2.8.2" {
  v="$(fastq-dump --version)"
  [[ "$v" =~ "2.8.2" ]]
}

@test "FASTX Toolkit 0.0.13" {
  v="$(fastq_quality_trimmer -h 2>&1 || true )"
  [[ "$v" =~ "FASTX Toolkit 0.0.13" ]]
}


@test "AWS CLI v1.11.146" {
  v="$(aws --version 2>&1)"
  [[ "$v" =~ "1.11.146" ]]
}


@test "Curl v7.47.0" {
  v="$(curl --version)"
  [[ "$v" =~ "7.47.0" ]]
}

@test "DIAMOND v0.9.10" {
  v="$(diamond --version)"
  [[ "$v" =~ "0.9.10" ]]
}

@test "Make sure the run script is in the PATH" {
  h="$(famli.py -h 2>&1)"

  [[ "$h" =~ "FAMLI" ]]
}

@test "Parse the alignments in batches" {
  h="$(python /usr/famli/famli/test_alignment_batchsize.py)"

  [[ "$h" =~ "PASSED TESTS" ]]
}

@test "Parse the alignments" {
  h="$(python /usr/famli/famli/test_parse_alignments.py)"

  [[ "$h" =~ "PASSED TESTS" ]]
}

@test "Quality trimming" {
  h="$(python /usr/famli/famli/test_quality_trim.py)"

  [[ "$h" =~ "PASSED TESTS" ]]
}

@test "FAMLI integration tests" {
  h="$(python /usr/famli/tests/integration_tests.py)"

  [[ "$h" =~ "PASSED TESTS" ]]
}
