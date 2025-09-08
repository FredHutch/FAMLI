#!/usr/bin/env bats

@test "FASTX Toolkit 0.0.13" {
  v="$(fastq_quality_trimmer -h 2>&1 || true )"
  [[ "$v" =~ "FASTX Toolkit 0.0.13" ]]
}


@test "AWS CLI" {
  v="$(aws --version 2>&1)"
  [[ "$v" =~ "aws-cli" ]]
}


@test "Curl v8.5.0" {
  v="$(curl --version)"
  [[ "$v" =~ "8.5.0" ]]
}

@test "DIAMOND v0.9.22" {
  v="$(diamond --version)"
  [[ "$v" =~ "0.9.22" ]]
}

@test "Make sure the run script is in the PATH" {
  h="$(famli -h 2>&1)"

  [[ "$h" =~ "FAMLI" ]]
}

@test "Parse the alignments in batches" {
  h="$(python /usr/famli/tests/test_alignment_batchsize.py)"

  [[ "$h" =~ "PASSED TESTS" ]]
}

@test "Parse the alignments" {
  h="$(python /usr/famli/tests/test_parse_alignments.py)"

  [[ "$h" =~ "PASSED TESTS" ]]
}

@test "Quality trimming" {
  h="$(python /usr/famli/tests/test_quality_trim.py)"

  [[ "$h" =~ "PASSED TESTS" ]]
}

@test "FAMLI integration tests" {
  h="$(python /usr/famli/tests/integration_tests.py)"

  [[ "$h" =~ "PASSED TESTS" ]]
}
