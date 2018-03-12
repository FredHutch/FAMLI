#!/usr/bin/env bats

@test "Aspera Connect" {
  v="$(ascp --version)"
  [[ "$v" =~ "3.7.4.147727" ]]
}

@test "Aspera Connect Key File" {
  [[ -s $ASPERA_KEY ]]
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
