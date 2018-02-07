#!/usr/bin/env bats

@test "SRA Toolkit v2.8.2" {
  v="$(fastq-dump --version)"
  [[ "$v" =~ "2.8.2" ]]
}


@test "AWS CLI v1.11.146" {
  v="$(aws --version 2>&1)"
  [[ "$v" =~ "1.11.146" ]]
}


@test "Curl v7.47.0" {
  v="$(curl --version)"
  [[ "$v" =~ "7.47.0" ]]
}

@test "Make sure Paladin is in the PATH" {
  h="$(paladin version 2>&1 || echo 'done')"

  [[ "$h" =~ "Version: 1.4.1" ]]
}

@test "Make sure the run script is in the PATH" {
  h="$(famli.py -h 2>&1)"

  [[ "$h" =~ "FAMLI" ]]
}

@test "Parse the alignments" {
  h="$(python /usr/famli/lib/test_parse_alignments.py)"
  echo $h

  [[ "$h" =~ "FAMLI" ]]
}