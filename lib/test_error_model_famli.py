#!/bin/bash

import json
from famli_helpers import ErrorModelFAMLI

error_model = ErrorModelFAMLI(0.001, genetic_code=11)

# Sanity check -- given an overall nucleotide error rate of 0.1%, the rate
# of amino acid changes (per codon) should be between 0.1% and 0.5%
# (aa error = nuc error * 3 * prop_non_synonymous)
assert error_model.aa_error_rate < 0.005
assert error_model.aa_error_rate > 0.001

# In a peptide of 50aa, the likelihood of AT LEAST 0 substitions is 100%
assert round(error_model.edit_dist_prob(0, 50), 2) == 1.0

# In a peptide of 50aa, the likelihood of AT LEAST 1 substitions is ~ 10%
assert round(error_model.edit_dist_prob(1, 50), 1) == 0.1

# In a peptide of 50aa, the likelihood of AT LEAST 2 substitions is ~ 0.6%
assert round(error_model.edit_dist_prob(2, 50), 3) == 0.006


# Test whether likelihood values are being added to each alignment
alignments = [
    {
        "ref": "ref_a",
        "score": 50,
        "len": 50,
    },
    {
        "ref": "ref_b",
        "score": 50,
        "len": 50,
    },
    {
        "ref": "ref_c",
        "score": 49,
        "len": 50,
    },
    {
        "ref": "ref_d",
        "score": 48,
        "len": 50,
    },
]

error_model.add_prob_to_alignments(alignments)
for a in alignments:
    assert "likelihood" in a

# The likelihood for the first two alignments should be equal, and the highest
assert alignments[0]["likelihood"] == alignments[1]["likelihood"]
assert alignments[0]["likelihood"] > alignments[2]["likelihood"]
assert alignments[0]["likelihood"] > alignments[3]["likelihood"]
assert alignments[2]["likelihood"] > alignments[3]["likelihood"]

# Now add alignments for a second read
error_model.add_prob_to_alignments([
    {
        "ref": "ref_a",
        "score": 50,
        "len": 50,
    },
])

# Now that the highest mass is for ref_a, an ambiguous read
# aligning equally well to ref_a and ref_b should be assigned to ref_a
ambiguous_alignment = [
    {
        "ref": "ref_a",
        "score": 50,
        "len": 50,
    },
    {
        "ref": "ref_b",
        "score": 50,
        "len": 50,
    },
]
error_model.add_prob_to_alignments(ambiguous_alignment)
final_pick = error_model.calc_max_likelihood(ambiguous_alignment)
assert len(final_pick) == 1
assert final_pick[0]["ref"] == "ref_a"

print("PASSED TESTS")
