#!/bin/bash 

set -e 

docker run --rm \
-v $PWD:/share \
-v ~/.aws/credentials:/root/.aws/credentials \
famli:latest \
python3 \
/share/diamond-tax.py \
--input s3://fh-pi-fredricks-d/lab/Sam_Minot/dbs/test.tax.db.dmnd/test.tax.db.query.fastp \
--ref-db s3://fh-pi-fredricks-d/lab/Sam_Minot/dbs/test.tax.db.dmnd/test.tax.db.dmnd \
--output-folder /share/ \
--blocks 1 \
--threads 1 \
--top-pct 1 \
--sample-name TEST
