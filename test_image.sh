#!/bin/bash

set -e

test_image(){
	
	img_tag=$1

	[[ "${#img_tag}" == "0" ]] && echo "Please specify image" && return

	[[ "$(docker run $img_tag echo True)" != "True" ]] && echo "Tag not found ($img_tag)" && return

	echo "TESTING SMALL LOCAL FILE"

	docker run \
		-v $PWD/tests:/share \
		--rm \
		$img_tag \
			famli.py \
				align \
				--input /share/example.fastq \
				--sample-name example \
				--ref-db /share/refdb.dmnd \
				--output-folder /share/ \
				--temp-folder /share/ \
				--batchsize 50000


	echo -e "\n\nTESTING SMALL SRA FILE\n\n"

	docker run \
		-v $PWD/tests:/share \
		--rm \
		$img_tag \
			famli.py \
				align \
				--input sra://SRR6818184 \
				--sample-name SRR6818184 \
				--ref-db /share/refdb.dmnd \
				--output-folder /share/ \
				--temp-folder /share/ \
				--batchsize 50000

}

test_image $1

