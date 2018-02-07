#!/bin/bash

set -e

test_image(){
	
	img_tag=$1

	[[ "${#img_tag}" == "0" ]] && echo "Please specify image" && return

	[[ "$(docker run $img_tag echo True)" != "True" ]] && echo "Tag not found ($img_tag)" && return

	docker run \
		-v $PWD/tests:/share \
		--rm \
		$img_tag \
			famli.py \
				--input /share/example.fastq \
				--ref-db /share/ \
				--output-folder /share/ \
				--temp-folder /share/
}

test_image $1

