#!/bin/bash

GRANAFILES=("sca" "hb" "def" "alias" "gr_wtdc.hst" "las_wtdc.hst" "hist.def")
for link in "${GRANAFILES[@]}"; do
	echo "ln -s GRAnalyzer/$link $link"
	ln -s GRAnalyzer/$link $link
done
