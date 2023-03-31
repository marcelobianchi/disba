#!/bin/bash

MODEL=$1

[ ! -f $MODEL ] && echo "No model found" && exit

tmp=$(mktemp -d )
HERE=$(pwd)

(
	cd $tmp/ || exit
	truncate -s 0 LOG
	
	ln -s $HERE/$MODEL ./
	ln -s $HERE/periods.txt ./

	# 1st - run sprep96
	sprep96 -M $MODEL -L -R -PARR periods.txt  2>&1 | tee -a LOG
	[ ! -f sdisp96.dat ] && echo "Failed sprep96" && exit

	# 2nd - run sdisp96
	sdisp96 -v 2>&1 | tee -a LOG
	[ ! -f sdisp96.lov ] && echo "Failed sdisp96 Love" && exit
	[ ! -f sdisp96.ray ] && echo "Failed sdisp96 Rayleigh" && exit

	# 3rd - run 
	slegn96 2>&1 | tee -a LOG
	[ ! -f slegn96.egn ] && echo "Failed slegn96 Love" && exit

	sregn96 2>&1 | tee -a LOG
	[ ! -f sregn96.egn ] && echo "Failed sregn96 Rayleigh" && exit

	sdpegn96 -R -U -C -ASC
	mv -iv SREGN.ASC $HERE/$(basename $MODEL .model96)-R.asc

	sdpegn96 -L -U -C -ASC
	mv -iv SLEGN.ASC $HERE/$(basename $MODEL .model96)-L.asc

)

rm -rfv $tmp
