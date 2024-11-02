#!/bin/bash

grep -oP '(?<=\#).*(?=at)' $1 > chain_number.dat
grep -oP '(?<=at ).*(?=rms)' $1 > percent_iteration.dat
grep -oP '(?<=\[).*(?=\]\[)' $1 > model_rms.dat
grep -oP '(?<=accepted:).*(?=\%)' $1 > acceptance_ratio.dat
