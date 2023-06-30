#!/bin/bash
echo "Running Test..."
#OUTPUT-FILE = ./outputs/AGAA_Output.txt
../../bin/RNACMB_debug -i inputfile.txt --parallel-lib-lengths 1-2-146-146-82 --parallel-lib-index 98-49-0-0-0 > outputs/runlog.dat
#../../bin/RNACMB_debug -i inputfile.txt 
