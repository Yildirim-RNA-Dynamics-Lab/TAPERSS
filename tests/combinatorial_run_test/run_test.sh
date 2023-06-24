#!/bin/bash
echo "Running Test..."
#OUTPUT-FILE = ./outputs/AGAA_Output.txt
../../bin/RNACMB_debug -i inputfile.txt --parallel-lib-lengths 33-127-146 --parallel-lib-index 0-0-0 -o outputs/AGAA_Output_01.txt 
../../bin/RNACMB_debug -i inputfile.txt --parallel-lib-lengths 33-127-146 --parallel-lib-index 1-0-0 -o outputs/AGAA_Output_02.txt 
../../bin/RNACMB_debug -i inputfile.txt --parallel-lib-lengths 33-127-146 --parallel-lib-index 2-0-0 -o outputs/AGAA_Output_03.txt 
../../bin/RNACMB_debug -i inputfile.txt --parallel-lib-lengths 33-127-146 --parallel-lib-index 3-0-0 -o outputs/AGAA_Output_04.txt 
#../../bin/RNACMB_debug -i inputfile.txt 
