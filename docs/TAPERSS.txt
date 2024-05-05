#####################################################################################################
#                                        TAPERSS                                                    #
#         (Theoretical Analyses, Prediction, and Evaluation of RNA Structures from Sequence)        #
#                                                                                                   #
#               Combinatorial Method for RNA structure Prediction. Written in C++.                  #
#             Developed by Ivan Riveros and Ilyas Yildirim in part using initial code               #
#                              from Dr. Kye Won Wang as a guide.                                    #
#                                                                                                   #
# To cite the code:                                                                                 #
# Riveros, I. I. and Yildirim, I., "Prediction of 3D RNA Structures from Sequence Using Energy      #
# Landscapes of RNA Dimers: Application to RNA Tetraloops", J. Chem. Theory Comput.,                #
# doi.org/10.1021/acs.jctc.4c00189, 2024.                                                           #
#####################################################################################################
# For potential issues and solutions, see the bottom of this document.

# Note:
Run time varies heavily based on sequence, overlap RMSD cutoff, and hardware. Lower RMSD cutoff values 
will result in few structures and shorter runtimes. Based on tests, a 6-long sequence with an overlap 
cutoff of 0.5 running on an Intel i7-8700K CPU takes approximately 2-4 hours to complete.
######################################################################################################
# To compile the code:
From the ./TAPERSS directory run 'make'.

You must specify the compiler to use by creating an environment variable "GCC", see example below. This 
code has only been tested with clang and g++. Compilation with clang has shown slightly better performance. 
Ensure the compiler used is up to date with the latest version, or at least supports the "std=c++17" option.

Some examples one can use to compile the code:

export GCC="gcc"; make
export GCC="g++"; make
export GCC="clang"; make

The GNU Scientific Library (GSL) must be installed to compile and run the program. If you are on Ubuntu or 
a debian based operating system, you may install GSL with:

	sudo apt-get install libgsl-dev

For other linux-based operating systems, you must use the coresponding package manager for your distribution.
***This code has only been tested on linux-based operating systems.***
######################################################################################################
# How to use:

To run with an input file: ./bin/TAPERSS -i [FILE]

#How to format the input file:
                                OPTION1 = [SELECTION1]
                                OPTION2 = [SELECTION2]
######################################################################################################
# OPTION1:

Example input shown in ./docs/ExampleInput.txt. All options are case-insensitive.

# NOTE:
Bash variables/Environment variables will NOT work inside the input file. If they must be used, use the 
command line options instead.

Below are the options for the input file.

Options:
    SEQUENCE          	Give the sequence to perform the specified action. 
			Case (uppercase or lowercase) must match that of the structure library file names

    OUTPUT-FILE       	Specify the output file name.

    WRITE-COORDINATES 	If equal to "TRUE", output file will be formatted as a PDB file with all coordinates written for each structure generated. 
		    	If "FALSE", only the index set, energy and structure filter RMSD values will be written for each structure generated.

    RUN-TYPE          Select action to be taken.
                      Available actions:
                      	"combinatorial"	- Perform combinatorial run
                        "from-index" 	- Create a single structure from sequence and given index
                        "idx-list" 	- Build several structures from an index list file, where the file contains several index sets (formatted as 1-1-1-1), 
                                     	  where each index set is separated by a newline. Index refers to the model in the fragment library.

    SINGLE-INDEX-LIST 				Indices in form "1-1-1-1-1" to be used with "from-index" action. 

    INDEX-LIST-FILE   				Specify the input file to be used with "idx-list" action.

    SECONDARY-STRUCTURE-FILTER (Optional)	Specifiy a structure filter in dot-bracket notation. Only accepts "(" and ".". 
						Example: to filter for a tetraloop with one closing base pair when studing a hexamer use notation: (....)

    OVERLAP-RMSD-LIMIT Set RMSD 		cutoff value for overlapped dinucleotides.

    WATSON-CRICK-OVERLAP-LIMIT 			Set RMSD cutoff value for the structure filter. 

    LIBRARY-PROTOTYPE 				Define a library file prototype. Give the location/filename of any library file but replace the dinucleotide sequence with XX. 
			       			The program will search for any files with that name as needed.
		                            	Example: ./LIBRARY_FILES/AA_library_combined.txt => ./LIBRARY_FILES/XX_library_combined.txt

    WATSON-CRICK-LIBARY-PROTOTYPE 		Same as the library option but for any Watson-Crick Libraries.

    STRUCTURE-COUNT-LIMIT-TYPE (Optional) 	Select the type of build limit (if any)
			       			Available Options:

                                                "blind"  - Once structure-count-limit structures are built, the calculation will stop. This option obeys the structure-filter,
                                                           so if it is set, it will stop once N structures which pass the filter are built.

                				***Work-in-Progress***: "energy" - All structures will be built as normal, but only the lowest energy structures will be saved to file,
							                     based on the structure-count-limit option.
									      					 Currently, if "WRITE-COORDINATES=TRUE" and STRUCTURE-COUNT-LIMIT <= 10, this option will work, but 
                                   may cause some memory issues.

    STRUCTURE-COUNT-LIMIT (Optional) 		Specifiy the total number of structures to by saved based on structure-count-limit-type option.
                                                If all structures are to be stored, then comment out line "#STRUCTURE-COUNT-LIMIT" and "#STRUCTURE-COUNT-LIMIT-TYPE".

    MEMORY-BUFFER-SIZE (Optional) 		Specifiy the size in Kilobytes (kB) of a memory buffer used to store output before writing to file. 
						This is to prevent slowdowns as a result of I/O processing. The larger the size, the larger the memory consumption 
						of the program, but the less I/O is performed. If not specified, a default size of 5MB is used.

######################################################################################################
# OPTION 2:

To run using flags:

./bin/TAPERSS -[flag] [option]

Below are the available flags. You may use both the inputfile and normal flags. Any conflicting options will be overwritten by whichever is specified last.

Flags:
    -i	[FILE]		Use to specify the input file (described above).
    -o  [FILE]		Specify output file to write to.
    -s  [SEQUENCE]	Sequence to use in performing action
    -c			Include to write PDB formatted output.
    --run-type [ACTION]				Specify the action to take. Same as RUN-TYPE described above
    --index [INDEX]    				Specify index to use with "idx" option selected.
    --index-list [FILE] 			Specify index file to use with "idx-list" option selected.
    --rmsd-lim [DECIMAL VALUE] 			Same as OVERLAP-RMSD-LIMIT.
    --wc-rmsd-lim [DECIMAL VALUE] 		Same as WATSON-CRICK-OVERLAP-LIMIT.
    --secondary-structure-filter [DOT-BRACKET]	Same as SECONDARY-STRUCTURE-FILTER.
    --lib-prototype [PROTOTYPE]			Same as LIBRARY-PROTOTYPE described above.
    --wc-lib-prototype [PROTOTYPE] 		Same as WATSON-CRICK-LIBRARY-PROTOTYPE described above.
    --str-count-limit-type [TYPE] 		Same as STRUCTURE-COUNT-LIMIT-TYPE.
    --str-count-limit [INTEGER] 		Same as STRUCTURE-COUNT-LIMIT.
######################################################################################################
# Potential problems and solutions:

The output files generated are very large: 
    As this is a combinatorial method, as the number of nucleotides in a sequence increases, the number of potential structures increases
    dramatically. It is not unexpected to generate over 1 billion structures for a 6-long nucleotide. To have more compact output, you may 
		try setting "WRITE-COORDINATES" to "FALSE" and use a low "OVERLAP-RMSD-LIMIT", less than 0.2 Angstroms. Using a lower "OVERLAP-RMSD-LIMIT" 
		will significantly reduce the number structures generated. A secondary structure filter will also reduce storage, but will only write those 
		which pass the filter to file.

The MEMORY-BUFFER-SIZE option may in some circumstances cause a "segmentation fault" error:
    If this is the case, try either removing or commenting out the MEMORY-BUFFER-OPTION to use the default size of 5MB, or attempt to increase/decrease the
		memory used manually.

The "energy" structure-build-limit-type option is causing segmentation faults, or strange outputs:
    This option is work in progress, and was implemented only to get a very small number of structures for a 6 nucleotide long sequence. 
		If you are to use this option, ensure the WRITE-COORDINATES option is set to "TRUE", and that the number of structures specified 
		in the "STRUCTURE-COUNT-LIMIT" is less than 10.


