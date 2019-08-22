SMART: Supermaximal Approximate Repeats Tool
===

SMART is a program which can be used in conjunction with any multiple sequence alignment program, to address this problem effectively and efficiently.

<b>Installation</b>: To compile SMART, please follow the instructions given in file INSTALL.

<b>INPUT</b>: A sequence in FASTA format. The input file is specified using the <b>-i</b> option.

<b>OUTPUT</b>: A BED file containing statistics of supermaximal repeats computed. The output file is specified using the <b>-o</b> option.

Columns 1 and 3: Length of supermaximal repeat. 

Column 2: Starting position of the first supermaximal repeat.

Column 4: Starting position of the second supermaximal repeat.

Column 5: Number of mismatches between the pair of supermaximal repeats.

Column 6: Identity score of the pair of supermaximal repeats.

Column 7: Raw score of the pair of supermaximal repeats.

Column 8: E-value of the pair of supermaximal repeats.

```
 Usage: smart <options>
 Available options:
	 -i <str>   Fasta input file
	 -o <str>   Output file
	 -l <int>   Minimum length of supermaximal repeats
	 -k <int>   Exact number of mismatches

	 -h         help
```
<b>Example</b>:

```
 ./smart -i ./data/ecoli.fasta -o ecoli.out -k 1 -l 1000
```
<b>License</b>: GNU GPLv3 License; Copyright (C) 2019 Lorraine A.K. Ayad, Panagiotis Charalampopoulos and Solon P. Pissis.
