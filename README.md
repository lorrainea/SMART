SMART: Supermaximal Approximate Repeats Tool
===

SMART is a program used to compute supermaximal k-mismtach repeats given a DNA or Protein sequence.

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
	 -i <file>  fasta/fastq input file
	 -o <file>  output file
	 -l <int>   minimum length of supermaximal repeats
	 -k <int>   maximum number of mismatches
	 -r <int>   0 for O(n) space for rmqs. 1 for O(nlogn) space for rmqs. Default: 1
	 -t <int>   0 not to trim to optimize e-value and 1 to trim to optimize e-value. Default: 1

	 -h         help
```
<b>Example</b>:

```
 ./smart -i ./data/ecoli.fasta -o ecoli.out -k 1 -l 1000
```
<b>Citation</b>
```
Lorraine A K Ayad, Panagiotis Charalampopoulos, Solon P Pissis
SMART: SuperMaximal Approximate Repeats Tool
Bioinformatics, 2019
https://doi.org/10.1093/bioinformatics/btz953
```
<b>License</b>: GNU GPLv3 License; Copyright (C) 2019 Lorraine A.K. Ayad, Panagiotis Charalampopoulos and Solon P. Pissis.
