# CodonAnalysis
Tool to analyze DNA sequences to identify relevant amino acids. 

Initially in python with potential for a version in a compiled language. 

## Inputs 
* path to directory of files containing DNA sequences
* prompt for sequence direction (forward/coding vs reverse/complement)
    * in general, given reverse, non-coding sequence. Need to reverse the recieved direction (from 5'-3'to 3'-5') then generate complementary base sequence
* base sequence of start tag/codon
    * default: ATG
* amino acid sequence of stop tag/codon
    * default: EQKLISEEDL <sup>[3]</sup>

## Data Structures
* Variables for end and start tag
* Running list of amino acids
* Mapping of codons to amino acids
    * Science Fact<sup>TM</sup>: one or more codons map to a each amino acid

## Procedure
Given a DNA sequence (~1k bases), a start tag/codon, and an stop tag/codon...
1. Read DNA sequence base-wise in the specified direction 
2. Identify location of start tag
3. Translate into amino acids backwards from end tag
4. Continue translation until stop tag/stop codon reached
5. Report following data points:
    * Original sequence
    * DNA sequence from start to stop
    * AA sequence from start to stop
    * In-frame? yes/no
  
## Resources
1. [DNA to Protein in Python 3](https://www.geeksforgeeks.org/dna-protein-python-3/)
2. [DNA and RNA codon tables](https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables)
3. [myc-tag](https://en.wikipedia.org/wiki/Myc-tag)
