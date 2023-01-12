# CodonAnalysis
Tool to analyze DNA sequences to identify relevant amino acids. 

Initially in python with potential for a version in a compiled language. 

## Inputs 
* path to directory of files containing DNA sequences
* prompt for sequence direction (forward/coding vs reverse/complement)
    * May need to reverse the sequence direction then generate complementary base sequence
* nucleotide sequence of start codon (or arbitrary string)
    * default: ATG
* amino acid sequence of tag/ stop codon
    * default: GAACAAAAGCTTATTTCTGAAGAGGACTTG <sup>[3]</sup>

## Procedure
Given a DNA sequence (~1k bases), a start codon, and a stop codon...
1. Read DNA sequence base-wise in the specified direction 
2. Identify location of start codon
3. Identify location of tag
4. Create open reading frame (ORF) between start and tag
5. Translate into amino acids backwards from tag
4. Continue translation until start codon reached
5. Report following data points:
    * Sample metadata
    * Original and transformed sequence
    * ORF sequence and nucleotide count
    * In frame? yes/no
    * AA sequence from start to stop
  
## Resources
1. [DNA to Protein in Python 3](https://www.geeksforgeeks.org/dna-protein-python-3/)
2. [DNA and RNA codon tables](https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables)
3. [myc-tag](https://en.wikipedia.org/wiki/Myc-tag)
