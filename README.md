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
* Mapping of codons to amino acids - [dictionary](https://www.geeksforgeeks.org/dna-protein-python-3/)
    * Science Fact<sup>TM</sup>: one or more codons map to a each amino acid
    * ```python
        dna2amino{
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        }
        ```

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
