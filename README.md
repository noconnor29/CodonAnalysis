# CodonAnalysis
Tool to analyze DNA sequences to identify relevant amino acids. 

Initially in python with potential for a version in a compiled language. 

## Inputs 
* path to CSV file with DNA sequence(s)
* amino acid sequence of end tag sequence (required)
    * default: EQKLISEEDL <sup>[3]</sup>
* base sequence of start tag sequence (optional)

## Data Structures
* Variables for end and start tag
* Running list of 30 last read bases (to identify end tag)
* Running list of amino acids
* Mapping of codons to amino acids - [dictionary](https://www.geeksforgeeks.org/dna-protein-python-3/)
    * Science Fact<sup>TM</sup>: one or more codons map to a each amino acid
    * ```python
        codon2amino{
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
Given a DNA sequence (~1k bases), an end tag, and an optional start tag...
1. Read DNA sequence base-wise downstream (5' to 3') until end tag sequence is identified
2. Move reading frame upsteam from end tag in 3-base (codon) increments
3. Within the frame, bases are read downstream to identify codon and/or amino acid (AA). Append codon/AA to beginning of list of codons/AA.
    * areas of low sequence fidelity ("uncalled bases") will require error handling
4. Repeat Steps 2 & 3 until...
    - start tag (if supplied) is reached (i.e. indicies [0-2] = start tag), or
    - beginning of sequence is reached
  
## Resources
1. [DNA to Protein in Python 3](https://www.geeksforgeeks.org/dna-protein-python-3/)
2. [DNA and RNA codon tables](https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables)
3. [myc-tag](https://en.wikipedia.org/wiki/Myc-tag)
