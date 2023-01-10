import sys, re, csv, os.path
import tkinter as tk       # sudo apt install python3-tk
from tkinter.filedialog import askdirectory
from pathlib import Path


### Functions
def findFirst(sequence, target, startIndex=0):
    position = sequence.find(target, startIndex)
    if(position>-1):
        return int(position)
    else:
        return "\'" + target + "\' not found in sequence."

def lookupAA(codon: str):
    aminoAcid = dna2amino.get(codon, 'Z') # Z not assoc. w/ AA
    return aminoAcid

def translate(dnaSequence, index=0):
    seqCodonList = [dnaSequence[i:i+3] for i in range((len(dnaSequence)-3), index-1, -3)]
    aaList = []
    for codon in seqCodonList:
        aaList.insert(0, lookupAA(''.join(codon)))
    return aaList

def reverse(input):
    return list(reversed(input))

def complement(input):
    # need error handling for N/uncalled bases
    inputComplement = []
    for base in input:
        inputComplement.append(baseComplement.get(base))
    return inputComplement

def transform(sequence):
    if isForward is False and isCoding is False: # most common case
        return complement(reverse(sequence))
    elif isForward is True and isCoding is False:
        return complement(sequence)
    elif isForward is False and isCoding is True:
        return reverse(sequence)
    else:
        return sequence

def inFrame(sequence):
    if 'not found' in sequence:
        return 'Error: End tag/stop codon not found in sequence.'
    elif len(sequence) % 3 == 0:
        return 'Sequence is in frame.'
    else:
        return 'Sequence is not in frame.'

def aa2dna(string):
    stringDNAlist = []
    key_list = list(dna2amino.keys())
    val_list = list(dna2amino.values())
    for i in string:
        val = val_list.index(i)
        stringDNAlist.append(key_list[val])
    return ''.join(stringDNAlist)

def createOutFile(path):
    fields = ['Sample_Name', 'Chromat_id', 'Read_id', 
              'Version', 'Length', 'Original Seq', 
              'Transformed Seq', 'ORF', 'ORF Base Count' 'In Frame', 
              'Amino Seq']
    global outfile
    outfile = (path + '/results_' + os.path.basename(path) + '.csv')
    with open(outfile, 'w') as csvfile:
        filewriter = csv.writer(csvfile)
        filewriter.writerow(fields)

def processFiles(path):
    source_files = Path(path).glob('*.seq')
    for file in source_files:
        with file.open('r') as f:
            data = f.readlines()
            # extract attributes from 1st line
            values = re.findall(r'(?<=\=)\w+', data[0])
            # process sequence from 2nd line, append to values
            #values.extend(analyze(testInput))
            #print((analyze(data[1])))
            values.extend(analyze(data[1])) # actual file content
            # write values to file
            with open(outfile, 'a') as f:
                csv.writer(f).writerow(values)

def analyze(sequence: str):
    ## Transform input into 5'-3' coding strand
    seqDNA = sequence.upper()
    seqDNAlist = list(sequence.upper())
    transformedDNA = transform(seqDNAlist)
    #seqAAlist = translate(transformedDNA)

    ## Find index of start codon
    startTagIndex = findFirst(''.join(transformedDNA), startTagDNA) # DNA
    #startTagIndex = findFirst(''.join(seqAAlist), startTag) # AA

    ## Create open reading frame, ORF
    stopTagIndex = findFirst(''.join(transformedDNA), stopTagDNA)
    if type(stopTagIndex) == int:
        orf = ''.join(transformedDNA)[startTagIndex:(stopTagIndex+len(stopTagDNA))]
    else:
        orf = stopTagIndex
    
    ##  Create encoded amino acid sequence
    if 'not found' in orf:
        aminos = 'Error: End tag/stop codon not found in sequence.'
    else:
        aminos = ''.join(translate(orf, startTagIndex))
    
    return [seqDNA, ''.join(transformedDNA), orf, len(orf), inFrame(orf), aminos]

### Variables
testInput = 'atgataggcatccatcggcttgctatgcttcgcgaccgtccctatgatcgcggcgccgaccaaaactgcgctcgtcaacgtatgcgtgccatgatgatgcaacatattaacttcagtacttggcccaactccatagcgagccattatgacattagtcatgtcccaagggctttccttgcaggtattgtatcttggacttcgcaggctaggagacatataagtcgcattcgtgtcgcggaaacaaaagatattgattgcttgatggctacgttcggctgtgacagggtcacaacttcaaccaaggcaattagcgccggactgcatttactcggatttcgggacaaaattcccgcagtgcgcagtcccttagtaccgccggataagcacggccgctcctttagctttcatgcttcggttgtactacctcaacataaactcccgtcccgcgtaggagtcgcaatcgtccgcatcgacccagcaccctgtgacgtggcaaacccgagtcgtttagaaaacgaagcaagtaatgagaagcagaactggattgcgtcatgcgcaggggcggacctttttcaagtgtcacaagcgacctgccttggtgcgtgtgcttccaataggttaatttcaggaagcactctgtattggattaggaggagttggatgaccagagacccgggtaatccattaaccctgagtgttgcaggattacgggttctggccacatatctagagaaaggggttttggagccaccaaaaatgagacgacatggtacgattgaggcaggctttttaaatgtaatgaaccccaaacaacagcgggtgaaacgatacacttcaagggctattttgcacgtgggatgttataaatcaggcctgcagggtttcgatagagtgtgtaccaacgccttccccccggccaagatcgtttccagggtattcgcagaacaaaaactcatctcagaagaggatctgtga'
dna2amino = {
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
baseComplement = {
    'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'
}
startTagAA = 'M' # DNA bp=ATG / AA=M
startTagDNA = 'ATG'
#startTag = input("Enter custom start tag:") # override default
stopTagAA = 'EQKLISEEDL'
stopTagDNA = 'GAACAAAAGCTTATTTCTGAAGAGGACTTG'
#stopTag = input("Enter custom end tag:") # override default
isForward = False  # from prompt
isCoding = False # from prompt


### Process
sourceLocation = askdirectory(title="Select Folder with Sequence Files", initialdir="/home/nick/projects/CodonAnalysis/samples")
createOutFile(sourceLocation)
processFiles(sourceLocation)

'''
To Do:
- test to ensure complement/reverse work
- should also sort files or output rows alpha
- should original name be preserved in spreadsheet
- finalize main function
- make UX improvements
- replace tk with builtin? something in os?
'''
