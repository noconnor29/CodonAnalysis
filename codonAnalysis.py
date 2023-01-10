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

def lookupAA(codon):
    aminoAcid = dna2amino.get(codon, 'No AA match')
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
    if len(sequence) % 3 == 0:
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
              'Transformed Seq', 'ORF', 'In Frame', 
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
            values.extend(analyze(testInput))
            #values.extend(analyze(data[1])) # actual file content
            # write values to file
            with open(outfile, 'a') as f:
                csv.writer(f).writerow(values)

def analyze(sequence: str):
    ## Transform input into 5'-3' coding strand
    seqDNA = sequence.upper()
    seqDNAlist = list(sequence.upper())
    transformedDNA = transform(seqDNAlist)
    seqAAlist = translate(transformedDNA)

    ## Find index of start codon
    startTagIndex = findFirst(''.join(seqDNAlist), startTagDNA) # DNA
    #startTagIndex = findFirst(''.join(seqAAlist), startTag) # AA
    #print('Start Tag Index: ' + str(startTagIndex))

    ## Create open reading frame, ORF
    stopTagIndex = findFirst(seqDNA, stopTagDNA) # encounter error here if stop tag not found
    # need desired workflow. If stoptag not found, read until end or error out?
    #print('Stop Tag Index: ' + str(stopTagIndex))
    orf = seqDNA[startTagIndex:(stopTagIndex+len(stopTagDNA))]

    ##  Create encoded amino acid sequence
    aminos = ''.join(translate(orf, startTagIndex))

    ## Return results
    return [startTagIndex, stopTagIndex, orf, aminos]


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
    'A':'T', 'C':'G', 'G':'C', 'T':'A'
}
startTag = 'M' # DNA bp=ATG / AA=M
startTagDNA = 'ATG'
#startTag = input("Enter custom start tag:") # override default
stopTag = 'EQKLISEEDL'
stopTagDNA = 'GAACAAAAACTCATCTCAGAAGAGGATCTG'
#stopTag = input("Enter custom end tag:") # override default
isForward = True  # from prompt
isCoding = True # from prompt


### Process
sourceLocation = askdirectory(title="Select Folder with Sequence Files")
createOutFile(sourceLocation)
processFiles(sourceLocation)

## Report Output
'''
print('\nFinal Results')
print('Original Sequence: ' + seqDNA)
print()
print('DNA ORF: ' + orf)
print('Stop Tag in frame?: ' + inFrame(orf))
print()
print('AA Seq: ' + aminos)
'''

'''
To Do:
- Develop error handling in stopTagIndex if not found
    - error due to "uncalled base" i.e. unknown/low fidelity
- finalize main function
- make UX improvements
- replace tk with builtin? something in os?
'''
