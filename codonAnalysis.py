import csv
import tkinter as tk       # sudo apt install python3-tk
from tkinter.filedialog import askopenfilename

### Import base seq from file
tk.Tk().withdraw()      # part of the import if not using other tkinter f(x)

#filename = askopenfilename()
#inputfile ="\path\to\DNA_sequence_placeholder.csv" # use for testing
#file = open(filename, "r")
#seq = file.read()

### Variables
startTag = 'ATG'
#startTag = input("Enter custom start tag:") # override default
endTag = 'EQKLISEEDL'
#endTag = input("Enter custom end tag:") # override default

rawInput = 'atgataggcatccatcggcttgctatgcttcgcgaccgtccctatgatcgcggcgccgaccaaaactgcgctcgtcaacgtatgcgtgccatgatgatgcaacatattaacttcagtacttggcccaactccatagcgagccattatgacattagtcatgtcccaagggctttccttgcaggtattgtatcttggacttcgcaggctaggagacatataagtcgcattcgtgtcgcggaaacaaaagatattgattgcttgatggctacgttcggctgtgacagggtcacaacttcaaccaaggcaattagcgccggactgcatttactcggatttcgggacaaaattcccgcagtgcgcagtcccttagtaccgccggataagcacggccgctcctttagctttcatgcttcggttgtactacctcaacataaactcccgtcccgcgtaggagtcgcaatcgtccgcatcgacccagcaccctgtgacgtggcaaacccgagtcgtttagaaaacgaagcaagtaatgagaagcagaactggattgcgtcatgcgcaggggcggacctttttcaagtgtcacaagcgacctgccttggtgcgtgtgcttccaataggttaatttcaggaagcactctgtattggattaggaggagttggatgaccagagacccgggtaatccattaaccctgagtgttgcaggattacgggttctggccacatatctagagaaaggggttttggagccaccaaaaatgagacgacatggtacgattgaggcaggctttttaaatgtaatgaaccccaaacaacagcgggtgaaacgatacacttcaagggctattttgcacgtgggatgttataaatcaggcctgcagggtttcgatagagtgtgtaccaacgccttccccccggccaagatcgtttccagggtattcgcagaacaaaaactcatctcagaagaggatctgtga'
seq = rawInput.upper()
### Read seq content to upper list
seqList = list(seq)


def firstFind(sequence, target, startIndex=0):
    position = sequence.find(target, startIndex)
    if(position>-1):
        return position
    else:
        return "\'" + target + "\' not found in sequence."

## Find index of start codon
startIndex = firstFind(seq, startTag)
print(startIndex)
#firstFind(seq, startTag)

## Starting with startTag, convert DNA to AA

## Read from start to stop codon
print(firstFind(seq, endTag, startIndex))