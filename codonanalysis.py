import sys, re, csv, os.path
from pathlib import Path
# pip install python-dotenv
from dotenv import load_dotenv

"""Functions"""
def main():
    global sourceLocation
    sourceLocation = (input('Enter path to .seq files: ').replace("'",""))
    print(sourceLocation)
    menu()
    global isCoding
    isCoding = getOrientation()
    global startSeqDNA
    startSeqDNA = getStart()
    global stopSeqDNA
    stopSeqDNA = getTag()
    createOutFile(sourceLocation)
    processFiles(sourceLocation)

def menu():
    print('Strand Orientation:')
    print('    1 -- Coding/Sense')
    print('    2 -- Non-Coding/Antisense\n')
    print('All strands in a batch are processed with same orientation.\n')

def getOrientation():
    while True:
        selection = input('Please select the number representing the orientation: ')
        if selection not in ('1', '2'):
            for i in range(5): print('\n')
            print('\n' + '>>> ' + str(selection) + ' is not a valid option. ' \
                'Please enter either 1 or 2. <<<\n')
            menu()
        elif selection == '1':
            return True
        elif selection == '2':
            return False

def getStart():
    tag = input("Enter Nucleotides of Start Sequence. Press [Enter] for 'ATG': ")
    if tag == '':
        return 'ATG'
    else:
        return str(tag)

def getTag():
    tag = input('Enter Nucleotides of Tag. Press [Enter] for Myc-tag: ')
    if tag == '':
        return 'GAACAAAAGCTTATTTCTGAAGAGGACTTG'
    else:
        return str(tag)

def findFirst(sequence, target, startIndex=0):
    position = sequence.find(target, startIndex)
    if(position>-1):
        return int(position)
    else:
        return "\'" + target + "\' not found in sequence."

def lookupAA(codon: str):
    aminoAcid = dna2amino.get(codon, 'Z')
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
    if not isCoding:
        return complement(reverse(sequence))
    else:
        return sequence

def inFrame(sequence):
    if 'not found' in sequence:
        return 'Error: Tag not found in sequence.'
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
    fields = ['Source File', 'Sample_Name', 'Chromat_id', 'Read_id',
              'Version', 'Length', 'Original Seq',
              'Transformed Seq', 'ORF', 'ORF Base Count', 'In Frame',
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
            values = []
            values.append(os.path.basename(file))
            data = f.readlines()
            """extract attributes from 1st line"""
            values.extend(re.findall(r'(?<=\=)\w+', data[0]))
            """process sequence from 2nd line, append to values"""
            values.extend(analyze(data[1]))
            """write values to file"""
            with open(outfile, 'a') as f:
                csv.writer(f).writerow(values)
    print('\nResults: ' + sourceLocation + os.sep + outfile)

def analyze(sequence: str):
    seqDNA = sequence.upper()
    seqDNAlist = list(sequence.upper())
    transformedDNA = transform(seqDNAlist)
    """Find index of start codon"""
    startIndex = findFirst(''.join(transformedDNA), startSeqDNA)
    """Create open reading frame, ORF"""
    stopIndex = findFirst(''.join(transformedDNA), stopSeqDNA)
    if type(stopIndex) == int:
        orf = ''.join(transformedDNA)[startIndex:(stopIndex+len(stopSeqDNA)+3)]
        lenORF = len(orf)
    else:
        orf = stopIndex
        lenORF = 'N/A'
    """Create encoded amino acid sequence"""
    if 'not found' in orf:
        aminos = 'Error: Tag codon not found in sequence.'
    else:
        aminos = ''.join(translate(orf, startIndex))
    return [seqDNA, ''.join(transformedDNA), orf, lenORF, inFrame(orf), aminos]

"""Variables"""
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
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}
baseComplement = {
    'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'
}
''' Unused
startTagAA = 'M'
stopTagAA = 'EQKLISEEDL'
'''

"""Process"""
if __name__ == "__main__":
    main()
