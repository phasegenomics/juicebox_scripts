#this function returns the contig name expressed by a contig+bin string, typically of the format <contig>_<bin>
def getContigName(contigBinString):
    ret = ""
    
    #split the input contig+bin to set the target variables
    splitInputContigBin = contigBinString.split("_")
    #if we couldn't split anything, whatever we had before is the contig name
    if len(splitInputContigBin) == 1:
        ret = splitInputContigBin[0]
    #if we only got two things in the split, the first one is the contig name because the second is the bin
    elif len(splitInputContigBin) == 2:
        ret = splitInputContigBin[0]
    #otherwise, the contig is everything but the last thing in the split, concatenated back together, with _ between tokens
    else :
        for item in splitInputContigBin[0:len(splitInputContigBin)-1]:
            ret += item + "_"
        #trim off the final _
        ret = ret[0:len(ret)-1]
    
    return ret

#this funtion returns the bin name expressed by a contig+bin string, typically of the format <contig>_<bin>
def getBinName(contigBinString):
    ret = ""
    
    #split the input contig+bin to set the target variables
    splitInputContigBin = contigBinString.split("_")
    #if we couldn't split anything, there is no bin name
    if len(splitInputContigBin) == 1:
        ret = ""
    #if we only got more than one thing in the split, the last one is the bin
    else :
        ret = splitInputContigBin[len(splitInputContigBin)-1]
    
    return ret

#this funtion returns true if line is a contig+bin line
def isContigBinLine(line):
    ret = 0
    if len(line) > 1:
        ret = (line[0:1] == ">")
    return ret

#this function returns the contig+bin portion of a contig+bin line
def getContigBinFromLine(line):
    ret = ""
    if isContigBinLine(line):
        #trim the >, split the line, and grab the first token, in case there is a long description on the line
        ret = line[1:].split()[0].strip()
    return ret

#this function returns the name of the contig+bin at index in fastafile
def getContigBinAtIndex(index, fastafile):
    #format of the read file is
    #
    #>contigName_bin<n><continuing_description>
    #<sequence_line_1>
    #<sequence_line_2>
    #...
    #<sequence_line_n>
    #
    #read the file line by line until we find the contig+bin we're looking for. Note that the present
    #implementation depends on the bins in the fasta file being sorted.

    #create a counter to track the index
    readCounter = 0
    #return variable
    ret = "none"
    
#    print "Searching for index " + str(index) + " in " + fastafile

    #read the fast file line-by-line, counting each contig+bin line, then return the contig+bin on the index line
    with open(fastafile) as file:
        for line in file:
            line = line.strip()
            #check if the line is a contig+bin line
            if isContigBinLine(line):
                if readCounter == index:
                    contigBin = getContigBinFromLine(line)
                    ret = contigBin
                    break
                else:
                    readCounter += 1
            else:
                continue

    return ret

def getReverseComplement(sequence):
    ret = ""
    for b in sequence[::-1]:
        if b == "A":
            ret += "T"
        elif b == "T":
            ret += "A"
        elif b == "G":
            ret += "C"
        elif b == "C":
            ret += "G"
        elif b == "a":
            ret += "t"
        elif b == "t":
            ret += "a"
        elif b == "g":
            ret += "c"
        elif b == "c":
            ret += "g"
        else:
            ret += b
    return ret

#this function returns the index of a contig+bin pair in fastafile
def getIndexOfContigBin(contigBin, fastafile):
    #format of the read file is
    #
    #>contigName_bin<n><continuing_description>
    #<sequence_line_1>
    #<sequence_line_2>
    #...
    #<sequence_line_n>
    #
    #read the file line by line until we find the contig+bin we're looking for. Note that the present
    #implementation depends on the bins in the fasta file being sorted.
    
    #return variable
    ret = -1

    #create a counter to track the index
    readCounter = 0

    #go through the fasta file line by line. when we encounter a contig line (the lines that start with ">"), check if
    #it has the correct contig+bin name. if it does, set the return variable to the index and break.
    with open(fastafile) as file:
        for line in file:
            line = line.strip()
            #check if the line is a contig+bin line
            if isContigBinLine(line):
                lineContigBin = getContigBinFromLine(line)
                if lineContigBin == contigBin:
                    ret = readCounter
                else:
                    readCounter += 1
            else:
                continue

    return ret

#this function gets an no-newlines contig string from a fasta file
def getUnbrokenContigFromFasta (fastaFile, contig, reverse):
    inTarget = 0
    sequence = ""
    
    with open(fastaFile) as file:
        for line in file:
            line = line.strip()
            if len(line) == 0:
                continue
            elif isContigBinLine(line):
                lineContig = getContigBinFromLine(line)
                if lineContig == contig:
                    inTarget = 1
                else:
                    if inTarget:
                        return sequence
                    inTarget = 0
            elif inTarget:
                if not reverse:
                    sequence = sequence + line
                else:
                    sequence = getReverseComplement(line) + sequence

    return sequence

#this function returns the sequence line length of a fasta
def getSequenceLineLengthFromFasta (fastaFile):
    if 1:
        return 80
    lineSize = 0
    with open(fastaFile) as file:
        for line in file:
            line = line.strip()
            if len(line) == 0:
                continue
            elif isContigBinLine(line):
                continue
            else:
                lineSize = len(line)
                break
    return lineSize

#this function prints a contig from a fasta file to stdout
def printContigFromFasta (fastaFile, contig, reverse):
    sequence = getUnbrokenContigFromFasta(fastaFile, contig, reverse)
    lineSize = getSequenceLineLengthFromFasta (fastaFile)

    printed = 0
#    print len(sequence), printed, ">",
    while len(sequence) - printed > lineSize:
        print sequence[printed:printed+lineSize]
        printed += lineSize
#    print printed,
    if len(sequence) - printed > 0:
        print sequence[printed:]
        printed += len(sequence[printed:])
#    print printed

#this function prints a contig from a fasta file to stdout
def printContigsFromFasta (fastaFile, contigs, reverse):
    sequence = ""
    i = 0
    for contig in contigs:
        sequence += getUnbrokenContigFromFasta(fastaFile, contig, i == 0)
        i+=1
    lineSize = getSequenceLineLengthFromFasta (fastaFile)
    
    printed = 0
    while len(sequence) - printed > lineSize:
        print sequence[printed:printed+lineSize]
        printed += lineSize
    if len(sequence) - printed > 0:
        print sequence[printed:]

#this function returns an array of all the contig names in a fasta file
def getAllContigsFromFasta (fastaFile):
    contigs = []
    with open(fastaFile, 'r') as file:
        for line in file:
            line = line.strip()
            if isContigBinLine(line):
                contig = getContigBinFromLine(line)
                contigs.append(contig)
    return contigs

#this function returns an array of all the contig names in a grouping file
def getAllContigsFromGrouping (groupingFile):
    contigs = []
    with open(groupingFile) as file:
        for line in file:
            line = line.strip()
            if line[0:1] == "#":
                continue
            else:
                #line format: contig_ID(local)	contig_name	contig_rc	orientation_Q_score	gap_size_after_contig
                contig = line.split()[1]
                contigs.append(contig)
    return contigs

#this function returns a dictionary of the contigs in a fasta file to their broken sequences
def getContigDictionary (fastaFile):
    contigs = {}
    with open(fastaFile, 'r') as file:
        contig = None
        sequence = None
        for line in file:
            line = line.strip()
            if isContigBinLine(line):
                if contig is not None:
                    contigs[contig] = sequence
                contig = getContigBinFromLine(line)
                sequence = ""
            else:
                sequence += line + "\n"
        if contig is not None:
            contigs[contig] = sequence
    return contigs

#this function returns a dictionary of the contigs in a fasta file to their unbroken sequences
def getUnbrokenContigDictionary (fastaFile):
    contigs = {}
    with open(fastaFile, 'r') as file:
        contig = None
        sequence = None
        for line in file:
            line = line.strip()
            if isContigBinLine(line):
                if contig is not None:
                    contigs[contig] = sequence
                contig = getContigBinFromLine(line)
                sequence = ""
            else:
                sequence += line
        if contig is not None:
            contigs[contig] = sequence
    return contigs





