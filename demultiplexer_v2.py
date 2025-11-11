import sys, os
import pandas as pd
import re
import distance
import argparse
import gzip


def countsin(inLoc):
    """
    Takes saved count file and reads it into a counts nested list.
    :param inLoc: counts file
    :return: nested list. counts nested list. [[read, total number, unique number],...]
    """
    countFile = open(inLoc, "r").readlines()
    counts=[]
    for i in range(1, len(countFile)):
        temp = countFile[i].rstrip().split(",")
        counts.append([temp[0][8:], temp[1], temp[2]])
    return counts

def CSVWriter (iterable, outLoc, header="", ):
    """
    Writes an iterable to a CSV file.
    :param iterable: List of list
    :param outLoc: file location. Where to place it.
    :param header: header of the CSV file
    :return: 1
    """
    if not iterable:
        print ("nothing to write")
        return 0

    out = open(outLoc, 'w')

    if header:
        out.write(header+'\n')

    #Only works if iterable is a nested list
    for member in iterable:
        for item in member:
            out.write(str(item)+',')
        out.write('\n')

    print("write to "+outLoc+" successful.")
    return 1

def reverseComplement(read):
    """
    :param str. a read in string format
    :return: str. the reverse complement of that read in string format
    """
    comp=[]

    for j in read:
        if j == 'A': comp.append('T')
        elif j == 'T': comp.append('A')
        elif j == 'C': comp.append('G')
        elif j == 'G': comp.append('C')
        else: comp.append('N')
    comp.reverse()
    comp=''.join(comp)
    return(comp)

def fastqParser(fileLoc, revComp = True):
    """
    Reads a fastq formatted file and returns a list of reads in string format.
    :param fileLoc: fastq file, reads are assumed to be represented as 4 lines, the second of which is the basecall.
    Can take .gz compressed files.
    :param revComp: default True, should the function return the reverse complement of the reads
    :return: list of reads in string format
    """

    if fileLoc[-3:] == ".gz": #Checks if file is compressed by .gz format
        fastq = gzip.open(fileLoc, 'rt').readlines()
    else:
        fastq = open(fileLoc, 'r').readlines()
    reads =[]

    for i in range(len(fastq)): #Assumes fastq files represent reads as 4 lines, the second of which are the basecalls
            if revComp:
                if (i - 1) % 4 == 0:
                    reads.append(reverseComplement(fastq[i].upper().rstrip()))
            else:
                if (i - 1) % 4 == 0:
                    reads.append(fastq[i].upper().rstrip())
    if not reads[-1]: reads = reads[:-1] #Sometimes adds empty string at the end of the list. This removes it.
    return reads

def uniqueFilter(readList, barcodeLength, barcodeMismatch=0):
    """
    Filters reads, collapses identical reads into a more compact object. Identifies duplicate reads by 3' barcode
    :param readList: list of reads
    :param barcodeLength: length of barcode at 3' end of read used to find unique reads. set to 0 if no barcode
    is used.
    :param barcodeMismatch: how similar can barcodes be and still call them duplicate reads. 0 means that they
    need to be exactly the same. 1 means if there is only one nt difference they will be called duplicates.
    Until the algorithm is sped up, anything other than 0 will be intensely slow. Barcodes are used to judge
    how many of those reads are unique
    :return: List. returns counts list of the structure [[read sequence, total reads, unique reads],...]
    """
    barcoded = []
    counts = []
    if not readList: return counts
    if barcodeLength:
        for read in readList:
            barcoded.append([read[:-barcodeLength],read[-barcodeLength:]])
        barcoded = sorted(barcoded, key = lambda x:x[0])
        tempBarcodes = []
        currRead = barcoded[0][0]
        for i in range(len(barcoded)):
            if i %100000 == 0: print ("{:.0%}".format((i/len(barcoded))))
            if barcoded[i][0] == currRead:
                tempBarcodes.append(barcoded[i][1])
            else:
                counts.append([currRead, len(tempBarcodes), uniqueNumber(tempBarcodes, mismatch=barcodeMismatch)])
                currRead = barcoded[i][0]
                tempBarcodes = [barcoded[i][1],]
        counts.append([currRead, len(tempBarcodes), uniqueNumber(tempBarcodes, mismatch=barcodeMismatch)])
        counts = sorted(counts, key = lambda x:x[2], reverse=True)
        return counts
    else:
        readList = sorted(readList)
        currRead = readList[0]
        n=0
        for i in range(len(readList)):
            if i%100000 == 0: print ("{:.0%}".format(i/len(readList)))
            if readList[i] ==  currRead:
                n+=1
            else:
                counts.append([currRead,n, n])
                currRead = readList[i]
                n = 0
        counts.append([currRead,n,n])
        counts = sorted(counts, key = lambda x:x[2], reverse = True)
        return counts

def uniqueNumber(list, mismatch=0):
    """
    Used to judge how many barcodes are unique. If mismatch is greater than 0, it uses hamming distance to
    judge how different items are. The algorithm is slower if the mismatch is greater than 0
    :param list: list of string items
    :param mismatch: how many differences do two items need to have before they are not called duplicates
    :return: int. How many of the items in that list are unique
    """
    if not list:
        return 0
    if not mismatch: #Much faster to do it this way if barcodes have to perfectly match
        return len(set(list))
    #This algorithm is slightly slower than above
    n=1
    list = sorted(list)
    for i in range(len(list)-1):
        if distance.hamming(list[i], list[i+1]) > mismatch: n+=1
    return n


def hitFinder2 (target, r1, r2, ham=0):
    hits = []
    for i in range(0, len(r2)):
        if len(target) > len(r2[i]): continue
        if distance.hamming(target, r2[i][:len(target)]) <= ham:
            hits.append(r1[i])
    return hits

def deMultiplexer(name, barcode, ranMerLen, r1, r2, umi_ham=0, bc_ham=0):
    """
    Takes the name of the sample, the barcode sequence, and the length of the randomMer. Finds the barcode in read1,
    returns the corresponding reads in read2, performs the read count operation and outputs it to a folder with the
    samplename.csv
    :param name: experiment name
    :param barcode:
    :param ranMerLen:
    :param r1:
    :param r2:
    :param mismatch:
    :return:
    """
    print(name, barcode, ranMerLen)
    hits = hitFinder2(barcode, r1, r2, bc_ham)

    counts = uniqueFilter(readList=hits, barcodeLength=ranMerLen, barcodeMismatch=umi_ham)
    return counts

def csv_from_excel(inLoc):
    xls = pd.ExcelFile(inLoc)
    df = xls.parse(sheetname="Sheet1", index_col=None, na_values=['NA'])
    df.to_csv('temp.csv', index=False)

def multipleLigationTrim(read, ranmer):
    patt = re.compile("AG[ATGC]{" + str(ranmer) + "}AGATCGGAAGAG")

    match = patt.search(read)

    if not match:
        return read
    else:
        return read[:match.start()]

########################################################################################################################

if __name__=="__main__":

    parser = argparse.ArgumentParser(description='Use inline barcodes to demultiplex reads')
    parser.add_argument('-r1', '--read1', type=str, help='Read1 location')
    parser.add_argument('-r2', '--read2', type=str, help='Read2 location')
    parser.add_argument('-m', '--manifest', type=str, help='Manifest location')
    parser.add_argument('-o', '--output', type=str, help='Output location')
    parser.add_argument('-s', '--sample', type=str, help='fastq sample names (e.g., S01)')
    parser.add_argument('-x', '--distance', type=int, help='inline barcode hamming distance', default=0)

    args = parser.parse_args()
    r1Loc = args.read1
    r2Loc = args.read2
    manifestLoc = args.manifest
    outFolder = args.output
    bc_ham = args.distance
    sample_name = args.sample

    r1 = fastqParser(r1Loc)
    r2 = fastqParser(r2Loc, revComp=False)
    assert len(r1) == len (r2)

    if not os.path.exists(outFolder):
        os.makedirs(outFolder)

    header = True  # Does the manifest have a header
    manifest = open(manifestLoc, 'r')
    if header: next(manifest) #skips header
    for line in manifest:
        line = line.split(",")
        primer_ID = line[0]
        barcode = line[2].rstrip() + line[3].rstrip()
        ranMerLen = int(line[4]) + 2
        counts = deMultiplexer(primer_ID, barcode, ranMerLen, r1, r2, umi_ham=0, bc_ham=bc_ham)
        counts2 = []
        for count in counts:
            temp_count = count[0][len(line[2].rstrip()):]
            temp_count = multipleLigationTrim(temp_count, ranMerLen - 2)
            temp_count = multipleLigationTrim(temp_count, ranMerLen - 2)
            temp_count = multipleLigationTrim(temp_count, ranMerLen - 2)
            counts2.append([temp_count,count[1],count[2]])

        outLoc1 = os.path.join(outFolder, f"{sample_name}_{primer_ID}_counts.csv")
        outLoc2 = os.path.join(outFolder, f"{sample_name}_demultiplex_dummy.txt")
        # This creates the real demultiplexed reads/counts files within each fastq sample
        CSVWriter(counts2, outLoc=outLoc1, header="Sequence, TotalReads, UniqueReads")
        # Create a dummy file per fastq sample (to avoid running the snakemake demultipler code multiple times)
        with open(outLoc2, "w") as f:
            f.write('This is to show the demultiplexing is done')

