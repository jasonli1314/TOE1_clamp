#Needs a CSV file with columns containing name, sequence, number of reads, number of unique reads
#Needs a FASTA file database with RNA genes containing exactly 50n tails (uses the 50n tail to identify mature end)
#Aligns reads with RNA genes from the FASTA file
#Returns a "Taildata" file listing sequences, counts, identified gene, 3' end positions, tail length and sequence.

from glob import glob
import csv
import argparse
import os

#Generates a dictionary of sequences with gene names as keys from a FASTA file (removes '>' in the name and '/n' and '/r' returns in the sequences)
def databasedict(fastafile):
    seqdict = {}
    seqs = fastafile.split('>')
    for seq in seqs[1:]:
        key = seq.splitlines()[0]
        seqdict[key] = ''.join(seq.splitlines(True)[1:]).replace('\n', '').replace('\r', '')
    return seqdict

def get_tail_info(RNAseqdict, datafile, seqcol, readcol):
        fcsv = open(datafile)
        Masterlist = [['Sequence', '#Unique Reads', 'Gene', '3-end', 'Tail length', 'Tail seq', 'Notes']]
        firstline = True

        for line in csv.reader(fcsv):
                if firstline:
                        firstline = False
                        continue

                # Generates the Masterlist for the Taildata file
                seq = line[seqcol]  # Grabs the sequence
                reads = line[readcol]  # Grabs the number of unique hammed reads

                # Identifies the best matching gene(s) for the current read
                linenumber = -1
                topgenescore = 0
                matchgeneseq = []
                matchgenelist = []
                matchlinenumber = []
                matchposition = []
                matchlength = [0]

                for name, geneseq in RNAseqdict.items():
                        linenumber += 1
                        x = 0
                        genescore = 0

                        while x <= (len(geneseq) - len(seq)):
                                y = 0  # counts the number of matches
                                m = 0  # number of mismatches
                                while y < (len(seq)) and m == 0:  # m==0, number of mismatches allowed
                                        if seq[y] == geneseq[x + y]:
                                                y += 1
                                        else:
                                                m += 1
                                if y > genescore:
                                        genescore = y
                                        geneposition = x
                                x += 1

                        if genescore == topgenescore and genescore > 0:
                                matchgeneseq.append(geneseq)
                                matchgenelist.append(name)
                                matchlinenumber.append(linenumber)
                                matchposition.append(geneposition)
                                matchlength.append(genescore)

                        if genescore > topgenescore:
                                matchgeneseq = [geneseq]
                                matchgenelist = [name]
                                matchlinenumber = [linenumber]
                                matchposition = [geneposition]
                                matchlength = [genescore]
                                topgenescore = genescore

                # Identifies the 3'end, tail length and tail sequence
                g = -1
                toss = 0
                errmsg = ''

                for geneseq in matchgeneseq:
                        g += 1
                        if matchlength[g] < len(seq) and toss == 0:
                                x = 1
                                match = 0
                                while (matchlength[g] + x) < len(seq) and (matchposition[g] + matchlength[g] + x) < len(geneseq):
                                        if seq[matchlength[g] + x] == geneseq[matchposition[g] + matchlength[g] + x]:
                                                match += 1
                                        x += 1

                                        if x == 3 and match == 2:  # Tosses if last 2n of seq match gene after mismatch
                                                toss = 1
                                                errmsg = 'ERR(Mismatch)'
                                        # Tosses the read if it matches the gene 75% at any point 3n after the mismatch
                                        if x > 3 and float(match)/(x - 1) >= 0.75:
                                                toss = 1
                                                errmsg = 'ERR(Mismatch)'

                                if toss == 0 and x > 1:  # Looks for deletion or insertion if the tail is longer than 1n
                                        y = 1
                                        match = 0
                                        while (matchlength[g] + y + 1) < len(seq) and (matchposition[g] + matchlength[g] + y) < len(geneseq):
                                                if seq[matchlength[g] + y + 1] == geneseq[matchposition[g] + matchlength[g] + y]:
                                                        match += 1
                                                y += 1
                                                # Tosses the read if it matches the gene 75% at any point 3n after the mismatch
                                                if y > 3 and float(match) / (y - 1) >= 0.75:
                                                        toss = 1
                                                        errmsg = 'ERR(Insert)'

                                        y = 1
                                        match = 0
                                        while (matchlength[g] + y) < len(seq) and (matchposition[g] + matchlength[g] + y + 1) < len(geneseq):
                                                if seq[matchlength[g] + y] == geneseq[matchposition[g] + matchlength[g] + y + 1]:
                                                        match += 1
                                                y += 1
                                                # Tosses the read if it matches the gene 75% at any point 3n after the mismatch
                                                if y > 3 and float(match) / (y - 1) >= 0.75:
                                                        toss = 1
                                                        errmsg = 'ERR(Deletion)'
                                # Tosses if tail is longer than sequenced part of gene
                                if toss == 0 and x > matchlength[0]:
                                        toss = 1
                                        errmsg = 'ERR(Tail longer than body)'

                        if matchlength[g] == len(seq):
                                x = 0

                if toss == 0 and matchlength[0] < 10:  # Requires 10 or more matches
                        toss = 1
                        errmsg = 'ERR(Short)'

                if toss == 0:
                        end = matchposition[0] + matchlength[0] + 50 - len(geneseq)
                        tailseq = seq[matchlength[0]:]
                        taillength = len(tailseq)

                else:
                        matchgenelist = []
                        end = errmsg
                        taillength = ''
                        tailseq = ''

                notes = ''
                Masterlist.append([seq, reads, matchgenelist, end, taillength, tailseq, notes])
        fcsv.close()
        return Masterlist

########################################################################################################################
if __name__=="__main__":

        parser = argparse.ArgumentParser(description='Analyze tail lengths')
        parser.add_argument('-i', '--input', type=str, help='input directory')
        parser.add_argument('-d', '--database', type=str, help='snRNA fasta sequence+50bp tail')
        parser.add_argument('-n', '--sample', type=str, help='name of the sample fastq file (e.g., S01)')
        parser.add_argument('-o', '--output', type=str, help='output directory')
        parser.add_argument('-s', '--seqcol', type=int, help='Sequence column (0-index)', default=0)
        parser.add_argument('-r', '--readcol', type=int, help='Reads column (0-index)', default=2)
        args = parser.parse_args()
        input = args.input
        fastafile = args.database
        sample = args.sample
        output = args.output
        seqcol = args.seqcol
        readcol = args.readcol

        if not os.path.exists(output):
                os.makedirs(output)

        #Generates a dictionary (RNAseqdict) containing gene names and sequences from the RNA FASTA database
        f = open(fastafile)
        RNAseqs = f.read() #loads the snRNA FASTA file into snRNAseqs
        f.close()
        RNAseqdict = databasedict(RNAseqs)

        # all files associated with a fastq sample (e.g., S01)
        datafiles = glob(os.path.join(input, sample+'*_counts.csv'))
        for datafile in datafiles:
                Masterlist = get_tail_info(RNAseqdict=RNAseqdict, datafile=datafile, seqcol=seqcol, readcol=readcol)
                output_file = os.path.join(output, os.path.basename(datafile))
                output_file = output_file.replace('_counts.csv', '_tail.csv')
                with open(output_file, 'w', newline='') as f:
                        writer = csv.writer(f)
                        writer.writerows(Masterlist)
        dummy_file = os.path.join(output, f"{sample}_analyzer_dummy.txt")
        with open(dummy_file, "w") as d:
                d.write('This is to show the tail analysis is done')