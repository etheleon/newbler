#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import errno
import re
import shutil
import subprocess


from collections   import defaultdict
import subprocess

from Bio           import SeqIO
from Bio.Seq       import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet  import SingleLetterAlphabet
import pandas as pd

class Alignment:
    '''
    Methods executed in this order generates a MDR centric assembly (round2).
        1. doPile
        2. getReadsFromPileup
        3. roundTwoAssembly
    Meant for further analysis after PADI

    Input:
        * contig sequence
        * MDR sequence with gaps
        * MDR loc
        * Reads in MDR
    Outputs:
        * read pileup (alignment.doPile())
        * reads in MDR (alignment.getReadsFromPileup())
    '''

    def __init__ (self, rootPath, ko):
        self.rootPath = rootPath
        self.ko = ko
        self.start       = None
        self.end         = None
        self.contigList  = {}
        self.readInfo    = {}
        self.outputRecords = []
        print("Processing %s:" % ko)

    def doPile(self):
        """
        Generates a full pileup for each of the contigs.
        To be used later for assembly
        """
        pileup = "%s/out/pileup/%s" % (self.rootPath, self.ko)
        try:
            os.makedirs(pileup)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise  # raises the error again
        print("Initializing pileup for %s" % self.ko)
        self.__getMSALOC()
        self.__readContigs()
        self.__readMSA()
        #just want to test out how contig000001 looks like
        self.__readStatusPair()
        self.__parseFastQ()

    def getReadsFromPileUP(self):
        self.__getMSALOC()
        self.__storeTAXAinfo()
        self.__cutMSA()

    def __cutMSA(self):
        #msa alignment can be Rev and ntRev
        try:
            outputDir = "%s/out/preNewbler/%s" % (self.rootPath, self.ko)
            os.makedirs(outputDir)
        except OSError as e:
            if e.errno != errno.EEXIST: #keep quiet if folder already exists
                raise  # raises the error again
        print("extracting reads")
        contigListFile = "%s/out/pAss03/%s.msa" % (self.rootPath, self.ko)
        print(contigListFile)
        seqIter = SeqIO.parse(contigListFile, 'fasta')
        for contig in seqIter:
            mdr = self.__getSeq(contig, self.start,self.end)
            contigInMDR = len(mdr) > 0
            if (contigInMDR):
                pileupFH = SeqIO.parse("%s/out/pileup/%s/%s-%s" % (self.rootPath, self.ko, self.ko, contig.id), 'fasta')
                fullContig = next(pileupFH)
                try:
                    indexVal = str(fullContig.seq).index(mdr)
                    #print("%s: %s" % (contig.id, indexVal))
                    self.__extractReads(indexVal, indexVal + len(mdr), pileupFH)
                except ValueError:
                    #print("%s : Cannot find MDR in 5' strand. Trying revcom of MDR" % contig.id)
                    try:
                        indexVal = str(fullContig.seq).index(str(Seq(mdr).reverse_complement()))
                        #print("%s: %s" % (contig.id, indexVal))
                        self.__extractReads(indexVal, indexVal + len(mdr), pileupFH)
                    except ValueError as err:
                        print("%s-%s has issues:"%(self.ko, contig.id))
                        print(err)
                SeqIO.write(self.outputRecords, "%s/out/preNewbler/%s/%s" % (self.rootPath, self.ko, self.ko), "fasta")
            else:
                print("%s is empty" % contig.id)

    def __getSeq(self, seqRecord, start, end):
        nt = str(seqRecord.seq[start : end]).upper().replace("-", "")
        return nt

    def __extractReads(self, indexVal, howLong, iterator):
        for record in iterator:
            read = self.__getSeq(record, indexVal, howLong)
            if len(read) > 0:
                readID, contig  =  re.match("^(\d+)-(contig\d+)$", record.id).groups()
                r               =  self.readInfo[readID]
                #>58526338-contig00001-33057/1;  KO:K00927       start: 575      offset: 287
                header = "%s-%s/%s\tKO:%s\tstart:%s\toffset:%s" % (record.id, r['taxa'], r['readnum'], self.ko, indexVal, howLong)
                newseq = Seq(str(record.seq).upper().translate({ord(i):None for i in '-'}), SingleLetterAlphabet())
                newrecord = SeqRecord(newseq, id=header, name="", description="")
                self.outputRecords.append(newrecord)
                self.readInfo[readID]['readnum'] += 1
        return

    def __storeTAXAinfo(self):
        print("processing input file...")
        fq1 = self.rootPath + "/out/newbler/" + self.ko + "/input/" + self.ko + ".1.fq"
        fq2 = self.rootPath + "/out/newbler/" + self.ko + "/input/" + self.ko + ".2.fq"
        for file in [fq1, fq2]:
            for record in SeqIO.parse(file, "fastq"):
                readID, taxa = re.search("^(\d+)\|(\d+)", record.description).groups()
                self.readInfo[readID] = {'taxa' : taxa , 'readnum' : 1}
        print("done")

    def __getMSALOC(self):
        """
        grab the MDR location from MDR file originating from pass
        """
        if (self.start == None and self.end == None):
            file        =  self.rootPath + '/out/pAss11/' + self.ko + ".fna"
            record      =  next(SeqIO.parse(file, "fasta"))
            theMatch    =  re.search("msaStart:(\d+) msaEND:(\d+)", record.description)
            self.start  =  int(theMatch.group(1))
            self.end    =  int(theMatch.group(2))
            print("MDR start:%s end:%s" % (self.start, self.end))

    def __readContigs(self):
        """
        Stores full length contigs
        """
        path = self.rootPath+'/out/newbler/'+self.ko+"/454AllContigs.fna"
        for record in SeqIO.parse(path, 'fasta'):
            self.contigList[record.id] = {'fullseq': record.seq.upper()}

    def __readMSA(self):
        """
        Parses the MSA for the MDR region.
        Outputs the portion of the contig sequence recorded from the 454 output which matches the sequences from the msa in the MDR
        """
        path = self.rootPath + '/out/pAss03/' + self.ko + ".msa"
        for record in SeqIO.parse(path, "fasta"):
            contig = self.contigList[record.id]['fullseq'] #might not have this contig in the newbler cause .... shet something's seriously not right
            recseq = str(record.seq)[self.start:self.end]

            shrunk = recseq.replace('-', '').upper()
            substrResult = contig.find(shrunk)
            #its a reverse strand
            if substrResult == -1:
                revshrunk = str(record.seq.reverse_complement())[self.start : self.end].replace('-', '').upper()
                substrResult = contig.find(revshrunk)
                self.contigList[record.id]['mdr'] = { 'start': substrResult,
                                                            'end' : len(revshrunk) + substrResult,
                                                            'seq' : revshrunk,
                                                            'direction' : 'ntRev'
                }
            #not a reverse strand
            else:
                self.contigList[record.id]['mdr'] = {
                    'start': substrResult,
                    'end' : len(shrunk) + substrResult,
                    'seq' : shrunk,
                    'direction' : 'ntRev'
                }


    def __readStatusPair(self):
        """
        Newbler interprets was given the command to intepret the reads as paired end reads.
        We got newbler to process them as a pair instead as single reads.
        NOTE: we should compare the results from pair and single read mapping

        Currently the function only takes into account:
            1. reads which assembled into the same contig,
            2.  assigns the reads to the starting position of the pair, not the individual read position,
                this should be compared with using the read level information

        Definition: Distance
        Distance:
        for reads that map to the same contig: the distance between the halves
        for reads that Link contigs into scaffolds: the sum of the distances from the position of each half to the end of the contig. So, the total distance between the halves for these pairs would be the distance mentioned in the

        Template        Status  Distance        Left Contig     Left Pos        Left Dir        Right Contig    Right Pos       Right Dir       Left Distance   RightDistance
        332|201096|478844-478945|s_1    SameContig      115     contig02355     137     -       contig02355     22      +
        2334|561|4092117-4092218|s_1    SameContig      151     contig01613     106     +       contig01613     257     -
        5924|237|2291358-2291459|s_1    SameContig      144     contig00035     2174    +       contig00035     2318    -
        7317|914|1906065-1906166|s_1    SameContig      143     contig00018     2418    +       contig00018     2561    -
        7989|914|1904223-1904324|s_1    SameContig      150     contig00018     575     +       contig00018     725     -
        8070|48736|485184-485285|s_1    SameContig      139     contig00379     413     +       contig00379     552     -
        9763|1763|762968-763069|s_1     SameContig      150     contig01020     332     -       contig01020     182     +
        10918|165696|411626-411727|s_1  SameContig      144     contig00150     265     +       contig00150     409     -
        12391|1091|1788625-1788726|s_1  SameContig      156     contig00405     553     +       contig00405     709     -

        """
        # 454Pairstatus
        df = pd.read_csv(self.rootPath + "/out/newbler/"+self.ko+"/454PairStatus.txt", sep="\t")

#isAss =  == ['SameContig', 'Link', 'OneUnmapped']
        def splitNstore(row):
            """
            newbler seems to be miss labelling based on status the reads: alright assemblies -> FalsePair
            False Pairs are Paired End reads whose ends either:
                * map to different reference sequences or alignin contigs that are in different scaffolds
                * map to locations on a single reference sequence or assembled contig with a distance outside the expected library span
                * have an unexpected relative orientation
            eg.
            Template        Status  Distance        Left Contig     Left Pos        Left Dir        Right Contig    Right Pos    Right Dir       Left Distance   RightDistance
            51910|221065|1544081-1544182|s_1        FalsePair       144     contig00235     156     +       contig00235 300      -
            """
            acceptedStatus = set(['SameContig', 'FalsePair'])
            if row['Status'] in acceptedStatus: #we should consider other statuses; yep come back to bite me in the butt
                splitNstore_sameContig(row)
            #switcher = {
                #'SameContig': splitNstore_sameContig(row)
                #'Link':
                #'OneUnmapped':
            #}
            #switcher.get(row['Status'], "not inside")

        def splitNstore_sameContig(row):
            readID = row['Template'].split("|")[0]
            isPos = row['Left Dir'] == '+'
            try:
                self.readInfo[readID] = {
                    'parent': row['Left Contig'],
                    'readone':  int(row['Left Pos']),
                    'readtwo':  row['Right Pos'],
                    'direction':  'forward' if isPos else 'reverse'
                }
            except:
                print("%s the parent:%s ; read1: %s; direction %s") % (readID, row['Left Contig'], row['Left Pos'], row['Left Dir'])
        df.apply(splitNstore, axis=1)

    def __guidedAlignment(self):
        cmd = "bwa mem"
        # incomplete: what i wanted to do is to have the mapping process passed onto bwa + 454ReadStatus, ie. get the read
        # which belong to

    def __parseFastQ(self):
        """
        Parses fastqfiles, stores then outputs the reads as fq pileups on the respective contigs.
        not really working. will begin work on new private method
        """
        def fixAlignment(rootPath, ko, contigID, debug=False):
            """
            temp fix for parseFastQ
            readStatus mapping only gives location for contig not read,
            ie.
               123456789
             --ATCGGGCAT  <contig> mapping position 1-4 (3 nts)
             CGATCG-----  <read>   maping  position 3-6 (3 nts)
             123456
            """
            #Part1: run muscle
            #"tests/out/pileup/K00927/K00927-contig000001"
            pileup = "%s/out/pileup/%s" % (rootPath, ko)
            inputFile = '%s/%s-%s' % (pileup, ko, contigID)
            cmd = "muscle -in %s -out %s.msa" % (inputFile, inputFile)
            print("Running muscle:\n%s"%cmd)
            if not debug:
                try:
                    subprocess.run(cmd, shell=True, check=True)
                except subprocess.CalledProcessError as err:
                    print("Error:\n", err.output)
            #Part2: Store MSA sequences
            msaed = {}
            seqIter = SeqIO.parse("%s.msa"%inputFile, 'fasta')
            for aligned in seqIter:
                #print(aligned.description)
                try:
                    readID = re.search("^(\S+)-\S+$", aligned.description).group(1)
                    msaed[readID] = str(aligned.seq)
                except AttributeError as err:
                    msaed[aligned.description] = str(aligned.seq)
            return msaed

        testing = False
        fq1 = self.rootPath + "/out/newbler/" + self.ko + "/input/" + self.ko + ".1.fq"
        fq2 = self.rootPath + "/out/newbler/" + self.ko + "/input/" + self.ko + ".2.fq"
        poshash = {} #defaultdict(list)
        for record in SeqIO.parse(fq1, "fastq"):
            readID = record.description.split("|")[0]
            if readID in self.readInfo:
                theParent = self.readInfo[readID]['parent']
                if theParent in self.contigList.keys(): #sometimes shorter contigs are not outputed so readInfo may not match the contig info. ie. read was assigned to a contig which was not printed
                    #print "yes Inside"
                    startPos = self.readInfo[readID]['readone']
                    direc = self.readInfo[readID]['direction']
                    back = "-" * int(len(self.contigList[theParent]['fullseq']) - startPos - len(record.seq))
                    if direc == 'reverse':
                        startPos = startPos -101
                        front = "-" * int(startPos-1)
                        seq = str(record.seq.reverse_complement())
                    else:
                        front = "-" * int(startPos-1)
                        seq = str(record.seq)
                    readSEQ = front+seq+back
                    readID = readID + "/1"
                    if theParent in poshash:
                        if testing:
                            if theParent == 'contig00001':
                                poshash[theParent][startPos].append((readID, readSEQ))
                        else:
                            poshash[theParent][startPos].append((readID, readSEQ))
                    else:
                        poshash[theParent] = defaultdict(list)
        for record in SeqIO.parse(fq2, "fastq"):
            readID = record.description.split("|")[0]
            if readID in self.readInfo:
                theParent = self.readInfo[readID]['parent']
                if theParent in self.contigList.keys(): #sometimes shorter contigs are not outputed so readInfo may not match the contig info. ie. read was assigned to a contig which was not printed
                    startPos = self.readInfo[readID]['readtwo']
                    direc = self.readInfo[readID]['direction']
                    if direc == 'reverse':
                        seq = str(record.seq)
                        front = "-" * int(startPos-1)
                    else:
                        seq = str(record.seq.reverse_complement())
                        startPos = startPos - 101
                        front = "-" * (int(startPos-1))
                    back = "-" * int(len(self.contigList[theParent]['fullseq']) - startPos - len(record.seq))
                    readSEQ = front+seq+back
                    readID = readID + "/2"
                    if theParent in poshash:
                        if testing:
                            if theParent == 'contig00001':
                                poshash[theParent][startPos].append((readID, readSEQ))
                        else:
                            poshash[theParent][startPos].append((readID, readSEQ))
                    else:
                        poshash[theParent] = defaultdict(list)
        #open one file for each contig
        #readAlignment
        pileup = "%s/out/pileup/%s" % (self.rootPath, self.ko)
        for contigID in self.contigList:
            with open('%s/%s-%s' % (pileup, self.ko, contigID), 'w') as f:
                    #print original contig and full sequence
                    f.write(">%s\n" % contigID)
                    f.write(str(self.contigList[contigID]['fullseq']+"\n"))
                    #print the reads in order
                    if contigID in poshash:
                        for key, info in sorted(poshash[contigID].items()):
                            f.write(">%s-%s\n%s\n" % (info[0][0], contigID, info[0][1]))
        #fixAlignment
        for contigID in self.contigList:
            msa = fixAlignment(self.rootPath, self.ko, contigID)
            newOutputFile = '%s/%s-%s.reAligned.msa' % (pileup, self.ko, contigID)
            with open(newOutputFile, 'w') as newOut:
                    #print original contig and full sequence
                    newOut.write(">%s\n" % contigID)
                    newOut.write(str(msa[contigID]) + "\n")
                    #print the reads in order
                    if contigID in poshash:
                        for key, info in sorted(poshash[contigID].items()):
                            newOut.write(">%s-%s\n%s\n" % (info[0][0], contigID, msa[info[0][0]]))

    def __readStatus(self):
        """
        stores read alignment information for use later
        """
        filePath = "%s/out/newbler/%s/454ReadStatus.txt"%(self.rootPath, self.ko)
        df = pd.read_csv(filePath, sep="\t")
        df.columns=['Accno','ReadStatus','5Contig','5Position','5Strand','3Contig','3Position','3Strand']

        """
        Accno   Read Status     5' Contig       5' Position     5' Strand       3' Contig       3' Position     3' Strand
        simuREAD_62|taxID|191767|loc|7076959-7077060|outpu      Assembled       contig00200     258     -       contig00200     157     +
        simuREAD_332|taxID|201096|loc|478844-478945|output      Assembled       contig02440     104     -       contig02440     8       +
        simuREAD_883|taxID|18|loc|841131-841232|output|s_1      Assembled       contig00107     211     +       contig00107     312     -
        simuREAD_2334|taxID|561|loc|4092117-4092218|output      Assembled       contig00767     304     +       contig00767     404     -
        """

        def splitNstore(row):
            readID = row['Accno'].split("|")[0]
            isAss = row['ReadStatus'] == 'Assembled'
            if isAss:
                isSame = row['5Contig'] == row['3Contig']
                isPos = row['5Strand'] == '+'
                if isSame :
                    self.readInfo[readID] = {
                        'parent': row['5Contig'],
                        'startPos':  int(row['5Position']) if isPos else int(row['3Position']),
                        'endPos'      : int(row['3Position']) if isPos else int(row['5Position']),
                        'direction':  'forward' if isPos else 'reverse'
                    }

        df.apply(splitNstore, axis=1)


