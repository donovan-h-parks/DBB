#!/usr/bin/env python

###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import os
import sys
import logging
import random
import multiprocessing as mp
from collections import defaultdict
from operator import itemgetter

import pysam
import numpy as np

from genomicSignatures import GenomicSignatures
from common import checkFileExists
from seqUtils import readFasta

class ReadLoader:
    """Callback for counting aligned reads with pysam.fetch"""
    def __init__(self, refLength):
        self.minAlignPer = 0.98
        self.maxEditDistPer = 0.02
        
        self.numReads = 0
        self.numMappedReads = 0
        self.numDuplicates = 0
        self.numSecondary = 0
        self.numFailedQC = 0
        self.numFailedAlignLen = 0
        self.numFailedEditDist = 0
        self.numFailedProperPair = 0

        self.coverage = np.zeros(refLength)

    def __call__(self, read):
        self.numReads += 1

        if read.is_unmapped:
            pass
        elif read.is_duplicate:
            self.numDuplicates += 1
        elif read.is_secondary:
            self.numSecondary += 1
        elif read.is_qcfail:
            self.numFailedQC += 1
        elif read.alen < self.minAlignPer*read.rlen:
            self.numFailedAlignLen += 1
        elif read.opt('NM') > self.maxEditDistPer*read.rlen:
            self.numFailedEditDist += 1
        elif not read.is_proper_pair:
            self.numFailedProperPair += 1
        else:
            self.numMappedReads += 1

            # Note: the alignment length (alen) is used instead of the
            # read length (rlen) as this bring the calculated coverage
            # in line with 'samtools depth' (at least when the min
            # alignment length and edit distance thresholds are zero).
            self.coverage[read.pos:read.pos+read.alen] += 1

class Preprocess(object):
    """ Preprocess contig data in order to calculate:
         1. Sequence statistics for contigs
         2. Sequence statistics for blocks (partitioned contigs)
         3. Coverage distribution
         4. Tetranucleotide signatures
    """
    def __init__(self):
        self.logger = logging.getLogger()
        
    def __calculateGC(self, seq):
        seq = seq.upper()

        gc = at = 0
        for base in seq:
            if base == 'G' or base == 'C':
                gc += 1
            elif base == 'A' or base == 'T':
                at += 1

        return float(gc) / (gc+at)
    
    def __partitions(self, seqLen, blockSize):
        partitions = []
        
        # break contig into partitions of a specific size      
        numBlocks = int(seqLen / blockSize)
        if numBlocks == 0 or numBlocks == 1:
            # sequence is too short to break into partitions
            partitions.append([0, seqLen])
        else:            
            blockStepSize = int(seqLen / numBlocks) 
            
            for blockIndex in xrange(0, numBlocks-1):
                partitions.append([blockStepSize*blockIndex, blockStepSize*(blockIndex+1)])
            
            # add last block and make sure it runs to the end of the sequence
            partitions.append([blockStepSize*(numBlocks-1), seqLen])
            
        return partitions

    def __workerThread(self, bamFile, contigs, genomicSig, minSeqLen, percent, blockSize, queueIn, queueOut):
        """Process each data item in parallel."""
        while True:
            contigIds, contigLens = queueIn.get(block=True, timeout=None)
            if contigIds == None:
                break

            bamfile = pysam.Samfile(bamFile, 'rb')

            for contigId, contigLen in zip(contigIds, contigLens):
                readLoader = ReadLoader(contigLen)
                bamfile.fetch(contigId, 0, contigLen, callback = readLoader)

                contigGC = self.__calculateGC(contigs[contigId])
                contigCoverage = np.mean(readLoader.coverage)
                contigTetraSig = genomicSig.seqSignature(contigs[contigId])
                
                # calculate sequence statistics for contig partitions
                partitions = self.__partitions(contigLen, blockSize)
                
                partitionSeqStats = {}
                if len(partitions) > 1:
                    for partitionIndex, partition in enumerate(partitions):
                        start, end = partition
                        partitionLen = end - start
                        partitionGC = self.__calculateGC(contigs[contigId][start:end])
                        partitionCoverage = np.mean(readLoader.coverage[start:end])
                        tetraSig = genomicSig.seqSignature(contigs[contigId][start:end])
                        partitionSeqStats[contigId + '_p' + str(partitionIndex)] = [partitionLen, partitionGC, partitionCoverage, tetraSig]
                else:
                    partitionSeqStats[contigId] = [contigLen, contigGC, contigCoverage, contigTetraSig]
                    
                # calculate coverage distribution parameters
                covDist = defaultdict(list)
                testLen = minSeqLen
                    
                while contigLen >= 4*testLen:
                    numTestPts = contigLen/testLen
                    
                    perErrors = []
                    for i in xrange(0, numTestPts):
                        testCoverage = np.mean(readLoader.coverage[testLen*i:testLen*(i+1)])
                        perErrors.append( (testCoverage - contigCoverage)*100.0 / contigCoverage)
                        
                    # limit contribution of any sequences to the coverage distribution, and
                    # ensure the number of test points over all sequences is kept to a reasonble number
                    if numTestPts > 1000:
                        covDist[testLen] = random.sample(perErrors, 1000)
                    else:
                        covDist[testLen] = perErrors
                      
                    testLen = int(testLen + testLen*percent)

                # save results for this contig
                contigSeqStats = [contigLen, contigGC, contigCoverage, contigTetraSig]
                readMappingStats = [readLoader.numReads, readLoader.numDuplicates, readLoader.numSecondary, 
                                    readLoader.numFailedQC, readLoader.numFailedAlignLen, readLoader.numFailedEditDist, 
                                    readLoader.numFailedProperPair, readLoader.numMappedReads]
                queueOut.put((contigId, contigSeqStats, partitionSeqStats, readMappingStats, covDist))

            bamfile.close()

    def __writerThread(self, genomicSig, blockSize, covDistFile, contigSeqStatsFile, contigTetraSigFile, partitionSeqStatsFile, partitionTetraSigFile, numRefSeqs, writerQueue):
        """Store or write results of worker threads in a single thread."""
        contigStatsOut = open(contigSeqStatsFile, 'w')
        contigStatsOut.write('Contig Id\tLength\tGC\tCoverage (mean depth)\tMapped reads\n')
        
        contigTetraSigOut = open(contigTetraSigFile, 'w')
        contigTetraSigOut.write('Sequence Id')
        for kmer in genomicSig.canonicalKmerOrder():
            contigTetraSigOut.write('\t' + kmer)
        contigTetraSigOut.write('\n')
        
        partitionStatsOut = open(partitionSeqStatsFile, 'w')
        partitionStatsOut.write('Partition Id\tLength\tGC\tCoverage (mean depth)\n')
        
        partitionTetraSigOut = open(partitionTetraSigFile, 'w')
        partitionTetraSigOut.write('Sequence Id')
        for kmer in genomicSig.canonicalKmerOrder():
            partitionTetraSigOut.write('\t' + kmer)
        partitionTetraSigOut.write('\n')

        totalReads = 0
        totalDuplicates = 0
        totalSecondary = 0
        totalFailedQC = 0
        totalFailedAlignLen = 0
        totalFailedEditDist = 0
        totalFailedProperPair = 0
        totalMappedReads = 0

        processedRefSeqs = 0
        processedPartitions = 0
        
        globalCovDist = defaultdict(list)
        
        while True:
            contigId, contigSeqStats, partitionSeqStats, readMappingStats, covDist = writerQueue.get(block=True, timeout=None)
            if contigId == None:
                break
            
            for testLen, testPts in covDist.iteritems():
                globalCovDist[testLen].extend(testPts)

            if self.logger.getEffectiveLevel() <= logging.INFO:
                processedRefSeqs += 1
                statusStr = '    Processed %d of %d (%.2f%%) contigs.' % (processedRefSeqs, numRefSeqs, float(processedRefSeqs)*100/numRefSeqs)
                sys.stdout.write('%s\r' % statusStr)
                sys.stdout.flush()

                totalReads += readMappingStats[0]
                totalDuplicates += readMappingStats[1]
                totalSecondary += readMappingStats[2]
                totalFailedQC += readMappingStats[3]
                totalFailedAlignLen += readMappingStats[4]
                totalFailedEditDist += readMappingStats[5]
                totalFailedProperPair += readMappingStats[6]
                totalMappedReads += readMappingStats[7]

            contigStatsOut.write(contigId + '\t' + str(contigSeqStats[0]) + '\t' + str(contigSeqStats[1]))
            contigStatsOut.write('\t' + str(contigSeqStats[2]) + '\t' + str(readMappingStats[7]) + '\n')
            
            contigTetraSigOut.write(contigId)
            contigTetraSigOut.write('\t' + '\t'.join(map(str, contigSeqStats[3])))
            contigTetraSigOut.write('\n')
            
            for partitionId, partitionStats in partitionSeqStats.iteritems():
                processedPartitions += 1
                partitionStatsOut.write(partitionId + '\t' + str(partitionStats[0]))
                partitionStatsOut.write('\t' + str(partitionStats[1]) + '\t' + str(partitionStats[2]) + '\n')
                
                partitionTetraSigOut.write(partitionId)
                partitionTetraSigOut.write('\t' + '\t'.join(map(str, partitionStats[3])))
                partitionTetraSigOut.write('\n')

        if self.logger.getEffectiveLevel() <= logging.INFO:
            statusStr = '    Processing %d of %d (%.2f%%) contigs.' % (processedRefSeqs, numRefSeqs, float(processedRefSeqs)*100/numRefSeqs)
            sys.stdout.write('%s\n' % statusStr)
            sys.stdout.flush()

            print ''
            print '  Read mapping statistics:'
            print '    # total reads: %d' % totalReads
            print '      # properly mapped reads: %d (%.1f%%)' % (totalMappedReads, float(totalMappedReads)*100/totalReads)
            print '      # duplicate reads: %d (%.1f%%)' % (totalDuplicates, float(totalDuplicates)*100/totalReads)
            print '      # secondary reads: %d (%.1f%%)' % (totalSecondary, float(totalSecondary)*100/totalReads)
            print '      # reads failing QC: %d (%.1f%%)' % (totalFailedQC, float(totalFailedQC)*100/totalReads)
            print '      # reads failing alignment length: %d (%.1f%%)' % (totalFailedAlignLen, float(totalFailedAlignLen)*100/totalReads)
            print '      # reads failing edit distance: %d (%.1f%%)' % (totalFailedEditDist, float(totalFailedEditDist)*100/totalReads)
            print '      # reads not properly paired: %d (%.1f%%)' % (totalFailedProperPair, float(totalFailedProperPair)*100/totalReads)
            
            print ''
            print '  Contigs longer than %d bp were partitioned.' % blockSize
            print '    Total sequences after partitioning: %d' % processedPartitions

        contigStatsOut.close()
        contigTetraSigOut.close()
        partitionStatsOut.close()
        partitionTetraSigOut.close()
        
        self.__calculateCoverageDistribution(globalCovDist, covDistFile)

    def __calculateCoverageDistribution(self, covDist, covDistFile):
        # calculate coverage distribution parameters for all test lengths
        # with more than 100 test points
        self.logger.info('')
        self.logger.info('  Calculating coverage distribution parameters.')
        t = {}
        for testLen in covDist.keys():
            testPts = covDist[testLen]
            if len(testPts) > 100:
                t[testLen] = {}
                for percentile in np.arange(0, 100+0.5, 0.5):
                    t[testLen][percentile] = np.percentile(testPts, percentile)
                    
        covDistOut = open(covDistFile, 'w')
        covDistOut.write(str(t))  
        covDistOut.close()
        
        self.logger.info('    Min. distribution length: %d, max. distribution length: %d' % (min(t.keys()), max(t.keys())))
        
    def run(self, contigFile, bamFile, minSeqLen, percent, blockSize, numThreads, outputDir, argsStr):
        checkFileExists(contigFile)
        checkFileExists(bamFile)
        
        # read contigs
        self.logger.info('  Reading contigs.')
        contigs = readFasta(contigFile)
        self.logger.info('    Contigs: %d' % len(contigs))

        # process reference sequences in parallel
        self.logger.info('')
        self.logger.info('  Reading BAM file.')
        
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        bamfile = pysam.Samfile(bamFile, 'rb')
        refSeqIds = bamfile.references
        refSeqLens = bamfile.lengths
        
        filteredRefSeqs = []
        for refSeqId, refLen in zip(refSeqIds, refSeqLens):
            if refLen >= minSeqLen:
                filteredRefSeqs.append([refSeqId, refLen])
                
        self.logger.info('    Reference contigs > %d bps: %d' % (minSeqLen, len(filteredRefSeqs)))

        # populate each thread with reference sequence to process
        # Note: to distribute the work load in a reasonably even manner
        # reference sequences are added in desending order of size and 
        # to the reference list with the fewest total base pairs
        refSeqLists = [[] for _ in range(numThreads)]
        refLenLists = [[] for _ in range(numThreads)]

        filteredRefSeqs.sort(key=itemgetter(1), reverse=True)   # sort in descending order of size
        
        totalRefSeqBP = [0]*len(refSeqLists)
        for refSeqId, refLen in filteredRefSeqs:
            threadIndex = np.argsort(totalRefSeqBP)[0]
            refSeqLists[threadIndex].append(refSeqId)
            refLenLists[threadIndex].append(refLen)
            totalRefSeqBP[threadIndex] += refLen

        for i in xrange(0, len(refSeqLists)):
            workerQueue.put((refSeqLists[i], refLenLists[i]))

        for _ in range(numThreads):
            workerQueue.put((None, None))

        gs = GenomicSignatures(4, 1)

        self.logger.info('')
        self.logger.info('  Processing reference contigs.')
        workerProc = [mp.Process(target = self.__workerThread, args = (bamFile, contigs, gs, minSeqLen, percent, blockSize, workerQueue, writerQueue)) for _ in range(numThreads)]
        
        covDistFile = os.path.join(outputDir, 'coverage_dist.txt') 
        contigSeqStatFile = os.path.join(outputDir, 'contigs.seq_stats.tsv')
        contigTetraSigFile = os.path.join(outputDir, 'contigs.tetra.tsv') 
        partitionSeqStatFile = os.path.join(outputDir, 'partitions.seq_stats.tsv')  
        partitionTetraSigFile = os.path.join(outputDir, 'partitions.tetra.tsv')   
        
        writeProc = mp.Process(target = self.__writerThread, args = (gs, blockSize, covDistFile, contigSeqStatFile, contigTetraSigFile, partitionSeqStatFile, partitionTetraSigFile, len(filteredRefSeqs), writerQueue))

        writeProc.start()

        for p in workerProc:
            p.start()

        for p in workerProc:
            p.join()

        writerQueue.put((None, None, None, None, None))
        writeProc.join()
        
        # report calculate parameters and statistics
        self.logger.info('')
        self.logger.info('  Contig statistics written to: ' + contigSeqStatFile)
        self.logger.info('  Tetranucleotide signatures for contigs written to: ' + contigTetraSigFile)
        self.logger.info('  Partition statistics written to: ' + partitionSeqStatFile)
        self.logger.info('  Tetranucleotide signatures for partitions written to: ' + partitionTetraSigFile)
        self.logger.info('  Coverage distribution parameters written to: ' + covDistFile)
        
        # write command line arguments to file
        # write arguments to file
        parameterFile = os.path.join(outputDir, 'parameters.txt')
        fout = open(parameterFile, 'w')
        fout.write(argsStr)
        fout.close()