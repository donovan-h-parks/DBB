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
import random
import logging
import multiprocessing as mp
from collections import defaultdict
from operator import itemgetter

import pysam
import numpy as np

from genomicSignatures import GenomicSignatures
from common import checkFileExists
from defaultValues import DefaultValues
from seqUtils import readSeq, baseCount
from splitScaffolds import SplitScaffolds


class ReadLoader:
    """Callback for counting aligned reads with pysam.fetch"""
    def __init__(self, refLength, minAlignPer, maxEditDistPer, bAllowImproperPairs):
        self.minAlignPer = minAlignPer
        self.maxEditDistPer = maxEditDistPer

        self.numReads = 0
        self.numMappedReads = 0
        self.numDuplicates = 0
        self.numSecondary = 0
        self.numFailedQC = 0
        self.numFailedAlignLen = 0
        self.numFailedEditDist = 0
        self.numFailedProperPair = 0

        self.bAllowImproperPairs = bAllowImproperPairs

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
        elif read.alen < self.minAlignPer * read.rlen:
            self.numFailedAlignLen += 1
        elif read.opt('NM') > self.maxEditDistPer * read.rlen:
            self.numFailedEditDist += 1
        elif not self.bAllowImproperPairs and not read.is_proper_pair:
            self.numFailedProperPair += 1
        else:
            self.numMappedReads += 1

            # Note: the alignment length (alen) is used instead of the
            # read length (rlen) as this bring the calculated coverage
            # in line with 'samtools depth' (at least when the min
            # alignment length and edit distance thresholds are zero).
            self.coverage[read.pos:read.pos + read.alen] += 1


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
        a, c, g, t = baseCount(seq)

        gc = g + c
        at = a + t

        return float(gc) / (gc + at)

    def __workerThread(self, bamFile, scaffolds, genomicSig, minSeqLen, percent, minN, minAlignLen, maxEditDist, bAllowImproperPairs, coverageType, queueIn, queueOut):
        """Process each data item in parallel."""

        splitScaffolds = SplitScaffolds()

        while True:
            scaffoldIds, scaffoldLens = queueIn.get(block=True, timeout=None)
            if scaffoldIds == None:
                break

            bamfile = pysam.Samfile(bamFile, 'rb')

            for scaffoldId, scaffoldLen in zip(scaffoldIds, scaffoldLens):
                readLoader = ReadLoader(scaffoldLen, minAlignLen, maxEditDist, bAllowImproperPairs)
                bamfile.fetch(scaffoldId, 0, scaffoldLen, callback=readLoader)

                # identify contigs in scaffold and calculate contig statistics
                scaffoldSeq = scaffolds[scaffoldId]
                contigs = splitScaffolds.splits(scaffoldId, scaffoldSeq, minN)

                # determine nucleotide mask for scaffold
                ntMask = []
                for contigId, split in contigs.iteritems():
                    ntMask += range(split[0], split[1])

                # calculate scaffold statistics
                scaffoldGC = self.__calculateGC(scaffoldSeq)
                if coverageType == 'mean':
                    scaffoldCoverage = np.mean(readLoader.coverage[ntMask])
                elif coverageType == 'trimmed_mean':
                    sortedValues = sorted(readLoader.coverage[ntMask])
                    lowerIndex = int(0.1 * len(sortedValues)) + 1
                    upperIndex = len(sortedValues) - int(0.1 * len(sortedValues)) - 1
                    scaffoldCoverage = np.mean(sortedValues[lowerIndex:upperIndex])
                else:
                    # unknown coverage type
                    self.logger.error('Unknown coverage type.')
                    sys.exit(-1)

                scaffoldTetraSig = genomicSig.seqSignature(scaffoldSeq)
                scaffoldSeqStats = [scaffoldLen, scaffoldGC, scaffoldCoverage, scaffoldTetraSig]

                # calculate contig statistics
                contigSeqStats = {}
                if len(contigs) > 1:
                    for contigId, split in contigs.iteritems():
                        contigLen = split[1] - split[0]
                        contigGC = self.__calculateGC(scaffoldSeq[split[0]:split[1]])
                        if coverageType == 'mean':
                            contigCoverage = np.mean(readLoader.coverage[split[0]:split[1]])
                        else:
                            sortedValues = sorted(readLoader.coverage[split[0]:split[1]])
                            lowerIndex = int(0.1 * len(sortedValues)) + 1
                            upperIndex = len(sortedValues) - int(0.1 * len(sortedValues)) - 1
                            contigCoverage = np.mean(sortedValues[lowerIndex:upperIndex])

                        tetraSig = genomicSig.seqSignature(scaffoldSeq[split[0]:split[1]])
                        contigSeqStats[contigId] = [contigLen, contigGC, contigCoverage, tetraSig, split[0], split[1]]
                else:
                    contigSeqStats[scaffoldId] = scaffoldSeqStats + [0, scaffoldSeq]

                # calculate coverage distribution parameters
                covDist = defaultdict(list)
                testLen = minSeqLen

                for contigId, contigStats in contigSeqStats.iteritems():
                    contigLen, contigGC, contigCoverage, tetraSig, start, _ = contigStats

                    if contigCoverage == 0:
                        continue

                    while contigLen >= 4 * testLen:
                        numTestPts = contigLen / testLen

                        perErrors = []
                        for i in xrange(0, numTestPts):
                            if coverageType == 'mean':
                                testCoverage = np.mean(readLoader.coverage[testLen * i + start:testLen * (i + 1) + start])
                            else:
                                sortedValues = sorted(readLoader.coverage[testLen * i + start:testLen * (i + 1) + start])
                                lowerIndex = int(0.1 * len(sortedValues)) + 1
                                upperIndex = len(sortedValues) - int(0.1 * len(sortedValues)) - 1
                                testCoverage = np.mean(sortedValues[lowerIndex:upperIndex])

                            testCoverage = np.mean(readLoader.coverage[testLen * i + start:testLen * (i + 1) + start])
                            perErrors.append((testCoverage - contigCoverage) * 100.0 / contigCoverage)

                        # limit contribution of any sequences to the coverage distribution, and
                        # ensure the number of test points over all sequences is kept to a reasonable number
                        if numTestPts > 1000:
                            covDist[testLen] = random.sample(perErrors, 1000)
                        else:
                            covDist[testLen] = perErrors

                        testLen = int(testLen + testLen * percent)

                # report mapping statistics
                readMappingStats = [readLoader.numReads, readLoader.numDuplicates, readLoader.numSecondary,
                                    readLoader.numFailedQC, readLoader.numFailedAlignLen, readLoader.numFailedEditDist,
                                    readLoader.numFailedProperPair, readLoader.numMappedReads]

                queueOut.put((scaffoldId, scaffoldSeqStats, contigSeqStats, readMappingStats, covDist))

            bamfile.close()

    def __writerThread(self, genomicSig, scaffoldSeqStatsFile, scaffoldTetraSigFile, contigSeqStatsFile, contigTetraSigFile, covDistFile, numRefSeqs, bAllowImproperPairs, writerQueue):
        """Store or write results of worker threads in a single thread."""

        scaffoldStatsOut = open(scaffoldSeqStatsFile, 'w')
        scaffoldStatsOut.write('Scaffold Id\tLength\tGC\tCoverage (mean depth)\tMapped reads\n')

        scaffoldTetraSigOut = open(scaffoldTetraSigFile, 'w')
        scaffoldTetraSigOut.write('Scaffold Id')
        for kmer in genomicSig.canonicalKmerOrder():
            scaffoldTetraSigOut.write('\t' + kmer)
        scaffoldTetraSigOut.write('\n')

        contigStatsOut = open(contigSeqStatsFile, 'w')
        contigStatsOut.write('Contig Id\tLength\tGC\tCoverage (mean depth)\n')

        contigTetraSigOut = open(contigTetraSigFile, 'w')
        contigTetraSigOut.write('Contig Id')
        for kmer in genomicSig.canonicalKmerOrder():
            contigTetraSigOut.write('\t' + kmer)
        contigTetraSigOut.write('\n')

        totalReads = 0
        totalDuplicates = 0
        totalSecondary = 0
        totalFailedQC = 0
        totalFailedAlignLen = 0
        totalFailedEditDist = 0
        totalFailedProperPair = 0
        totalMappedReads = 0

        processedRefSeqs = 0
        processedContigs = 0
        contigSizeDist = defaultdict(int)

        globalCovDist = defaultdict(list)

        while True:
            scaffoldId, scaffoldSeqStats, contigSeqStats, readMappingStats, covDist = writerQueue.get(block=True, timeout=None)
            if scaffoldId == None:
                break

            for testLen, testPts in covDist.iteritems():
                globalCovDist[testLen].extend(testPts)

            if self.logger.getEffectiveLevel() <= logging.INFO:
                processedRefSeqs += 1
                statusStr = '    Processed %d of %d (%.2f%%) scaffolds.' % (processedRefSeqs, numRefSeqs, float(processedRefSeqs) * 100 / numRefSeqs)
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

            scaffoldStatsOut.write(scaffoldId + '\t' + str(scaffoldSeqStats[0]) + '\t' + str(scaffoldSeqStats[1]))
            scaffoldStatsOut.write('\t' + str(scaffoldSeqStats[2]) + '\t' + str(readMappingStats[7]) + '\n')

            scaffoldTetraSigOut.write(scaffoldId)
            scaffoldTetraSigOut.write('\t' + '\t'.join(map(str, scaffoldSeqStats[3])))
            scaffoldTetraSigOut.write('\n')

            for contigId, contigStats in contigSeqStats.iteritems():
                processedContigs += 1
                contigStatsOut.write(contigId + '\t' + str(contigStats[0]))
                contigStatsOut.write('\t' + str(contigStats[1]) + '\t' + str(contigStats[2]) + '\n')

                if contigStats[0] >= 10000:
                    contigSizeDist[10000] += 1
                if contigStats[0] >= 5000:
                    contigSizeDist[5000] += 1
                if contigStats[0] >= 2000:
                    contigSizeDist[2000] += 1
                if contigStats[0] >= 1000:
                    contigSizeDist[1000] += 1

                contigTetraSigOut.write(contigId)
                contigTetraSigOut.write('\t' + '\t'.join(map(str, contigStats[3])))
                contigTetraSigOut.write('\n')

        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stdout.write('\n')

            print ''
            print '  Read mapping statistics:'
            print '    # total reads: %d' % totalReads
            print '      # properly mapped reads: %d (%.1f%%)' % (totalMappedReads, float(totalMappedReads) * 100 / totalReads)
            print '      # duplicate reads: %d (%.1f%%)' % (totalDuplicates, float(totalDuplicates) * 100 / totalReads)
            print '      # secondary reads: %d (%.1f%%)' % (totalSecondary, float(totalSecondary) * 100 / totalReads)
            print '      # reads failing QC: %d (%.1f%%)' % (totalFailedQC, float(totalFailedQC) * 100 / totalReads)
            print '      # reads failing alignment length: %d (%.1f%%)' % (totalFailedAlignLen, float(totalFailedAlignLen) * 100 / totalReads)
            print '      # reads failing edit distance: %d (%.1f%%)' % (totalFailedEditDist, float(totalFailedEditDist) * 100 / totalReads)

            if not bAllowImproperPairs:
                print '      # reads not properly paired: %d (%.1f%%)' % (totalFailedProperPair, float(totalFailedProperPair) * 100 / totalReads)

            print ''
            print '  Total contigs within scaffolds: %d' % processedContigs
            for contigLen in sorted(contigSizeDist.keys()):
                print '    Contigs >= %d bp: %d' % (contigLen, contigSizeDist[contigLen])

        scaffoldStatsOut.close()
        scaffoldTetraSigOut.close()
        contigStatsOut.close()
        contigTetraSigOut.close()

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

                qTest = np.arange(0, 100 + 0.5, 0.5).tolist()
                percentiles = np.percentile(testPts, qTest)
                for q, p in zip(qTest, percentiles):
                    t[testLen][q] = p

        covDistOut = open(covDistFile, 'w')
        covDistOut.write(str(t))
        covDistOut.close()

        self.logger.info('    Min. distribution length: %d, max. distribution length: %d' % (min(t.keys()), max(t.keys())))

    def run(self, scaffoldFile, bamFile, minSeqLen, percent, minN, minAlignLen, maxEditDist, bAllowImproperPairs, coverageType, numThreads, outputDir, argsStr):
        checkFileExists(scaffoldFile)
        checkFileExists(bamFile)

        # read scaffolds
        self.logger.info('  Reading scaffolds.')
        scaffolds = {}
        numScaffolds = 0
        for seqId, seq in readSeq(scaffoldFile):
            numScaffolds += 1
            if len(seq) >= minSeqLen:
                scaffolds[seqId] = seq
        self.logger.info('    Scaffolds: %d' % numScaffolds)

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

        self.logger.info('    Reference scaffolds >= %d bps: %d' % (minSeqLen, len(filteredRefSeqs)))

        # populate each thread with reference sequence to process
        # Note: to distribute the work load in a reasonably even manner
        # reference sequences are added in descending order of size and
        # to the reference list with the fewest total base pairs
        refSeqLists = [[] for _ in range(numThreads)]
        refLenLists = [[] for _ in range(numThreads)]

        filteredRefSeqs.sort(key=itemgetter(1), reverse=True)  # sort in descending order of size

        totalRefSeqBP = [0] * len(refSeqLists)
        for refSeqId, refLen in filteredRefSeqs:
            threadIndex = np.argsort(totalRefSeqBP)[0]
            refSeqLists[threadIndex].append(refSeqId)
            refLenLists[threadIndex].append(refLen)
            totalRefSeqBP[threadIndex] += refLen

        for i in xrange(0, len(refSeqLists)):
            workerQueue.put((refSeqLists[i], refLenLists[i]))

        for _ in range(numThreads):
            workerQueue.put((None, None))

        gs = GenomicSignatures(DefaultValues.DEFAULT_KMER_SIZE, 1)

        try:
            self.logger.info('')
            self.logger.info('  Processing reference scaffolds.')
            workerProc = [mp.Process(target=self.__workerThread, args=(bamFile, scaffolds, gs, minSeqLen, percent, minN, minAlignLen, maxEditDist, bAllowImproperPairs, coverageType, workerQueue, writerQueue)) for _ in range(numThreads)]

            scaffoldSeqStatFile = os.path.join(outputDir, 'scaffolds.seq_stats.tsv')
            scaffoldTetraSigFile = os.path.join(outputDir, 'scaffolds.tetra.tsv')
            contigSeqStatFile = os.path.join(outputDir, 'contigs.seq_stats.tsv')
            contigTetraSigFile = os.path.join(outputDir, 'contigs.tetra.tsv')
            covDistFile = os.path.join(outputDir, 'coverage_dist.txt')

            writeProc = mp.Process(target=self.__writerThread, args=(gs, scaffoldSeqStatFile, scaffoldTetraSigFile, contigSeqStatFile, contigTetraSigFile, covDistFile, len(filteredRefSeqs), bAllowImproperPairs, writerQueue))

            writeProc.start()

            for p in workerProc:
                p.start()

            for p in workerProc:
                p.join()

            writerQueue.put((None, None, None, None, None))
            writeProc.join()
        except:
            # make sure all processes are terminated
            for p in workerProc:
                p.terminate()

            writeProc.terminate()

        # write command line arguments to file
        parameterFile = os.path.join(outputDir, 'parameters.txt')
        fout = open(parameterFile, 'w')
        fout.write(argsStr)
        fout.close()

        # report calculate parameters and statistics
        self.logger.info('')
        self.logger.info('  Scaffold statistics written to: ' + scaffoldSeqStatFile)
        self.logger.info('  Tetranucleotide signatures for scaffolds written to: ' + scaffoldTetraSigFile)
        self.logger.info('  Contig statistics written to: ' + contigSeqStatFile)
        self.logger.info('  Tetranucleotide signatures for contigs written to: ' + contigTetraSigFile)
        self.logger.info('  Coverage distribution parameters written to: ' + covDistFile)
        self.logger.info('  Preprocessing parameters written to: ' + parameterFile)
