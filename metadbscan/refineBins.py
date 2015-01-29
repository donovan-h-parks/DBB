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
import re
import itertools
import multiprocessing as mp
from collections import defaultdict

from common import checkFileExists
from defaultValues import DefaultValues
from seqUtils import readSeqStats
from distributions import Distributions, readDistributions

from resolveConflicts import ResolveConflicts
from greedy import Greedy, Contig, Core, reportBins, readContigIdtoBinId
from genomicSignatures import GenomicSignatures


class RefineBins(object):
    def __init__(self, bDebug=False):
        self.logger = logging.getLogger()

    def __workerThread(self, unbinnedContigs, cores, distributions, genomicSig, workerQueue, writerQueue):
        """Process each unbinned contig in parallel."""

        # assign unbinned contig to a core bin if and only if:
        #   1) the core is valid: within the defined GC, TD and coverage distribution cutoffs
        #   2) the core is the closest in GC, TD, and coverage space of all valid cores

        while True:
            index = workerQueue.get(block=True, timeout=None)
            if index == None:
                break

            unbinnedContig = unbinnedContigs[index]

            # find closest core in GC, TD, and coverage space
            closestBin = [None] * 3
            minDistances = [1e9] * 3

            for binId, core in cores.iteritems():
                if not distributions.withinDistGC(unbinnedContig, core):
                    continue

                if not distributions.withinDistCov(unbinnedContig, core):
                    continue

                tdDistance = genomicSig.distance(unbinnedContig.tetraSig, core.tetraSig)
                if not distributions.witinDistTD(tdDistance, unbinnedContig):
                    continue

                gcDistance = distributions.gcDistance(unbinnedContig, core)
                covDistance = distributions.covDistance(unbinnedContig, core)

                distances = [gcDistance, tdDistance, covDistance]
                for i, distance in enumerate(distances):
                    if distance < minDistances[i]:
                        minDistances[i] = distance
                        closestBin[i] = binId

            # check if unbinned contig is closest to a single core is all three spaces
            bAssigned = False
            bFailedDistanceTest = False
            bInvalidWithAllCores = False
            if closestBin[0] != None:
                if closestBin[0] == closestBin[1] and closestBin[0] == closestBin[2]:
                    bAssigned = True
                    newBinNum = cores[closestBin[0]].binId
                else:
                    bFailedDistanceTest = True
                    newBinNum = Greedy.UNBINNED
            else:
                bInvalidWithAllCores = True
                newBinNum = Greedy.UNBINNED

            # allow results to be processed or written to file
            writerQueue.put((index, newBinNum, unbinnedContig.length, bAssigned, bFailedDistanceTest, bInvalidWithAllCores))

    def __writerThread(self, newAssignments, numUnbinnedContigs, writerQueue):
        """Store results of worker threads in a single thread."""

        numAssigned = 0
        numInvalidWithAllCores = 0
        numFailedDistanceTest = 0
        assignedBP = 0
        unassignedBP = 0

        processedSeqs = 0
        while True:
            index, newBinNum, unbinnedLen, bAssigned, bFailedDistanceTest, bInvalidWithAllCores = writerQueue.get(block=True, timeout=None)
            if index == None:
                break

            newAssignments[index] = newBinNum

            if self.logger.getEffectiveLevel() <= logging.INFO:
                processedSeqs += 1
                statusStr = '    Processed %d of %d (%.2f%%) unbinned contigs.' % (processedSeqs, numUnbinnedContigs, float(processedSeqs) * 100 / numUnbinnedContigs)
                sys.stdout.write('%s\r' % statusStr)
                sys.stdout.flush()

            # check if unbinned contig is closest to a single core is all three spaces
            if bAssigned:
                    numAssigned += 1
                    assignedBP += unbinnedLen
            elif bFailedDistanceTest:
                    numFailedDistanceTest += 1
                    unassignedBP += unbinnedLen
            elif bInvalidWithAllCores:
                numInvalidWithAllCores += 1
                unassignedBP += unbinnedLen
            else:
                self.logger.error('  [Error] Unknown assignment for unbinned contig.')
                sys.exit()

        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stdout.write('\n')

        numUnassigned = numUnbinnedContigs - numAssigned
        self.logger.info('    Assigned %d of %d (%.2f%%) unbinned contigs to a core bin.' % (numAssigned, numUnbinnedContigs, numAssigned * 100.0 / numUnbinnedContigs))
        self.logger.info('      - %d of %d (%.2f%%) unassigned contigs were invalid with all cores' % (numInvalidWithAllCores, numUnassigned, numInvalidWithAllCores * 100.0 / numUnassigned))
        self.logger.info('      - %d of %d (%.2f%%) unassigned contigs failed distance test' % (numFailedDistanceTest, numUnassigned, numFailedDistanceTest * 100.0 / numUnassigned))
        self.logger.info('    Assigned %.2f of %.2f Mbps (%.2f%%) of the unbinned bases.' % (float(assignedBP) / 1e6, float(assignedBP + unassignedBP) / 1e6, float(assignedBP) * 100 / (assignedBP + unassignedBP)))

    def __refineBins(self, unbinnedContigs, cores, distributions, genomicSig, numThreads):
        """Process each unbinned contig in parallel."""

        # process each unbinned sequence independently
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        # populate worker queue with data to process
        for i in xrange(0, len(unbinnedContigs)):
            workerQueue.put(i)

        for _ in range(numThreads):
            workerQueue.put(None)

        newAssignments = mp.Manager().list([None] * len(unbinnedContigs))

        workerProc = [mp.Process(target=self.__workerThread, args=(unbinnedContigs, cores, distributions, genomicSig, workerQueue, writerQueue)) for _ in range(numThreads)]
        writeProc = mp.Process(target=self.__writerThread, args=(newAssignments, len(unbinnedContigs), writerQueue))

        writeProc.start()

        for p in workerProc:
            p.start()

        for p in workerProc:
            p.join()

        writerQueue.put((None, None, None, None, None, None))
        writeProc.join()

        for i, binNum in enumerate(newAssignments):
            unbinnedContigs[i].binId = binNum

    def __resolveConflictingContigs(self, unbinnedContigs, binnedContigs):
        # map contigs to their originating scaffolds
        scaffoldsToContigs = defaultdict(list)
        for contig in itertools.chain(unbinnedContigs, binnedContigs):
            match = re.search('_c[0-9]*$', contig.contigId)
            if match:
                scaffoldId = contig.contigId[0:match.start()]
                scaffoldsToContigs[scaffoldId].append(contig)

        rc = ResolveConflicts()
        rc.resolve(scaffoldsToContigs, Greedy.UNBINNED)

        if rc.numPartitionedSeqs != 0:
            self.logger.info('')
            self.logger.info('  Resolving scaffolds in conflict: ')
            self.logger.info('    - %d of %d (%.2f%%) partitioned scaffolds completely within a single bin' % (rc.numProperBin, rc.numPartitionedSeqs, rc.numProperBin * 100.0 / rc.numPartitionedSeqs))
            self.logger.info('    - %d of %d (%.2f%%) partitioned scaffolds were completely unbinned' % (rc.numNoise, rc.numPartitionedSeqs, rc.numNoise * 100.0 / rc.numPartitionedSeqs))

            self.logger.info('    - %d of %d (%.2f%%) partitioned scaffolds were partially unbinned' % (rc.numPartialNoise, rc.numPartitionedSeqs, rc.numPartialNoise * 100.0 / rc.numPartitionedSeqs))
            self.logger.info('      - %d (%.2f%%) resolved by reassignment to a single bin' % (rc.numPartialReassignment, rc.numPartialReassignment * 100.0 / max(rc.numPartialNoise, 1)))
            self.logger.info('      - %d (%.2f%%) resolved by marking all partitions as noise' % (rc.numPartialToNoise, rc.numPartialToNoise * 100.0 / max(rc.numPartialNoise, 1)))

            self.logger.info('    - %d of %d (%.2f%%) partitioned scaffolds were assigned to multiple bins' % (rc.numConflicts, rc.numPartitionedSeqs, rc.numConflicts * 100.0 / rc.numPartitionedSeqs))
            self.logger.info('      - %d (%.2f%%) resolved by reassignment to a single bin' % (rc.numConflictsReassignment, rc.numConflictsReassignment * 100.0 / max(rc.numConflicts, 1)))
            self.logger.info('      - %d (%.2f%%) resolved by marking all partitions as noise' % (rc.numConflictsToNoise, rc.numConflictsToNoise * 100.0 / max(rc.numConflicts, 1)))

    def __reportBins(self, refinedBinFile, unbinnedContigs, binnedContigs, argStr):
        # determine contigs in each bin
        binIdToContigs = defaultdict(set)
        for contig in unbinnedContigs:
            binIdToContigs[contig.binId].add(contig)

        for contig in binnedContigs:
            binIdToContigs[contig.binId].add(contig)

        # determine size of each bin in base pairs
        binSizeBP = defaultdict(int)
        totalBinedBP = 0
        totalUnbinedBP = 0
        totalBinedContigs = 0
        for binId in sorted(binIdToContigs):
            binSizeBP = 0
            for contig in binIdToContigs[binId]:
                binSizeBP += contig.length
                if binId == Greedy.UNBINNED:
                    totalUnbinedBP += contig.length
                else:
                    totalBinedBP += contig.length

            if binId != Greedy.UNBINNED:
                self.logger.info('    Contigs in bin %d: %d (%.2f Mbp)' % (binId, len(binIdToContigs[binId]), float(binSizeBP) / 1e6))
                totalBinedContigs += len(binIdToContigs[binId])

        self.logger.info('')
        self.logger.info('    Total binned contigs: %d (%.2f Mbp)' % (totalBinedContigs, float(totalBinedBP) / 1e6))
        self.logger.info('    Total unbinned contigs: %d (%.2f Mbp)' % (len(binIdToContigs[Greedy.UNBINNED]), float(totalUnbinedBP) / 1e6))

        # save bining to file
        fout = open(refinedBinFile, 'w')
        fout.write('# ' + argStr + '\n')
        fout.write('Contig Id\tBin Id\tLength\tGC\tCoverage\n')
        for contigId, contigs in binIdToContigs.iteritems():
            for contig in contigs:
                fout.write(contig.contigId + '\t%d\t%d\t%.3f\t%.3f\n' % (contigId, contig.length, contig.GC, contig.coverage))

        for contig in binIdToContigs[Greedy.UNBINNED]:
            fout.write(contig.contigId + '\t%s\t%d\t%.3f\t%.3f\n' % ('unbinned', contig.length, contig.GC, contig.coverage))

        fout.close()

        self.logger.info('')
        self.logger.info('    Binning information written to: ' + refinedBinFile)

    def run(self, preprocessDir, binningFile, minSeqLen, gcDistPer, tdDistPer, covDistPer, numThreads, refinedBinFile, argStr):
        # verify inputs
        contigStatsFile = os.path.join(preprocessDir, 'contigs.seq_stats.tsv')
        checkFileExists(contigStatsFile)

        contigTetraFile = os.path.join(preprocessDir, 'contigs.tetra.tsv')
        checkFileExists(contigTetraFile)

        # read contig stats
        self.logger.info('  Reading contig statistics.')
        contigStats = readSeqStats(contigStatsFile)

        # read tetranucleotide signatures
        self.logger.info('  Reading contig tetranucleotide signatures.')
        genomicSig = GenomicSignatures(DefaultValues.DEFAULT_KMER_SIZE, 1)
        tetraSigs = genomicSig.read(contigTetraFile)

        # read bin assignments
        self.logger.info('  Reading core bin assignments.')
        contigIdToBinId = readContigIdtoBinId(binningFile)

        # calculate statistics of core bins and unbinned contigs
        self.logger.info('')
        self.logger.info('  Calculating statistics of core bins.')

        cores = {}
        unbinnedContigs = []
        binnedContigs = []
        for contigId, contigStat in contigStats.iteritems():
            binId = contigIdToBinId.get(contigId, Greedy.UNBINNED)

            contigLen, contigGC, contigCov = contigStat
            contigTetraSig = tetraSigs[contigId]

            contig = Contig(contigId, contigLen, contigGC, contigCov, contigTetraSig)
            contig.binId = binId

            if binId == Greedy.UNBINNED:
                if contigLen > minSeqLen:
                    unbinnedContigs.append(contig)
            else:
                binnedContigs.append(contig)

                # build statistics for core bins
                if binId not in cores:
                    cores[binId] = Core(binId, 0, 0, 0, [0] * genomicSig.numKmers())

                core = cores[binId]

                core.length += contigLen
                weight = float(contigLen) / core.length
                core.GC = contigGC * weight + core.GC * (1.0 - weight)
                core.coverage = contigCov * weight + core.coverage * (1.0 - weight)
                for i in xrange(0, len(core.tetraSig)):
                    core.tetraSig[i] = contigTetraSig[i] * weight + core.tetraSig[i] * (1.0 - weight)

        self.logger.info('    Identified %d binned contigs from %d core bins.' % (len(binnedContigs), len(cores)))
        self.logger.info('    Identified %d unbinned contigs >= %d bps.' % (len(unbinnedContigs), minSeqLen))

        # read GC, TD and coverage distributions
        self.logger.info('')
        self.logger.info('  Reading GC, TD, and coverage distributions.')

        gcDist, tdDist, covDist = readDistributions(preprocessDir, gcDistPer, tdDistPer, covDistPer)

        # refine bins
        self.logger.info('')
        self.logger.info('  Refining bins:')
        distributions = Distributions(gcDist, gcDistPer, tdDist, tdDistPer, covDist, covDistPer)
        self.__refineBins(unbinnedContigs, cores, distributions, genomicSig, numThreads)

        # resolve conflicts between contigs originating on the same scaffold
        # and use the conflicts as a measure of the binning quality
        self.__resolveConflictingContigs(unbinnedContigs, binnedContigs)

        # report and save results of binning
        self.logger.info('')
        self.logger.info('  Refined bins:')
        reportBins(refinedBinFile, unbinnedContigs + binnedContigs, argStr)

