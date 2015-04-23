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
import ast
import logging
import re
from collections import defaultdict

from common import checkFileExistsfrom defaultValues import DefaultValues

from distributions import Distributions

from resolveConflicts import ResolveConflicts
from genomicSignatures import GenomicSignatures

from coverageNeighbours import CoverageNeighbours
from gcNeighbours import GcNeighbours
from tdNeighbours import TdNeighbours

from prettytable import PrettyTable

import numpy as np


def readSeqStatsForBins(binningFile):
    '''Read file with sequence bin and stat information.'''
    seqStatsForBins = defaultdict(list)
    with open(binningFile) as f:
        f.readline()  # parameter line
        f.readline()  # header line
        for line in f:
            lineSplit = line.split('\t')
            binId = Greedy.UNBINNED if lineSplit[1] == 'unbinned' else int(lineSplit[1])
            seqStatsForBins[binId].append([lineSplit[0], int(lineSplit[2]), float(lineSplit[3]), float(lineSplit[4])])

    return seqStatsForBins


def readContigIdtoBinId(binningFile):
    '''Determing bin id of each contig.'''
    contigIdToBinId = {}
    with open(binningFile) as f:
        f.readline()  # parameter line
        f.readline()  # header line
        for line in f:
            lineSplit = line.split('\t')
            binId = Greedy.UNBINNED if lineSplit[1] == 'unbinned' else lineSplit[1]
            contigIdToBinId[lineSplit[0]] = binId

    return contigIdToBinId


def reportBins(binningFile, contigs, argStr):
    ''' Report bins and write to file.'''
    # get contigs in each core
    cores = defaultdict(list)
    totalUnbinedContigs = 0
    totalUnbinedBP = 0
    unbinnedContigs = []
    for contig in contigs:
        if contig.binId != Greedy.UNBINNED:
            cores[contig.binId].append(contig)
        else:
            totalUnbinedContigs += 1
            totalUnbinedBP += contig.length
            unbinnedContigs.append(contig)

    # report cores
    logger = logging.getLogger()
    if logger.getEffectiveLevel() <= logging.INFO:
        table = PrettyTable(['Bin Id', '# contigs', 'Size (Mbps)', 'GC (%)', 'Coverage'])

        # add cores to table
        totalBinedContigs = 0
        totalBinedBP = 0
        for binId, contigs in cores.iteritems():
            length = 0
            GC = 0
            coverage = 0
            for contig in contigs:
                length += contig.length
                weight = float(contig.length) / length
                GC = contig.GC * weight + GC * (1.0 - weight)
                coverage = contig.coverage * weight + coverage * (1.0 - weight)

            table.add_row([binId, len(contigs), float(length) / 1e6, GC * 100, coverage])
            totalBinedContigs += len(contigs)
            totalBinedBP += length

        table.float_format['Size (Mbps)'] = '.2'
        table.float_format['GC (%)'] = '.1'
        table.float_format['Coverage'] = '.1'
        table.align = 'r'
        print table.get_string(sortby='Bin Id', border=False, padding_width=4)

        logger.info('')
        logger.info('    Total binned contigs: %d (%.2f Mbp)' % (totalBinedContigs, float(totalBinedBP) / 1e6))
        logger.info('    Total unbinned contigs: %d (%.2f Mbp)' % (totalUnbinedContigs, float(totalUnbinedBP) / 1e6))

    # save bins to file
    fout = open(binningFile, 'w')
    fout.write('# ' + argStr + '\n')
    fout.write('Contig Id\tBin Id\tLength\tGC\tCoverage\n')
    for binId, contigs in cores.iteritems():
        for contig in contigs:
            fout.write(contig.contigId + '\t%s\t%d\t%.3f\t%.3f\n' % (binId, contig.length, contig.GC, contig.coverage))

    for contig in unbinnedContigs:
        fout.write(contig.contigId + '\t%s\t%d\t%.3f\t%.3f\n' % ('unbinned', contig.length, contig.GC, contig.coverage))

    fout.close()

    logger.info('')
    logger.info('    Binning information written to: ' + binningFile)


class Contig(object):
    def __init__(self, contigId, length, GC, coverage, tetraSig):
        self.binId = Greedy.UNBINNED

        self.contigId = contigId
        self.length = length
        self.GC = GC
        self.coverage = coverage
        self.tetraSig = tetraSig

    def __str__(self):
        return 'Contig id: %s, bin id: %s, length: %d, GC: %.3f, coverage: %.1f' % (self.contigId, self.binId, self.length, self.GC, self.coverage)


class Core(object):
    def __init__(self, binId, length, GC, coverage, tetraSig):
        self.binId = binId

        self.length = length
        self.GC = GC
        self.coverage = coverage
        self.tetraSig = tetraSig

        self.contigs = None

    def __str__(self):
        return 'Bin id: %d, length: %d, GC: %.3f, coverage: %.1f, # contigs: %d' % (self.binId, self.length, self.GC, self.coverage, len(self.contigs))


class Greedy(object):
    UNBINNED = -1

    def __init__(self, bDebug=False):
        self.logger = logging.getLogger()
    def __greedy(self, sortedContigs, minBinSize, buildDistPer, gcDist, tdDist, covDist, genomicSig, numThreads):
        """Greedily bin contigs with closest neighbour within distributional parameters."""

        distributions = Distributions(gcDist, buildDistPer, tdDist, buildDistPer, covDist, buildDistPer)

        for contig in sortedContigs:
            contig.bProcessed = False

        numCores = 0
        cores = []
        for contigIndex, contig in enumerate(sortedContigs):
            if self.logger.getEffectiveLevel() <= logging.INFO:
                statusStr = '    Processed %d of %d (%.2f%%) contigs.' % (contigIndex + 1, len(sortedContigs), float(contigIndex + 1) * 100 / len(sortedContigs))
                sys.stdout.write('%s\r' % statusStr)
                sys.stdout.flush()

            if contig.bProcessed:
                continue

            contig.bProcessed = True
            contig.binId = numCores
            contigsInCurCore = [contig]

            curCore = Core(numCores, contig.length, contig.GC, contig.coverage, contig.tetraSig)

            # add suitable contigs to core
            bCoreExpanded = True
            while bCoreExpanded:
                # find closest contig in coverage space that is within
                # GC, TD, and coverage distributions
                closestContig = None
                closestAbsCovDiff = 1e12
                for curContig in sortedContigs[contigIndex + 1:]:
                    if curContig.bProcessed:
                        continue

                    absCovDiff = abs(curContig.coverage - curCore.coverage)
                    if absCovDiff > closestAbsCovDiff:
                        continue

                    if not distributions.withinDistGC(curContig, curCore):
                        continue

                    if not distributions.withinDistCov(curContig, curCore):
                        continue

                    tdDistance = genomicSig.distance(curContig.tetraSig, curCore.tetraSig)
                    if not distributions.witinDistTD(tdDistance, curContig):
                        continue

                    closestAbsCovDiff = absCovDiff
                    closestContig = curContig

                # check if core should be expanded
                if closestContig != None:
                    # add contig to the core, updating all core parameters
                    bCoreExpanded = True

                    closestContig.bProcessed = True
                    closestContig.binId = numCores
                    contigsInCurCore.append(closestContig)

                    curCore.length += closestContig.length
                    weight = float(closestContig.length) / curCore.length
                    curCore.GC = closestContig.GC * weight + curCore.GC * (1.0 - weight)
                    curCore.coverage = closestContig.coverage * weight + curCore.coverage * (1.0 - weight)
                    for i in xrange(0, len(curCore.tetraSig)):
                        curCore.tetraSig[i] = closestContig.tetraSig[i] * weight + curCore.tetraSig[i] * (1.0 - weight)
                else:
                    # no suitable contigs to add to the core
                    bCoreExpanded = False

                    if curCore.length > minBinSize:
                        numCores += 1
                        curCore.contigs = contigsInCurCore
                        cores.append(curCore)
                    else:
                        for contig in contigsInCurCore:
                            contig.bProcessed = False
                            contig.binId = Greedy.UNBINNED

        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stdout.write('\n')

        return cores

    def __merge(self, cores, mergeDistPer, gcDist, tdDist, covDist, genomicSig, numThreads):
        """Merge similar cores."""

        distributions = Distributions(gcDist, mergeDistPer, tdDist, mergeDistPer, covDist, mergeDistPer)

        for core in cores:
            core.bProcessed = False

        numMergedCores = 0
        mergedCores = []
        for coreIndex, core in enumerate(cores):
            if self.logger.getEffectiveLevel() <= logging.INFO:
                statusStr = '    Processed %d of %d (%.2f%%) cores.' % (coreIndex + 1, len(cores), float(coreIndex + 1) * 100 / len(cores))
                sys.stdout.write('%s\r' % statusStr)
                sys.stdout.flush()

            if core.bProcessed:
                continue

            core.bProcessed = True

            curMergedCore = Core(numMergedCores, core.length, core.GC, core.coverage, core.tetraSig)
            curMergedCore.contigs = core.contigs

            # merge all similar cores
            bCoreMerged = True
            while bCoreMerged:
                # find closest core in coverage space that is within
                # GC, TD, and coverage distributions
                closestCore = None
                closestAbsCovDiff = 1e12
                for curCore in cores[coreIndex + 1:]:
                    if curCore.bProcessed:
                        continue

                    absCovDiff = abs(curCore.coverage - curMergedCore.coverage)
                    if absCovDiff > closestAbsCovDiff:
                        continue

                    if not distributions.withinDistGC(curCore, curMergedCore):
                        continue

                    if not distributions.withinDistCov(curCore, curMergedCore):
                        continue

                    tdDistance = genomicSig.distance(curCore.tetraSig, curMergedCore.tetraSig)
                    if not distributions.witinDistTD(tdDistance, curCore):
                        continue

                    closestAbsCovDiff = absCovDiff
                    closestCore = curCore

                # check if core should be merged
                if closestCore != None:
                    # merge cores and updating all parameters
                    bCoreMerged = True

                    closestCore.bProcessed = True
                    closestCore.binId = numMergedCores
                    curMergedCore.contigs += closestCore.contigs

                    curMergedCore.length += closestCore.length
                    weight = float(closestCore.length) / curMergedCore.length
                    curMergedCore.GC = closestCore.GC * weight + curMergedCore.GC * (1.0 - weight)
                    curMergedCore.coverage = closestCore.coverage * weight + curMergedCore.coverage * (1.0 - weight)
                    for i in xrange(0, len(curCore.tetraSig)):
                        curMergedCore.tetraSig[i] = closestCore.tetraSig[i] * weight + curMergedCore.tetraSig[i] * (1.0 - weight)
                else:
                    # no suitable core to merge
                    bCoreMerged = False
                    numMergedCores += 1
                    mergedCores.append(curMergedCore)

        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stdout.write('\n')

        # assign contigs to cores
        for core in mergedCores:
            for contig in core.contigs:
                contig.binId = core.binId

        return mergedCores

    def __resolveConflictingContigs(self, contigs, gcDist, tdDist, covDist, mergeDistPer):
        # map contigs to their originating scaffolds
        scaffoldsToContigs = defaultdict(list)
        for contig in contigs:
            match = re.search('_c[0-9]*$', contig.contigId)
            if match:
                scaffoldId = contig.contigId[0:match.start()]
                scaffoldsToContigs[scaffoldId].append(contig)

        rc = ResolveConflicts()
        # rc.merge(scaffoldsToContigs, gcDist, tdDist, covDist, mergeDistPer, Greedy.UNBINNED)
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

    def __filterBins(self, contigs, minBinSize, numCores):
        # get contigs in each core
        cores = defaultdict(list)
        for contig in contigs:
            cores[contig.binId].append(contig)

        # determine size of each core in base pairs
        coreSizeBP = defaultdict(int)
        for i in xrange(Greedy.UNBINNED, numCores):
            for contig in cores[i]:
                coreSizeBP[i] += contig.length

        # filter out core of insufficient size and
        # adjust bin numbers so they are still consecutive
        adjCoreNum = 0
        for i in xrange(0, numCores):
            if coreSizeBP[i] < minBinSize:
                for contig in cores[i]:
                    contig.binId = Greedy.UNBINNED
            else:
                for contig in cores[i]:
                    contig.binId = adjCoreNum

                adjCoreNum += 1

        return adjCoreNum

    def __sortContigs(self, contigs, minSeqLen, buildDistPer, gcDist, tdDist, covDist, numThreads):
        """Sort contigs by the size (in base pairs) of their neighbourhood."""

        # calculate neighbours of each contig under a fixed sequence length
        self.logger.info('    Calculating GC neighbours:')
        gcNeighbours = GcNeighbours(minSeqLen)
        gcNeighbours.run(contigs, gcDist, buildDistPer, numThreads)

        self.logger.info('    Calculating TD neighbours:')
        tdNeighbours = TdNeighbours(minSeqLen)
        tdNeighbours.run(contigs, tdDist, buildDistPer, numThreads)

        self.logger.info('    Calculating coverage neighbours:')
        coverageNeighbours = CoverageNeighbours(minSeqLen)
        coverageNeighbours.run(contigs, covDist, buildDistPer, numThreads)

        # determine size of each contigs neighbourhood in base pairs
        neighbourhoodSizes = []
        for contig in contigs:
            covNeighbours = np.where(contig.covNeighbours == 1)[0]
            gcNeighbours = np.where(contig.gcNeighbours == 1)[0]
            tdNeighbours = np.where(contig.tdNeighbours == 1)[0]

            intersectCovGC = np.intersect1d(covNeighbours, gcNeighbours, assume_unique=True)
            neighbourIndices = np.intersect1d(tdNeighbours, intersectCovGC, assume_unique=True)

            totalNeighbourBPs = sum(contigs[i].length for i in neighbourIndices)
            neighbourhoodSizes.append(totalNeighbourBPs)

        # sort contigs by size of neighbourhood
        sortedIndices = np.argsort(neighbourhoodSizes)[::-1]
        contigs = np.array(contigs)
        contigs = contigs[sortedIndices]

    def run(self, preprocessDir, minSeqLen, minBinSize, buildDistPer, mergeDistPer, numThreads, binningFile, argStr):
        # verify inputs
        gcDistFile = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), '..', 'data', 'gc_dist.txt')
        checkFileExists(gcDistFile)

        tdDistFile = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), '..', 'data', 'td_dist.txt')
        checkFileExists(tdDistFile)

        coverageDistFile = os.path.join(preprocessDir, 'coverage_dist.txt')
        checkFileExists(coverageDistFile)

        contigStatsFile = os.path.join(preprocessDir, 'contigs.seq_stats.tsv')
        checkFileExists(contigStatsFile)

        contigTetraFile = os.path.join(preprocessDir, 'contigs.tetra.tsv')
        checkFileExists(contigTetraFile)

        # read tetranucleotide signatures
        self.logger.info('  Reading contig tetranucleotide signatures.')
        genomicSig = GenomicSignatures(DefaultValues.DEFAULT_KMER_SIZE, 1)
        tetraSigs = genomicSig.read(contigTetraFile)

        # create contig
        self.logger.info('  Reading contig statistics.')
        contigs = []
        contigIdsOfInterest = set()
        with open(contigStatsFile) as f:
            f.readline()

            for line in f:
                lineSplit = line.split('\t')

                contigId = lineSplit[0]
                contigLen = int(lineSplit[1])
                GC = float(lineSplit[2])
                coverage = float(lineSplit[3])

                if contigLen > minSeqLen:
                    contig = Contig(contigId, contigLen, GC, coverage, tetraSigs[contigId])
                    contigs.append(contig)
                    contigIdsOfInterest.add(contigId)

        self.logger.info('    Contigs > %d bps: %d' % (minSeqLen, len(contigs)))

        # read GC, TD and coverage distributions
        self.logger.info('')
        self.logger.info('  Reading GC, TD and coverage distributions.')

        with open(gcDistFile, 'r') as f:
            s = f.read()
            gcDist = ast.literal_eval(s)

        with open(tdDistFile, 'r') as f:
            s = f.read()
            tdDist = ast.literal_eval(s)

        with open(coverageDistFile) as f:
            s = f.read()
            covDist = ast.literal_eval(s)

        # sort contigs in ascending order of neighbour size
        self.logger.info('')
        self.logger.info('  Sorting contigs in descending order of neighbourhood size.')
        self.__sortContigs(contigs, minSeqLen, buildDistPer, gcDist, tdDist, covDist, numThreads)

        # perform greedy binning
        self.logger.info('')
        self.logger.info('  Running greedy binning.')
        cores = self.__greedy(contigs, minBinSize, buildDistPer, gcDist, tdDist, covDist, genomicSig, numThreads)
        numOriginalCores = len(cores)
        self.logger.info('    Found %d cores.' % numOriginalCores)

        # sort contigs in ascending order of neighbour size
        self.logger.info('')
        self.logger.info('  Sorting cores in descending order of size.')
        cores.sort(key=lambda x: x.length, reverse=True)

        # merge similar cores
        self.logger.info('')
        self.logger.info('  Merging similar cores.')
        cores = self.__merge(cores, mergeDistPer, gcDist, tdDist, covDist, genomicSig, numThreads)
        numMergedCores = len(cores)
        self.logger.info('    Merged %d cores.' % (numOriginalCores - numMergedCores))

        # resolve conflicts between contigs originating on the same scaffold
        # and use the conflicts as a measure of the binning quality
        self.__resolveConflictingContigs(contigs, gcDist, tdDist, covDist, mergeDistPer)

        # filter bins sufficiently small that they may belong to another bin
        # and from which there is limited data too calculate robust statistics
        self.logger.info('')
        self.logger.info('  Filtering cores with < %d bps.' % minBinSize)
        numFilteredCores = self.__filterBins(contigs, minBinSize, numMergedCores)

        if numMergedCores != 0:
            self.logger.info('    Retained %d of %d (%.2f%%) bins.' % (numFilteredCores, numMergedCores, numFilteredCores * 100.0 / numMergedCores))

        # report and save results of binning
        self.logger.info('')
        self.logger.info('  Core bins:')
        reportBins(binningFile, contigs, argStr)
