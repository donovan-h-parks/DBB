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
from collections import defaultdict

from common import checkFileExists
from defaultValues import DefaultValues
from seqUtils import readSeqStats

from greedy import Greedy, Contig, Core, readContigIdtoBinId
from genomicSignatures import GenomicSignatures
from distributions import Distributions, readDistributions


class Unbinned(object):
    def __init__(self, bDebug=False):
        self.logger = logging.getLogger()

    def run(self, preprocessDir, binningFile, gcDistPer, tdDistPer, covDistPer, minScaffoldLen, minContigLen, allScaffolds, outputFile):
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

        # read distributions
        self.logger.info('  Reading GC, TD, and coverage distributions.')

        gcDist, tdDist, covDist = readDistributions(preprocessDir, gcDistPer, tdDistPer, covDistPer)
        distributions = Distributions(gcDist, gcDistPer, tdDist, tdDistPer, covDist, covDistPer)

        # calculate statistics of core bins and unbinned contigs
        self.logger.info('')
        self.logger.info('  Calculating statistics of bins.')

        bins = {}
        unbinnedContigs = defaultdict(list)
        for contigId, contigStat in contigStats.iteritems():
            binId = contigIdToBinId.get(contigId, Greedy.UNBINNED)

            scaffoldId = contigId
            match = re.search('_c[0-9]*$', contigId)
            if match:
                scaffoldId = contigId[0:match.start()]

            contigLen, contigGC, contigCov = contigStat
            contigTetraSig = tetraSigs[contigId]

            contig = Contig(contigId, contigLen, contigGC, contigCov, contigTetraSig)
            contig.binId = binId

            if contigLen < minContigLen:
                continue

            if binId == Greedy.UNBINNED:
                unbinnedContigs[scaffoldId].append(contig)
            else:
                # build statistics for core bins
                if binId not in bins:
                    bins[binId] = Core(binId, 0, 0, 0, [0] * genomicSig.numKmers())

                bin = bins[binId]

                bin.length += contigLen
                weight = float(contigLen) / bin.length
                bin.GC = contigGC * weight + bin.GC * (1.0 - weight)
                bin.coverage = contigCov * weight + bin.coverage * (1.0 - weight)
                for i in xrange(0, len(bin.tetraSig)):
                    bin.tetraSig[i] = contigTetraSig[i] * weight + bin.tetraSig[i] * (1.0 - weight)

                if allScaffolds:
                    unbinnedContigs[scaffoldId].append(contig)

        # compare contigs to bin
        processedScaffolds = 0

        fout = open(outputFile, 'w')
        fout.write('Scaffold Id\tLength\tGC\tCoverage\tBin Id\tGC\tCoverage\tContigs in Scaffold\tWithin GC\tWithin coverage\tWithin tetra\tClosest GC\tClosest coverage\tClosest tetra\tBinning score\n')
        for scaffoldId, contigs in unbinnedContigs.iteritems():
            if self.logger.getEffectiveLevel() <= logging.INFO:
                processedScaffolds += 1
                statusStr = '    Processed %d of %d (%.2f%%) unbinned scaffolds.' % (processedScaffolds, len(unbinnedContigs), float(processedScaffolds) * 100 / len(unbinnedContigs))
                sys.stdout.write('%s\r' % statusStr)
                sys.stdout.flush()

            for queryBinId, queryBin in bins.iteritems():
                numWithinGC = 0
                numWithinCov = 0
                numWithinTetra = 0
                closestGC = 0
                closestCov = 0
                closestTetra = 0
                binningScore = 0

                scaffoldGC = 0
                scaffoldCov = 0
                scaffoldLen = 0
                for contig in contigs:
                    scaffoldLen += contig.length
                    weight = float(contig.length) / scaffoldLen
                    scaffoldGC = contig.GC * weight + scaffoldGC * (1.0 - weight)
                    scaffoldCov = contig.coverage * weight + scaffoldCov * (1.0 - weight)

                    tdDistanceQuery = genomicSig.distance(contig.tetraSig, queryBin.tetraSig)

                    bInDists = True
                    if distributions.withinDistGC(contig, queryBin):
                        numWithinGC += 1
                    else:
                        bInDists = False

                    if distributions.withinDistCov(contig, queryBin):
                        numWithinCov += 1
                    else:
                        bInDists = False

                    if distributions.witinDistTD(tdDistanceQuery, contig):
                        numWithinTetra += 1
                    else:
                        bInDists = False

                    closestBin = [None] * 3
                    minDistances = [1e9] * 3

                    if bInDists:
                        for binId, bin in bins.iteritems():
                            if not distributions.withinDistGC(contig, bin):
                                continue

                            if not distributions.withinDistCov(contig, bin):
                                continue

                            tdDistance = genomicSig.distance(contig.tetraSig, bin.tetraSig)
                            if not distributions.witinDistTD(tdDistance, contig):
                                continue

                            gcDistance = distributions.gcDistance(contig, bin)
                            covDistance = distributions.covDistance(contig, bin)

                            distances = [gcDistance, tdDistance, covDistance]
                            for i, distance in enumerate(distances):
                                if distance < minDistances[i]:
                                    minDistances[i] = distance
                                    closestBin[i] = binId

                    if closestBin[0] == queryBinId:
                        closestGC += 1
                    if closestBin[1] == queryBinId:
                        closestCov += 1
                    if closestBin[2] == queryBinId:
                        closestTetra += 1

                    if closestBin[1] == queryBinId and closestBin[2] == queryBinId:
                        binningScore += 1

                if scaffoldLen >= minScaffoldLen:
                    fout.write('%s\t%d\t%.1f\t%.1f\t%s\t%.1f\t%.1f\t%d' % (scaffoldId, scaffoldLen, scaffoldGC * 100, scaffoldCov, queryBinId, queryBin.GC * 100, queryBin.coverage, len(contigs)))
                    fout.write('\t%d\t%d\t%d\t%d\t%d\t%d' % (numWithinGC, numWithinCov, numWithinTetra, closestGC, closestCov, closestTetra))
                    fout.write('\t%.2f\n' % (float(binningScore) / len(contigs)))

        fout.close()

        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stdout.write('\n')
