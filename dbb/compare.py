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
import logging
import re

from common import checkFileExists
from defaultValues import DefaultValues
from seqUtils import readSeqStats

from greedy import Greedy, Contig, Core, readContigIdtoBinId
from genomicSignatures import GenomicSignatures
from distributions import Distributions, readDistributions

from numpy import mean


class Compare(object):
    def __init__(self, bDebug=False):
        self.logger = logging.getLogger()

    def run(self, preprocessDir, binningFile, queryBinId, queryScaffoldId, gcDistPer, tdDistPer, covDistPer):
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

        queryContigs = []
        queryBinContigs = []
        bins = {}
        for contigId, contigStat in contigStats.iteritems():
            binId = contigIdToBinId.get(contigId, Greedy.UNBINNED)
            contigLen, contigGC, contigCov = contigStat
            contigTetraSig = tetraSigs[contigId]

            scaffoldId = contigId
            match = re.search('_c[0-9]*$', contigId)
            if match:
                scaffoldId = contigId[0:match.start()]

            if scaffoldId == queryScaffoldId:
                contig = Contig(contigId, contigLen, contigGC, contigCov, contigTetraSig)
                contig.binId = binId

                queryContigs.append(contig)

            # calculate bin statistics
            if binId not in bins:
                bins[binId] = Core(binId, 0, 0, 0, [0] * genomicSig.numKmers())

            bin = bins[binId]
            bin.length += contigLen
            weight = float(contigLen) / bin.length
            bin.GC = contigGC * weight + bin.GC * (1.0 - weight)
            bin.coverage = contigCov * weight + bin.coverage * (1.0 - weight)
            for i in xrange(0, len(bin.tetraSig)):
                bin.tetraSig[i] = contigTetraSig[i] * weight + bin.tetraSig[i] * (1.0 - weight)

            if binId == queryBinId:
                queryBinContigs.append(Contig(contigId, contigLen, contigGC, contigCov, contigTetraSig))

        queryBin = bins[queryBinId]

        # calculate mean change in tetra signature
        deltaTetra = []
        for contig in queryBinContigs:
            deltaTetra.append(genomicSig.distance(contig.tetraSig, queryBin.tetraSig))

        meanDeltaTetra = mean(deltaTetra)

        # compare contigs to bin
        print ''
        print 'Contig Id\tGC\tCoverage\tdelta GC\tdelta Coverage\tdelta tetra\tWithin GC\tWithin coverage\tWithin tetra\tClosest GC\tClosest coverage\tClosest tetra'
        for contig in queryContigs:
            tdDistanceQuery = genomicSig.distance(contig.tetraSig, queryBin.tetraSig)

            bWithinGC = distributions.withinDistGC(contig, queryBin)
            bWithinCov = distributions.withinDistCov(contig, queryBin)
            bWithinTetra = distributions.witinDistTD(tdDistanceQuery, contig)

            closestBin = [None] * 3
            minDistances = [1e9] * 3

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

                distances = [gcDistance, covDistance, tdDistance]
                for i, distance in enumerate(distances):
                    if distance < minDistances[i]:
                        minDistances[i] = distance
                        closestBin[i] = binId

            print '%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.2f\t%s\t%s\t%s\t%s\t%s\t%s' % (contig.contigId,
                                                                    contig.GC * 100, contig.coverage, (contig.GC - queryBin.GC) * 100, contig.coverage - queryBin.coverage, tdDistanceQuery,
                                                                    str(bWithinGC), str(bWithinCov), str(bWithinTetra),
                                                                    str(closestBin[0]), str(closestBin[1]), str(closestBin[2]))

        print ''
        print 'Bin Id = %s, mean GC = %.1f, mean coverage = %.1f, mean delta tetra = %.2f' % (queryBinId, queryBin.GC * 100, queryBin.coverage, meanDeltaTetra)
