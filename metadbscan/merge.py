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
from collections import defaultdict

from common import checkFileExists
from defaultValues import DefaultValues
from seqUtils import readSeqStats

from greedy import Greedy, Contig, Core, readContigIdtoBinId
from genomicSignatures import GenomicSignatures
from distributions import Distributions, readDistributions

class Merge(object):
    def __init__(self, bDebug=False):
        self.logger = logging.getLogger()
               
    def run(self, preprocessDir, binningFile, gcDistPer, tdDistPer, covDistPer, outputFile):
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
        binnedContigs = defaultdict(list)
        for contigId, contigStat in contigStats.iteritems():  
            binId = contigIdToBinId.get(contigId, Greedy.UNBINNED)
            if binId != Greedy.UNBINNED: 
                contigLen, contigGC, contigCov = contigStat
                contigTetraSig = tetraSigs[contigId]
                
                contig = Contig(contigId, contigLen, contigGC, contigCov, contigTetraSig)
                contig.binId = binId       
            
                binnedContigs[binId].append(contig)
 
                # build statistics for core bins
                if binId not in bins:
                    bins[binId] = Core(binId, 0, 0, 0, [0]*genomicSig.numKmers())
                    
                bin = bins[binId]

                bin.length += contigLen
                weight = float(contigLen)/bin.length
                bin.GC = contigGC*weight + bin.GC*(1.0-weight)
                bin.coverage = contigCov*weight + bin.coverage*(1.0-weight)
                for i in xrange(0, len(bin.tetraSig)):
                    bin.tetraSig[i] = contigTetraSig[i]*weight + bin.tetraSig[i]*(1.0-weight)
                
        # compare contigs to bin
        fout = open(outputFile, 'w')
        fout.write('Src Bin Id\tGC\tCoverage\tDest Bin Id\tGC\tCoverage\tContigs in Src Bin\tWithin GC\tWithin coverage\tWithin tetra\tClosest GC\tClosest coverage\tClosest tetra\tMerging score\n')
        for binId1, bin1 in bins.iteritems():
            contigs = binnedContigs[binId1]
            
            for binId2, bin2 in bins.iteritems():
                if binId1 == binId2:
                    continue
                    
                fout.write('%s\t%.1f\t%.1f\t%s\t%.1f\t%.1f\t%d' % (binId1, bin1.GC*100, bin1.coverage, binId2, bin2.GC*100, bin2.coverage, len(contigs)))
                numWithinGC = 0
                numWithinCov = 0
                numWithinTetra = 0
                closestGC = 0
                closestCov = 0
                closestTetra = 0
                for contig in contigs:
                    tdDistanceQuery = genomicSig.distance(contig.tetraSig, bin2.tetraSig)
                    
                    if distributions.withinDistGC(contig, bin2):
                        numWithinGC += 1
                        
                    if distributions.withinDistCov(contig, bin2):
                        numWithinCov += 1
                        
                    if distributions.witinDistTD(tdDistanceQuery, contig):
                        numWithinTetra += 1
                     
                    closestBin = [None]*3
                    minDistances = [1e9]*3
                    
                    for binId, bin in bins.iteritems():
                        if binId == binId1:
                            continue
                        
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
                                
                    if closestBin[0] == binId2:
                        closestGC += 1
                    if closestBin[1] == binId2:
                        closestCov += 1
                    if closestBin[2] == binId2:
                        closestTetra += 1
                
                fout.write('\t%d\t%d\t%d\t%d\t%d\t%d' % (numWithinGC, numWithinCov, numWithinTetra, closestGC, closestCov, closestTetra))
                fout.write('\t%.2f\n' % (float(closestGC + closestCov + closestTetra) / (3*len(contigs))))
                    
        fout.close()
             