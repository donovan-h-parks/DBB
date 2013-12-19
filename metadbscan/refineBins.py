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
import ast
import logging
import re
import itertools
from collections import defaultdict

from common import checkFileExists
from seqUtils import readSeqStats
from distributions import readDistribution, findNearest, inRange

from resolveConflicts import ResolveConflicts
from dbscan import DBSCAN, Sequence
from genomicSignatures import GenomicSignatures

class RefineBins(object):
    def __init__(self, bDebug=False):
        self.logger = logging.getLogger()
        
    def __withinDistGC(self, gcDist, unbinnedContig, core):
        closestMeanGC = findNearest(gcDist.keys(), core.GC)
        closestSeqLen = findNearest(gcDist[closestMeanGC].keys(), unbinnedContig.seqLen)
        lowerBound, upperBound = gcDist[closestMeanGC][closestSeqLen]
        
        gcDiff = unbinnedContig.GC - core.GC
        return inRange(gcDiff, lowerBound, upperBound)
        
    def __witinDistTD(self, tdDist, tetraDist, unbinnedContig):
        return tetraDist < tdDist[findNearest(tdDist.keys(), unbinnedContig.seqLen)]
        
    def __withinDistCov(self, covDist, covDistPer, unbinnedContig, core):
        closestSeqLen = findNearest(covDist.keys(), unbinnedContig.seqLen)
        lowerBound = covDist[closestSeqLen][int((100 - covDistPer)/2)]
        upperBound = covDist[closestSeqLen][int((100 + covDistPer)/2)]
        
        covPerDiff = (unbinnedContig.coverage - core.coverage)*100.0 / core.coverage
        return inRange(covPerDiff, lowerBound, upperBound)
                    
    def __gcDistance(self, unbinnedContig, core):
        return abs(unbinnedContig.GC - core.GC)
        
    def __covDistance(self, unbinnedContig, core):
        return abs(unbinnedContig.coverage - core.coverage)
        
    def __refineBins(self, unbinnedContigs, cores, gcDist, tdDist, covDist, covDistPer, genomicSig):
        # assign unbinned contig to a core bin if and only if:
        #   1) the core is valid: within the defined GC, TD and coverage distribution cutoffs
        #   2) the core is the closest in GC, TD, and coverage space of all valid cores
        #   
   
        numAssigned = 0
        numInvalidWithAllCores = 0
        numFailingDistanceTest = 0
        assignedBP = 0
        unassignedBP = 0
        for processedSeqs, unbinnedContig in enumerate(unbinnedContigs):
            # find closest core in GC, TD, and coverage space
            closestCluster = [None]*3
            minDistances = [1e9]*3
            
            if self.logger.getEffectiveLevel() <= logging.INFO:
                processedSeqs += 1
                statusStr = '    Processed %d of %d (%.2f%%) unbinned contigs.' % (processedSeqs, len(unbinnedContigs), float(processedSeqs)*100/len(unbinnedContigs))
                sys.stdout.write('%s\r' % statusStr)
                sys.stdout.flush()

            for clusterId, core in cores.iteritems():
                if not self.__withinDistGC(gcDist, unbinnedContig, core):
                    continue
                    
                if not self.__withinDistCov(covDist, covDistPer, unbinnedContig, core):
                    continue
                    
                tdDistance = genomicSig.distance(unbinnedContig.tetraSig, core.tetraSig)
                if not self.__witinDistTD(tdDist, tdDistance, unbinnedContig):
                    continue
                    
                gcDistance = self.__gcDistance(unbinnedContig, core)
                covDistance = self.__covDistance(unbinnedContig, core)
                
                distances = [gcDistance, tdDistance, covDistance]
                for i, distance in enumerate(distances):
                    if distance < minDistances[i]:
                        minDistances[i] = distance
                        closestCluster[i] = clusterId
                        
            # check if unbinned contig is closest to a single core is all three spaces
            if closestCluster[0] != None:
                if closestCluster[0] == closestCluster[1] and closestCluster[0] == closestCluster[2]:
                    numAssigned += 1
                    unbinnedContig.clusterNum = cores[closestCluster[0]].clusterNum
                    assignedBP += unbinnedContig.seqLen
                else:
                    numFailingDistanceTest += 1
                    unassignedBP += unbinnedContig.seqLen
            else:
                numInvalidWithAllCores += 1
                unassignedBP += unbinnedContig.seqLen
                
        if self.logger.getEffectiveLevel() <= logging.INFO:
            statusStr = '    Processed %d of %d (%.2f%%) unbinned contigs.' % (processedSeqs, len(unbinnedContigs), float(processedSeqs)*100/len(unbinnedContigs))
            sys.stdout.write('%s\n' % statusStr)
            sys.stdout.flush()
                
        numUnassigned = len(unbinnedContigs)-numAssigned
        self.logger.info('    Assigned %d of %d (%.2f%%) unbinned contigs to a core bin.' % (numAssigned, len(unbinnedContigs), numAssigned*100.0/len(unbinnedContigs)))
        self.logger.info('      - %d of %d (%.2f%%) unassigned contigs were invalid with all cores' % (numInvalidWithAllCores, numUnassigned, numInvalidWithAllCores*100.0/numUnassigned))
        self.logger.info('      - %d of %d (%.2f%%) unassigned contigs failed distance test' % (numFailingDistanceTest, numUnassigned, numFailingDistanceTest*100.0/numUnassigned))
        self.logger.info('    Assigned %.2f of %.2f Mbps (%.2f%%) of the unbinned bases.' % (float(assignedBP)/1e6, float(assignedBP+unassignedBP)/1e6, float(assignedBP)*100/(assignedBP+unassignedBP)))
        
    def __resolveConflictingContigs(self, unbinnedContigs, binnedContigs):
        # map contigs (or the partitions comprising a contig) 
        # to their originating scaffolds
        scaffoldsToContigs = defaultdict(list)
        for contig in itertools.chain(unbinnedContigs, binnedContigs):
            match = re.search('_c[0-9]', contig.seqId)
            if match:
                scaffoldId = contig.seqId[0:match.start()]
                scaffoldsToContigs[scaffoldId].append(contig)
                
        rc = ResolveConflicts()
        rc.run(scaffoldsToContigs, DBSCAN.NOISE)
        
        self.logger.info('')
        self.logger.info('  Resolving scaffolds in conflict: ')
        self.logger.info('    - %d of %d (%.2f%%) partitioned scaffolds completely within a single bin' % (rc.numProperBin, rc.numPartitionedSeqs, rc.numProperBin*100.0/rc.numPartitionedSeqs))
        self.logger.info('    - %d of %d (%.2f%%) partitioned scaffolds were completely unbinned' % (rc.numNoise, rc.numPartitionedSeqs, rc.numNoise*100.0/rc.numPartitionedSeqs))
        
        self.logger.info('    - %d of %d (%.2f%%) partitioned scaffolds were partially unbinned' % (rc.numPartialNoise, rc.numPartitionedSeqs, rc.numPartialNoise*100.0/rc.numPartitionedSeqs))
        self.logger.info('      - %d (%.2f%%) resolved by reassignment to a single bin' % (rc.numPartialReassignment, rc.numPartialReassignment*100.0/max(rc.numPartialNoise, 1)))
        self.logger.info('      - %d (%.2f%%) resolved by marking all partitions as noise' % (rc.numPartialToNoise, rc.numPartialToNoise*100.0/max(rc.numPartialNoise, 1)))
        
        self.logger.info('    - %d of %d (%.2f%%) partitioned scaffolds were assigned to multiple bins' % (rc.numConflicts, rc.numPartitionedSeqs, rc.numConflicts*100.0/rc.numPartitionedSeqs))
        self.logger.info('      - %d (%.2f%%) resolved by reassignment to a single bin' % (rc.numConflictsReassignment, rc.numConflictsReassignment*100.0/max(rc.numConflicts, 1)))
        self.logger.info('      - %d (%.2f%%) resolved by marking all partitions as noise' % (rc.numConflictsToNoise, rc.numConflictsToNoise*100.0/max(rc.numConflicts, 1)))

    def __reportClusters(self, refinedBinFile, unbinnedContigs, binnedContigs, argStr):
        # determine contigs in each cluster
        clusterIdToContigs = defaultdict(set)
        for contig in unbinnedContigs:
            clusterIdToContigs[contig.clusterNum].add(contig)
            
        for contig in binnedContigs:
            clusterIdToContigs[contig.clusterNum].add(contig)
        
        # determine size of each bin in base pairs
        clusterSizeBP = defaultdict(int)
        totalClusteredBP = 0
        totalUnclusteredBP = 0
        totalClusteredContigs = 0
        for clusterId in sorted(clusterIdToContigs):
            clusterSizeBP = 0
            for contig in clusterIdToContigs[clusterId]:
                clusterSizeBP += contig.seqLen
                if clusterId == DBSCAN.NOISE:
                    totalUnclusteredBP += contig.seqLen
                else:
                    totalClusteredBP += contig.seqLen
            
            if clusterId != DBSCAN.NOISE:
                self.logger.info('    Contigs in bin %d: %d (%.2f Mbp)' % (clusterId, len(clusterIdToContigs[clusterId]), float(clusterSizeBP)/1e6))
                totalClusteredContigs += len(clusterIdToContigs[clusterId])
     
        self.logger.info('')
        self.logger.info('    Total binned contigs: %d (%.2f Mbp)' % (totalClusteredContigs, float(totalClusteredBP)/1e6))
        self.logger.info('    Total unbinned contigs: %d (%.2f Mbp)' % (len(clusterIdToContigs[DBSCAN.NOISE]), float(totalUnclusteredBP)/1e6))
        
        # save clustering to file
        fout = open(refinedBinFile, 'w')
        fout.write('# ' + argStr + '\n')
        fout.write('Contig Id\tCluster Id\tLength\tGC\tCoverage\n')
        for contigId, contigs in clusterIdToContigs.iteritems():
            for contig in contigs:
                fout.write(contig.seqId + '\t%d\t%d\t%.3f\t%.3f\n' % (contigId, contig.seqLen, contig.GC, contig.coverage))
                
        for contig in clusterIdToContigs[DBSCAN.NOISE]:
            fout.write(contig.seqId + '\t%s\t%d\t%.3f\t%.3f\n' % ('unbinned', contig.seqLen, contig.GC, contig.coverage))
        
        fout.close()
        
        self.logger.info('')
        self.logger.info('    Binning information written to: ' + refinedBinFile)
               
    def run(self, preprocessDir, binningFile, minSeqLen, gcDistPer, tdDistPer, covDistPer, refinedBinFile, argStr):
        # verify inputs
        contigStatsFile = os.path.join(preprocessDir, 'contigs.seq_stats.tsv')
        checkFileExists(contigStatsFile)
        
        contigTetraFile = os.path.join(preprocessDir, 'contigs.tetra.tsv')
        checkFileExists(contigTetraFile)
        
        coverageDistFile = os.path.join(preprocessDir, 'coverage_dist.txt')
        checkFileExists(coverageDistFile)
        
        # read tetranucleotide signatures
        self.logger.info('  Reading contig statistics.')
        contigStats = readSeqStats(contigStatsFile)
 
        # read tetranucleotide signatures
        self.logger.info('  Reading contig tetranucleotide signatures.')
        genomicSig = GenomicSignatures(4, 1)
        tetraSigs = genomicSig.read(contigTetraFile)
        
        # read bin assignments
        self.logger.info('  Reading core bin assignments.')
        contigIdToClusterId = {}
        with open(binningFile) as f:
            f.readline() # parameter line 
            f.readline() # header line
            for line in f:
                lineSplit = line.split('\t')
                clusterId = DBSCAN.NOISE if lineSplit[1] == 'unbinned' else int(lineSplit[1])
                contigIdToClusterId[lineSplit[0]] = clusterId
    
        # calculate statistics of core bins and unbinned contigs
        self.logger.info('')
        self.logger.info('  Calculating statistics of core bins.')

        cores = {}
        unbinnedContigs = set()
        binnedContigs = set()
        for contigId, contigStat in contigStats.iteritems():  
            clusterId = contigIdToClusterId.get(contigId, DBSCAN.NOISE)
            
            contigLen, contigGC, contigCov = contigStat
            contigTetraSig = tetraSigs[contigId]
            
            seq = Sequence(0, contigId, contigLen, contigGC, contigCov)
            seq.clusterNum = clusterId
            seq.tetraSig = contigTetraSig
            
            if clusterId == DBSCAN.NOISE: 
                if contigLen > minSeqLen:
                    unbinnedContigs.add(seq)
            else:
                binnedContigs.add(seq)
 
                # build statistics for core bins
                if clusterId not in cores:
                    cores[clusterId] = Sequence(clusterId, 'core_' + str(clusterId), 0, 0, 0)
                    cores[clusterId].clusterNum = clusterId
                    cores[clusterId].tetraSig = [0]*genomicSig.numKmers()
                    
                core = cores[clusterId]

                core.seqLen += contigLen
                weight = float(contigLen)/core.seqLen
                core.GC = contigGC*weight + core.GC*(1.0-weight)
                core.coverage = contigCov*weight + core.coverage*(1.0-weight)
                for i in xrange(0, len(core.tetraSig)):
                    core.tetraSig[i] = contigTetraSig[i]*weight + core.tetraSig[i]*(1.0-weight)
        
        self.logger.info('    Identified %d unbinned contigs >= %d bps.' % (len(unbinnedContigs), minSeqLen))
        self.logger.info('    Identified %d binned contigs from %d core bins.' % (len(binnedContigs), len(cores)))

        # read GC, TD and coverage distributions
        self.logger.info('')
        self.logger.info('  Reading GC, TD, and coverage distributions.')
        gcDist = readDistribution('gc_dist', gcDistPer)
        tdDist = readDistribution('td_dist', tdDistPer)
        
        with open(coverageDistFile) as f:
            s = f.read()
            covDist = ast.literal_eval(s)
            
        # refine bins
        self.logger.info('')
        self.logger.info('  Refining bins:')
        self.__refineBins(unbinnedContigs, cores, gcDist, tdDist, covDist, covDistPer, genomicSig)
                
        # resolve conflicts between contigs originating on the same scaffold
        # and use the conflicts as a measure of the binning quality
        self.__resolveConflictingContigs(unbinnedContigs, binnedContigs)
           
        # report and save results of binning
        self.logger.info('')
        self.logger.info('  Refined bins:')
        self.__reportClusters(refinedBinFile, unbinnedContigs, binnedContigs, argStr)
        
