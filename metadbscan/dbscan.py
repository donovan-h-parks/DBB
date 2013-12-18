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

__prog_desc__ = 'genome aware DBSCAN'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import ast
import logging
import re
from collections import defaultdict, Counter

import numpy as np

from common import checkFileExists
from seqUtils import readSeqStats
from distributions import readDistribution

from resolveConflicts import ResolveConflicts

from coverageNeighbours import CoverageNeighbours
from gcNeighbours import GcNeighbours
from tdNeighbours import TdNeighbours

def readSeqStatsForClusters(binningFile):
    seqStatsForClusters = defaultdict(list)
    with open(binningFile) as f:
        f.readline()
        for line in f:
            lineSplit = line.split('\t')
            clusterId = DBSCAN.NOISE if lineSplit[1] == 'unbinned' else int(lineSplit[1])
            seqStatsForClusters[clusterId].append([lineSplit[0], int(lineSplit[2]), float(lineSplit[3]), float(lineSplit[4])])
    
    return seqStatsForClusters

class Sequence(object):
    def __init__(self, seqIndex, seqId, seqLen, GC, coverage):
        self.bVisited = False
        self.clusterNum = DBSCAN.NOISE
        self.seqIndex = seqIndex
        
        self.seqId = seqId
        self.seqLen = seqLen
        self.GC = GC
        self.coverage = coverage
        
        self.covNeighbours = None
        self.gcNeighbours = None
        self.tdNeighbours = None
        
class DBSCAN(object):
    NOISE = -1
    
    def __init__(self, bDebug=False):
        self.logger = logging.getLogger()
        
        self.numClusters = 0  
        self.sortedSeqs = None
    
    def __neighbourIndices(self, seq):
        covNeighbours = np.where(seq.covNeighbours == 1)[0]
        gcNeighbours = np.where(seq.gcNeighbours == 1)[0]
        tdNeighbours = np.where(seq.tdNeighbours == 1)[0]
        
        intersectCovGC = np.intersect1d(covNeighbours, gcNeighbours, assume_unique=True)
        return np.intersect1d(tdNeighbours, intersectCovGC, assume_unique=True)
        
    def __neighbours(self, neighbourIndices):
        neighbourSeqs = set()
        for i in neighbourIndices:
            neighbourSeqs.add(self.sortedSeqs[i])
        return neighbourSeqs
    
    def __dbscan(self, seqs, minPts, minCoreLen):
        for seq in seqs:
            if seq.bVisited:
                continue

            seq.bVisited = True
            neighbourIndices = self.__neighbourIndices(seq)

            if len(neighbourIndices) >= minPts and seq.seqLen >= minCoreLen:
                # this point is core so expand it's neighbourhood
                neighbours = self.__neighbours(neighbourIndices)
                self.__expandCluster(seq, neighbours, self.numClusters, minPts, minCoreLen)
                self.numClusters += 1
         
    def __expandCluster(self, seq, neighbours, clusterNum, minPts, minCoreLen):
        seq.clusterNum = clusterNum   
            
        while neighbours:
            curSeq = neighbours.pop()
            
            if not curSeq.bVisited:
                curSeq.bVisited = True
                neighbourIndices = self.__neighbourIndices(curSeq)
                if len(neighbourIndices) >= minPts and curSeq.seqLen >= minCoreLen:
                    # this is a core point so add it's neighbours to the neighbourhood list
                    neighbours.update(self.__neighbours(neighbourIndices))
                
            if curSeq.clusterNum == DBSCAN.NOISE:
                curSeq.clusterNum = clusterNum
                
    def __filterBins(self, minBinSize):
        # get partitions in each cluster
        clusters = defaultdict(list)
        for seq in self.sortedSeqs:
            clusters[seq.clusterNum].append(seq)
        
        # determine size of each cluster in base pairs
        clusterSizeBP = defaultdict(int)
        for i in xrange(DBSCAN.NOISE, self.numClusters):
            for seq in clusters[i]:
                clusterSizeBP[i] +=  seq.seqLen
                
        # filter out bins of insufficient size and
        # adjust clustering numbers so they are still consecutive
        adjClusterNum = 0
        for i in xrange(0, self.numClusters):
            if clusterSizeBP[i] < minBinSize:
                for seq in clusters[i]:
                    seq.clusterNum = DBSCAN.NOISE
            else:
                for seq in clusters[i]:
                    seq.clusterNum = adjClusterNum
                    
                adjClusterNum += 1
                
        self.numClusters = adjClusterNum
              
    def __resolveConflictingPartitions(self):
        # map partitions to their originating contigs
        contigsToPartitions = defaultdict(list)
        for seq in self.sortedSeqs:
            match = re.search('_p[0-9]', seq.seqId)
            if match:
                contigId = seq.seqId[0:match.start()]
                contigsToPartitions[contigId].append(seq)
 
        rc = ResolveConflicts()
        rc.run(contigsToPartitions, DBSCAN.NOISE)
                      
        self.logger.info('')
        self.logger.info('  Resolving contigs in conflict: ')
        self.logger.info('    - %d of %d (%.2f%%) partitioned contigs completely within a single bin' % (rc.numProperBin, rc.numPartitionedSeqs, rc.numProperBin*100.0/rc.numPartitionedSeqs))
        self.logger.info('    - %d of %d (%.2f%%) partitioned contigs were completely unbinned' % (rc.numNoise, rc.numPartitionedSeqs, rc.numNoise*100.0/rc.numPartitionedSeqs))
        
        self.logger.info('    - %d of %d (%.2f%%) partitioned contigs were at least partially unbinned' % (rc.numPartialNoise, rc.numPartitionedSeqs, rc.numPartialNoise*100.0/rc.numPartitionedSeqs))
        self.logger.info('      - %d (%.2f%%) resolved by reassignment to a single bin' % (rc.numPartialReassignment, rc.numPartialReassignment*100.0/max(rc.numPartialNoise, 1)))
        self.logger.info('      - %d (%.2f%%) resolved by marking all partitions as noise' % (rc.numPartialToNoise, rc.numPartialToNoise*100.0/max(rc.numPartialNoise, 1)))
         
        self.logger.info('    - %d of %d (%.2f%%) partitioned contigs were assigned to multiple bins' % (rc.numConflicts, rc.numPartitionedSeqs, rc.numConflicts*100.0/rc.numPartitionedSeqs))
        self.logger.info('      - %d (%.2f%%) resolved by reassignment to a single bin' % (rc.numConflictsReassignment, rc.numConflictsReassignment*100.0/max(rc.numConflicts, 1)))
        self.logger.info('      - %d (%.2f%%) resolved by marking all partitions as noise' % (rc.numConflictsToNoise, rc.numConflictsToNoise*100.0/max(rc.numConflicts, 1)))
           
    def __resolveConflictingContigs(self):
        # map contigs (or the partitions comprising a contig) 
        # to their originating scaffolds
        scaffoldsToSeqs = defaultdict(list)
        for seq in self.sortedSeqs:
            match = re.search('_c[0-9]', seq.seqId)
            if match:
                scaffoldId = seq.seqId[0:match.start()]
                scaffoldsToSeqs[scaffoldId].append(seq)
 
        rc = ResolveConflicts()
        rc.run(scaffoldsToSeqs, DBSCAN.NOISE)
                      
        self.logger.info('')
        self.logger.info('  Resolving scaffolds in conflict: ')
        self.logger.info('    - %d of %d (%.2f%%) partitioned scaffolds completely within a single bin' % (rc.numProperBin, rc.numPartitionedSeqs, rc.numProperBin*100.0/rc.numPartitionedSeqs))
        self.logger.info('    - %d of %d (%.2f%%) partitioned scaffolds were completely unbinned' % (rc.numNoise, rc.numPartitionedSeqs, rc.numNoise*100.0/rc.numPartitionedSeqs))
        
        self.logger.info('    - %d of %d (%.2f%%) partitioned scaffolds were at least partially unbinned' % (rc.numPartialNoise, rc.numPartitionedSeqs, rc.numPartialNoise*100.0/rc.numPartitionedSeqs))
        self.logger.info('      - %d (%.2f%%) resolved by reassignment to a single bin' % (rc.numPartialReassignment, rc.numPartialReassignment*100.0/max(rc.numPartialNoise, 1)))
        self.logger.info('      - %d (%.2f%%) resolved by marking all partitions as noise' % (rc.numPartialToNoise, rc.numPartialToNoise*100.0/max(rc.numPartialNoise, 1)))
        
        self.logger.info('    - %d of %d (%.2f%%) partitioned scaffolds were assigned to multiple bins' % (rc.numConflicts, rc.numPartitionedSeqs, rc.numConflicts*100.0/rc.numPartitionedSeqs))
        self.logger.info('      - %d (%.2f%%) resolved by reassignment to a single bin' % (rc.numConflictsReassignment, rc.numConflictsReassignment*100.0/max(rc.numConflicts, 1)))
        self.logger.info('      - %d (%.2f%%) resolved by marking all partitions as noise' % (rc.numConflictsToNoise, rc.numConflictsToNoise*100.0/max(rc.numConflicts, 1)))
          
    def __reportClusters(self, binningFile, contigStats):
        # get partitions in each cluster
        clusters = defaultdict(list)
        for seq in self.sortedSeqs:
            clusters[seq.clusterNum].append(seq)
        
        # determine contigs in each cluster
        totalClusteredBP = 0
        totalUnclusteredBP = 0
        clusterIdToContigIds = defaultdict(set)
        for i in xrange(DBSCAN.NOISE, self.numClusters):
            clusterSizeBP = 0
            for seq in clusters[i]:
                clusterSizeBP += seq.seqLen
                if i == DBSCAN.NOISE:
                    totalUnclusteredBP += seq.seqLen
                else:
                    totalClusteredBP += seq.seqLen
                    
                match = re.search('_p[0-9]', seq.seqId)
                contigId = seq.seqId[0:match.start()] if match else seq.seqId
                clusterIdToContigIds[i].add(contigId)
                
            if i != DBSCAN.NOISE:
                self.logger.info('    Contigs in bin %d: %d (%.2f Mbp)' % (i, len(clusterIdToContigIds[i]), float(clusterSizeBP)/1e6))
                
        totalClusteredContigs = sum([len(contigIds) for clusterId, contigIds in clusterIdToContigIds.iteritems() if clusterId != DBSCAN.NOISE])
           
        self.logger.info('')
        self.logger.info('    Total binned contigs: %d (%.2f Mbp)' % (totalClusteredContigs, float(totalClusteredBP)/1e6))
        self.logger.info('    Total unbinned contigs: %d (%.2f Mbp)' % (len(clusterIdToContigIds[DBSCAN.NOISE]), float(totalUnclusteredBP)/1e6))
        
        # save clustering to file
        fout = open(binningFile, 'w')
        fout.write('Contig Id\tCluster Id\tLength\tGC\tCoverage\n')
        for clusterId, contigIds in clusterIdToContigIds.iteritems():
            clusterName = str(clusterId) if clusterId != DBSCAN.NOISE else 'unbinned'
                
            for contigId in contigIds:
                contigLen, contigGC, contigCov = contigStats[contigId]
                fout.write(contigId + '\t%s\t%d\t%.3f\t%.3f\n' % (clusterName, contigLen, contigGC, contigCov))
        
        fout.close()
        
        self.logger.info('')
        self.logger.info('    Binning information written to: ' + binningFile)
               
    def run(self, preprocessDir, minSeqLen, minBinSize, gcDistPer, tdDistPer, covDistPer, minPts, minCoreLen, numThreads, binningFile):
        # verify inputs
        seqStatsFile = os.path.join(preprocessDir, 'partitions.seq_stats.tsv')
        checkFileExists(seqStatsFile)
        
        tetraFile = os.path.join(preprocessDir, 'partitions.tetra.tsv')
        checkFileExists(tetraFile)
        
        coverageDistFile = os.path.join(preprocessDir, 'coverage_dist.txt')
        checkFileExists(coverageDistFile)
        
        contigStatsFile = os.path.join(preprocessDir, 'contigs.seq_stats.tsv')
        checkFileExists(contigStatsFile)

        # create sequence
        self.logger.info('  Reading partitioned contig statistics.')
        self.sortedSeqs = []
        self.seqIdsOfInterest = set()
        numPotentialCores = 0
        with open(seqStatsFile) as f:
            f.readline()
            
            seqIndex = 0
            for line in f:
                lineSplit = line.split('\t')
                
                seqId = lineSplit[0]
                seqLen = int(lineSplit[1])
                GC = float(lineSplit[2])
                coverage = float(lineSplit[3])
                
                if seqLen > minSeqLen:
                    self.sortedSeqs.append(Sequence(seqIndex, seqId, seqLen, GC, coverage))
                    self.seqIdsOfInterest.add(seqId)
                    seqIndex += 1
                    
                    if seqLen > minCoreLen:
                        numPotentialCores += 1
                    
        self.logger.info('    Partitioned contig > %d bps: %d' % (minSeqLen, len(self.sortedSeqs)))
        self.logger.info('    Potential cores (partitioned contig > %d bps): %d' % (minCoreLen, numPotentialCores))
        
        # sort sequences in ascending order of length
        self.logger.info('')
        self.logger.info('  Sorting partitioned contig in ascending order of length.')
        self.sortedSeqs.sort(key=lambda x: x.seqLen, reverse=True)
                
        # read GC, TD and coverage distributions
        self.logger.info('')
        self.logger.info('  Reading GC, TD and coverage distributions.')
        gcDist = readDistribution('gc_dist', gcDistPer)
        tdDist = readDistribution('td_dist', tdDistPer)

        with open(coverageDistFile) as f:
            s = f.read()
            covDist = ast.literal_eval(s)
                
        # calculate coverage, GC, and tetranucleotide neighbours
        self.logger.info('')
        self.logger.info('  Calculating coverage neighbours:')
        coverageNeighbours = CoverageNeighbours()
        coverageNeighbours.run(self.sortedSeqs, covDist, covDistPer, numThreads)
        
        self.logger.info('  Calculating GC neighbours:')
        gcNeighbours = GcNeighbours()
        gcNeighbours.run(self.sortedSeqs, gcDist, numThreads)
        
        self.logger.info('  Calculating TD neighbours:')
        tdNeighbours = TdNeighbours()
        tdNeighbours.run(self.seqIdsOfInterest, self.sortedSeqs, tetraFile, tdDist, numThreads)

        # perform clustering via genome-aware DBSCAN
        self.logger.info('')
        self.logger.info('  Running DBSCAN.')
        self.__dbscan(self.sortedSeqs, minPts, minCoreLen)
        self.logger.info('    Found %d clusters.' % self.numClusters)
        
        # resolve conflicts between partitions originating on the same contig
        # and use the conflicts as a measure of the binning quality
        self.__resolveConflictingPartitions()
        
        # resolve conflicts between contigs originating on the same scaffold
        # and use the conflicts as a measure of the binning quality
        self.__resolveConflictingContigs()
        
        # filter bins sufficiently small that they may belong to another bin
        # and from which there is limited data too calculate robust statistics
        self.logger.info('')
        self.logger.info('  Filtering bins with < %d bps.' % minBinSize)
        self.__filterBins(minBinSize)
        self.logger.info('    Retained %d bins.' % self.numClusters)
        
        # report and save results of binning
        self.logger.info('')
        self.logger.info('  Core bins:')
        contigStats = readSeqStats(contigStatsFile)
        self.__reportClusters(binningFile, contigStats)
        
