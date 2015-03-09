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

import logging
from collections import Counter, defaultdict

class ResolveConflicts(object):
    def __init__(self):
        self.logger = logging.getLogger()
        
        self.numPartitionedSeqs = 0
        
        self.numConflicts = 0
        self.numConflictsReassignment = 0
        self.numConflictsToNoise = 0 
        
        self.numPartialNoise = 0
        self.numPartialReassignment = 0
        self.numPartialToNoise = 0
        
        self.numNoise = 0
        self.numProperBin = 0
        
    def merge(self, seqsToPartitions, gcDist, tdDist, covDist, mergeDistPer, unbinnedId):
        '''Identify partitions spanning two or more bins and determine if they should be merged.'''
        pass
           
    def resolve(self, seqsToPartitions, unbinnedId):
        ''' Identify and resolve conflicting bin assignments between sequence partitions (e.g., contigs from a scaffold)''' 
        self.numPartitionedSeqs = 0
        
        self.numConflicts = 0
        self.numConflictsReassignment = 0
        self.numConflictsToNoise = 0 
        
        self.numPartialNoise = 0
        self.numPartialReassignment = 0
        self.numPartialToNoise = 0
        
        self.numNoise = 0
        self.numProperBin = 0
        #mergeVotes = defaultdict(int)
        for _, partitions in seqsToPartitions.iteritems():
            if len(partitions) == 1:
                continue
            
            self.numPartitionedSeqs += 1
            binVotes = []
            noiseVotes = 0
            for contig in partitions:
                if contig.binId != unbinnedId:
                    binVotes.append(contig.binId)
                else:
                    noiseVotes += 1
                    
            counter = Counter(binVotes) 
            if len(counter) == 1 and noiseVotes == 0:
                # all partitions are in a single bin
                self.numProperBin += 1
            elif len(binVotes) == 0:
                # all partitions assigned as noise
                self.numNoise += 1
            elif len(counter) == 1 and noiseVotes > 0:
                # partitions are partially assigned to a single bin
                self.numPartialNoise += 1
                
                # determine most common assignment for partitions
                binId, partitionsInTopBin = counter.most_common(1)[0]

                # resolve conflict by assigning partitions to the bin if
                # AT LEAST 50% of the partitions are binned
                if partitionsInTopBin >= 0.5*len(partitions):
                    self.numPartialReassignment += 1
                    newBinId = binId
                else:
                    self.numPartialToNoise += 1
                    newBinId = unbinnedId
                    
                for p in partitions:
                    p.binId = newBinId
            else:
                # partitions span multiple bins
                self.numConflicts += 1
                
                # determine most common assignment for partitions
                topBin, partitionsInTopBin = counter.most_common(1)[0]
                
                # resolve conflict by assigning partitions to the majority bin if
                # GREATER THAN 50% of the partitions are within the majority bin
                if partitionsInTopBin > 0.5*len(partitions):
                    # move all partitions into the most likely bin
                    self.numConflictsReassignment += 1
                    newBinId = topBin
                else:
                    # mark partitions as noise as their correct binning is unclear
                    self.numConflictsToNoise += 1
                    newBinId = unbinnedId
                    
                    #secondBin, _ = counter.most_common(2)[1]
                    #mergeString = str(min(topBin, secondBin)) + '-' + str(max(topBin, secondBin))
                    #mergeVotes[mergeString] += 1
                    
                for p in partitions:
                    p.binId = newBinId
                    
        #print mergeVotes
