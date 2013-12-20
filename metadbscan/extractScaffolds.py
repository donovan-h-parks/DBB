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

from common import checkFileExists, isDirEmpty
from seqUtils import readFasta

from dbscan import DBSCAN, readSeqStatsForClusters
       
class ExtractScaffolds(object):
    def __init__(self):
        self.logger = logging.getLogger()
        
    def __saveClusters(self, scaffoldFile, binningFile, outputDir):
        # read bin assignment of each contig partition
        seqStatsForClusters = readSeqStatsForClusters(binningFile)
        
        # get cluster for each scaffold
        clusterIdToScaffoldIds = defaultdict(set)
        for clusterId, contigStats in seqStatsForClusters.iteritems():
            if clusterId == DBSCAN.NOISE:
                continue
            
            for contigStats in contigStats:
                contigId, _, _, _ = contigStats
                match = re.search('_c[0-9]', contigId)
                scaffoldId = contigId[0:match.start()] if match else contigId
                    
                clusterIdToScaffoldIds[clusterId].add(scaffoldId)

        # read contigs
        scaffolds = readFasta(scaffoldFile)
        
        totalBPs = 0
        for scaffoldId, scaffold in scaffolds.iteritems():
            totalBPs += len(scaffold)
        
        # write bins to directory
        self.logger.info('')
        self.logger.info('  Bin scaffold statistics:')
        totalBinnedBPs = 0
        totalBinnedScaffolds = 0
        for clusterId, scaffoldIds in clusterIdToScaffoldIds.iteritems():
            binnedBPs = 0
            
            fout = open(os.path.join(outputDir, 'bin_' + str(clusterId) + '.fna'), 'w')
            for scaffoldId in scaffoldIds:
                fout.write('>' + scaffoldId + '\n')
                fout.write(scaffolds[scaffoldId] + '\n')
                binnedBPs += len(scaffolds[scaffoldId])   
            fout.close()
            
            self.logger.info('    Bin %d: %d scaffolds (%.2f Mbp)' % (clusterId, len(scaffoldIds), float(binnedBPs)/1e6))
            
            totalBinnedBPs += binnedBPs
            totalBinnedScaffolds += len(scaffoldIds)
             
        self.logger.info('')
        self.logger.info('    Total binned scaffolds: %d (%.2f Mbp)' % (totalBinnedScaffolds, float(totalBinnedBPs)/1e6))
        self.logger.info('    Total unbinned scaffolds: %d (%.2f Mbp)' % (len(scaffolds) - totalBinnedScaffolds, float(totalBPs - totalBinnedBPs)/1e6))
                 
    def run(self, scaffoldFile, binningFile, binDir):
        # verify inputs
        checkFileExists(scaffoldFile)
        checkFileExists(binningFile)
              
        if not isDirEmpty(binDir):
            self.logger.warning('  [Warning] Specified output path for bins must be empty: ' + binDir + '\n')
            sys.exit()
        
        # report and save results of binning
        self.__saveClusters(scaffoldFile, binningFile, binDir)
