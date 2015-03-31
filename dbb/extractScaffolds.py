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
from seqUtils import readFasta, readSeqStats

from greedy import Greedy, readSeqStatsForBins


class ExtractScaffolds(object):
    def __init__(self):
        self.logger = logging.getLogger()

    def __saveBins(self, preprocessDir, scaffoldFile, binningFile, outputDir):
        # verify inputs
        scaffoldStatsFile = os.path.join(preprocessDir, 'scaffolds.seq_stats.tsv')
        checkFileExists(scaffoldStatsFile)

        # read scaffold stats
        self.logger.info('    Reading scaffold statistics.')
        scaffoldStats = readSeqStats(scaffoldStatsFile)

        # read bin assignment of each contig partition
        self.logger.info('    Reading bin assignments.')
        seqStatsForBins = readSeqStatsForBins(binningFile)

        # get bin for each scaffold
        binIdToScaffoldIds = defaultdict(set)
        for binId, contigStats in seqStatsForBins.iteritems():
            if binId == Greedy.UNBINNED:
                continue

            for contigStats in contigStats:
                contigId, _, _, _ = contigStats
                match = re.search('_c[0-9]*$', contigId)
                scaffoldId = contigId[0:match.start()] if match else contigId

                binIdToScaffoldIds[binId].add(scaffoldId)

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
        for binId, scaffoldIds in binIdToScaffoldIds.iteritems():
            binnedBPs = 0

            fout = open(os.path.join(outputDir, 'bin_' + str(binId) + '.fna'), 'w')
            for scaffoldId in scaffoldIds:
                stats = scaffoldStats[scaffoldId]
                fout.write('>' + scaffoldId + ' [length=%d,GC=%.2f,coverage=%.1f]\n' % (stats[0], stats[1], stats[2]))
                fout.write(scaffolds[scaffoldId] + '\n')
                binnedBPs += len(scaffolds[scaffoldId])
            fout.close()

            self.logger.info('    Bin %d: %d scaffolds (%.2f Mbp)' % (binId, len(scaffoldIds), float(binnedBPs) / 1e6))

            totalBinnedBPs += binnedBPs
            totalBinnedScaffolds += len(scaffoldIds)

        self.logger.info('')
        self.logger.info('    Total binned scaffolds: %d (%.2f Mbp)' % (totalBinnedScaffolds, float(totalBinnedBPs) / 1e6))
        self.logger.info('    Total unbinned scaffolds: %d (%.2f Mbp)' % (len(scaffolds) - totalBinnedScaffolds, float(totalBPs - totalBinnedBPs) / 1e6))

    def run(self, preprocessDir, scaffoldFile, binningFile, binDir):
        # verify inputs
        checkFileExists(scaffoldFile)
        checkFileExists(binningFile)

        if not isDirEmpty(binDir):
            self.logger.warning('  [Warning] Specified output path for bins must be empty: ' + binDir + '\n')
            sys.exit()

        # report and save results of binning
        self.__saveBins(preprocessDir, scaffoldFile, binningFile, binDir)
