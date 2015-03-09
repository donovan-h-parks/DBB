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

import re
import logging

from greedy import Greedy, readSeqStatsForBins


class Expand(object):
    def __init__(self, bDebug=False):
        self.logger = logging.getLogger()

    def run(self, unbinnedStatsFile, binningFile, scoreThreshold, outputFile):
        # read bin assignments
        self.logger.info('  Reading current bin statistics.')
        seqStatsForBins = readSeqStatsForBins(binningFile)

        # get unbinned scaffolds to bin
        self.logger.info('  Determining unbinned scaffolds to bin.')
        scaffoldsToBin = {}
        with open(unbinnedStatsFile) as f:
            f.readline()  # skip header

            for line in f:
                lineSplit = line.split('\t')
                scaffoldId = lineSplit[0]
                binId = lineSplit[4]
                score = float(lineSplit[14])

                if score >= scoreThreshold:
                    scaffoldsToBin[scaffoldId] = binId

        # create new binning file
        self.logger.info('')
        self.logger.info('  Writing new bin file.')
        fout = open(outputFile, 'w')
        fout.write('Contig Id\tBin Id\tLength\tGC\tCoverage\n')

        numNewlyBinned = 0
        newlyBinnedBPs = 0
        for binId, scaffoldStats in seqStatsForBins.iteritems():
            for stats in scaffoldStats:
                contigId, length, GC, coverage = stats

                scaffoldId = contigId
                match = re.search('_c[0-9]*$', contigId)
                if match:
                    scaffoldId = contigId[0:match.start()]

                newBinId = binId
                if scaffoldId in scaffoldsToBin:
                    if binId == Greedy.UNBINNED:
                        numNewlyBinned += 1
                        newlyBinnedBPs += length

                    newBinId = scaffoldsToBin[scaffoldId]

                if newBinId != Greedy.UNBINNED:
                    fout.write('%s\t%s\t%d\t%.3f\t%.3f\n' % (contigId, newBinId, length, GC, coverage))
                else:
                    fout.write('%s\t%s\t%d\t%.3f\t%.3f\n' % (contigId, 'unbinned', length, GC, coverage))

        fout.close()

        self.logger.info('    Newly binned scaffolds: %d' % numNewlyBinned)
        self.logger.info('    Newly binned base pairs (Mbps): %.2f' % (float(newlyBinnedBPs) / 1e6))
