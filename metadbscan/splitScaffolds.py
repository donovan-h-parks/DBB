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

import sys
import logging

from seqUtils import readFasta

class SplitScaffolds(object):
    def __init__(self):
        self.logger = logging.getLogger()

    def run(self, seqFile, minN, outputFile):
        # read scaffolds
        self.logger.info('  Reading scaffolds.')
        seqs = readFasta(seqFile)
        self.logger.info('    Read %d scaffolds.' % len(seqs))

        # extract contigs from scaffolds
        self.logger.info('  Extracting contigs.')

        fout = open(outputFile, 'w')
        processedSeqs = 0
        numContigs = 0
        for seqId, seq in seqs.iteritems():
            contigCount = 0
            nCount = 0
            startIndex = 0
            for curIndex, ch in enumerate(seq):
                if ch == 'N':
                    if nCount == 0:
                        nStart = curIndex
                    nCount += 1
                else:
                    if nCount > minN:
                        contig = seq[startIndex:nStart]
                        
                        fout.write('>' + seqId + '_c' + str(contigCount) + '\n')
                        fout.write(contig + '\n')
                        
                        contigCount += 1
                        numContigs += 1
                        startIndex = curIndex

                    nCount = 0

            if startIndex != len(seq)-1:
                fout.write('>' + seqId + '_c' + str(contigCount) + '\n')
                if nCount == 0:
                    contig = seq[startIndex:]
                else:
                    contig = seq[startIndex:nStart]

                fout.write(contig + '\n')
                numContigs += 1

            if self.logger.getEffectiveLevel() <= logging.INFO:
                processedSeqs += 1
                statusStr = '     Finished processing %d of %d (%.2f%%) scaffolds.' % (processedSeqs, len(seqs), float(processedSeqs)*100/len(seqs))
                sys.stdout.write('%s\r' % statusStr)
                sys.stdout.flush()

        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stdout.write('\n')
            
        fout.close()
        
        self.logger.info('  Extracted %d contigs.' % numContigs)
