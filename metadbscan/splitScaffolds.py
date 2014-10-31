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
        
    def splits(self, scaffoldId, scaffold, minN):
        contigs = {}
        
        contigCount = 0
        nCount = 0
        startIndex = 0
        for curIndex, ch in enumerate(scaffold):
            if ch == 'N' or ch == 'n':
                if nCount == 0:
                    nStart = curIndex
                nCount += 1
            else:
                if nCount > minN:
                    contigs[scaffoldId + '_c' + str(contigCount)] = [startIndex, nStart]
                    
                    contigCount += 1
                    startIndex = curIndex

                nCount = 0

        if startIndex != len(scaffold)-1:
            if nCount == 0:
                contigs[scaffoldId + '_c' + str(contigCount)] = [startIndex, len(scaffold)]
            else:
                contigs[scaffoldId + '_c' + str(contigCount)] = [startIndex, nStart]
                
        return contigs

    def splitFile(self, seqFile, minN, outputFile):
        # read scaffolds
        self.logger.info('  Reading scaffolds.')
        scaffolds = readFasta(seqFile)
        self.logger.info('    Read %d scaffolds.' % len(scaffolds))

        # extract contigs from scaffolds
        self.logger.info('  Extracting contigs.')

        fout = open(outputFile, 'w')
        processedSeqs = 0
        numContigs = 0
        for scaffoldId, scaffold in scaffolds.iteritems():
            contigs = self.splits(scaffoldId, scaffold, minN)
            for contigId, split in contigs.iteritems():
                fout.write('>' + contigId + '\n')
                fout.write(scaffolds[split[0]:split[1]])
                numContigs += 1

            if self.logger.getEffectiveLevel() <= logging.INFO:
                processedSeqs += 1
                statusStr = '     Finished processing %d of %d (%.2f%%) scaffolds.' % (processedSeqs, len(scaffolds), float(processedSeqs)*100/len(scaffolds))
                sys.stdout.write('%s\r' % statusStr)
                sys.stdout.flush()

        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stdout.write('\n')
            
        fout.close()
        
        self.logger.info('  Extracted %d contigs.' % numContigs)
