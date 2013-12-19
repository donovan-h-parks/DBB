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

import sys
import logging
import multiprocessing as mp

from distributions import findNearest, inRange

import numpy as np

class GcNeighbours(object):
    def __init__(self):
        self.logger = logging.getLogger()

    def __workerThread(self, seqs, gcDist, gcDistPer, queueIn, queueOut):
        """Process each data item in parallel."""
        
        sampleMeanGC = gcDist.keys()[0]
        sampleSeqLen = gcDist[sampleMeanGC].keys()[0]
        d = gcDist[sampleMeanGC][sampleSeqLen]
        lowerBoundKey = findNearest(d.keys(), (100 - gcDistPer)/2.0)
        upperBoundKey = findNearest(d.keys(), (100 + gcDistPer)/2.0)
            
        while True:
            index, seqI = queueIn.get(block=True, timeout=None)
            if index == None:
                break

            row = np.zeros(len(seqs))
            for j, seqJ in enumerate(seqs):
                if seqI.seqLen > seqJ.seqLen:
                    meanGC = seqI.GC
                    gcDiff = seqJ.GC - meanGC
                    seqLen = seqJ.seqLen
                else:
                    meanGC = seqJ.GC
                    gcDiff = seqI.GC - meanGC
                    seqLen = seqI.seqLen

                closestMeanGC = findNearest(gcDist.keys(), meanGC)
                closestSeqLen = findNearest(gcDist[closestMeanGC].keys(), seqLen)

                d = gcDist[closestMeanGC][closestSeqLen]
                lowerBound = d[lowerBoundKey]
                upperBound = d[upperBoundKey]
                
                if inRange(gcDiff, lowerBound, upperBound):
                    row[j] = 1

            # allow results to be processed or written to file
            queueOut.put((index, row))
    
    def __writerThread(self, neighbours, numDataItems, writerQueue):
        """Store or write results of worker threads in a single thread."""
        processedItems = 0
        while True:
            index, gcNeighbours = writerQueue.get(block=True, timeout=None)
            if index == None:
                break
            
            if self.logger.getEffectiveLevel() <= logging.INFO:
                processedItems += 1
                statusStr = '    Finished processing %d of %d (%.2f%%) sequences.' % (processedItems, numDataItems, float(processedItems)*100/numDataItems)
                sys.stdout.write('%s\r' % statusStr)
                sys.stdout.flush()
                
            neighbours[index] = gcNeighbours
                
        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stdout.write('\n')
    
    def run(self, seqs, gcDist, gcDistPer, threads):  
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()
        
        # populate worker queue with data to process
        for i, seq in enumerate(seqs):
            workerQueue.put((i, seq))
        
        for _ in range(threads):
            workerQueue.put((None, None))
        
        gcNeighbours = mp.Manager().list([None]*len(seqs))

        workerProc = [mp.Process(target = self.__workerThread, args = (seqs, gcDist, gcDistPer, workerQueue, writerQueue)) for _ in range(threads)]
        writeProc = mp.Process(target = self.__writerThread, args = (gcNeighbours, len(seqs), writerQueue))
        
        writeProc.start()
        
        for p in workerProc:
            p.start()
        
        for p in workerProc:
            p.join()
        
        writerQueue.put((None, None))
        writeProc.join()
        
        for i, row in enumerate(gcNeighbours):
            seqs[i].gcNeighbours = row
