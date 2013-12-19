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

from distributions import findNearest

import numpy as np

class TdNeighbours(object):
    def __init__(self):
        self.logger = logging.getLogger()
    
    def __workerThread(self, seqs, tetra, tdDist, tdDistPer, queueIn, queueOut):
        """Process each data item in parallel."""
        
        distKey = findNearest(tdDist[tdDist.keys()[0]].keys(), tdDistPer)
  
        while True:
            index, seqI = queueIn.get(block=True, timeout=None)
            if index == None:
                break

            curTetraSig = tetra[seqI.seqId]
            
            row = np.zeros(len(seqs))
            for j, seqJ in enumerate(seqs):  
                tetraSig = tetra[seqJ.seqId]        
                dist = np.sum(np.abs(curTetraSig - tetraSig))
                
                minLen = min(seqI.seqLen, seqJ.seqLen)
                if dist < tdDist[findNearest(tdDist.keys(), minLen)][distKey]:
                    row[j] = 1
  
            # allow results to be processed or written to file
            queueOut.put((index, row))
    
    def __writerThread(self, neighbours, numDataItems, writerQueue):
        """Store or write results of worker threads in a single thread."""
        processedItems = 0
        while True:
            index, tdNeighbours = writerQueue.get(block=True, timeout=None)
            if index == None:
                break
            
            if self.logger.getEffectiveLevel() <= logging.INFO:
                processedItems += 1
                statusStr = '    Finished processing %d of %d (%.2f%%) sequences.' % (processedItems, numDataItems, float(processedItems)*100/numDataItems)
                sys.stdout.write('%s\r' % statusStr)
                sys.stdout.flush()
                
            neighbours[index] = tdNeighbours
                
        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stdout.write('\n')
    
    def run(self, seqIdsOfInterest, seqs, tetraFile, tdDist, tdDistPer, threads):  
        # read in tetranucleotide signatures for each sequence
        self.logger.info('    Reading tetranucleotide signatures.')
        tetra = {}
        with open(tetraFile) as f:
            f.readline()
            for line in f:
                lineSplit = line.split('\t')
                seqId = lineSplit[0]
                if seqId in seqIdsOfInterest:
                    tetra[seqId] = np.array([float(x) for x in lineSplit[1:]])
                
        # populate worker queue with data to process
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()
        
        for i, seq in enumerate(seqs):
            workerQueue.put((i, seq))
        
        for _ in range(threads):
            workerQueue.put((None, None))
            
        tdNeighbours = mp.Manager().list([None]*len(seqs))
        
        workerProc = [mp.Process(target = self.__workerThread, args = (seqs, tetra, tdDist, tdDistPer, workerQueue, writerQueue)) for _ in range(threads)]
        writeProc = mp.Process(target = self.__writerThread, args = (tdNeighbours, len(seqs), writerQueue))
        
        writeProc.start()
        
        for p in workerProc:
            p.start()
        
        for p in workerProc:
            p.join()
        
        writerQueue.put((None, None))
        writeProc.join()

        for i, row in enumerate(tdNeighbours):
            seqs[i].tdNeighbours = row

