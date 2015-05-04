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
import ast
import logging
import numpy as np

from common import checkFileExists


def readDistributions(preprocessDir, gcDistPer, tdDistPer, covDistPer):
    gcDistFile = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'data', 'gc_dist.txt')
    checkFileExists(gcDistFile)

    tdDistFile = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'data', 'td_dist.txt')
    checkFileExists(tdDistFile)

    coverageDistFile = os.path.join(preprocessDir, 'coverage_dist.txt')
    checkFileExists(coverageDistFile)

    with open(gcDistFile, 'r') as f:
        s = f.read()
        gcDist = ast.literal_eval(s)

    with open(tdDistFile, 'r') as f:
        s = f.read()
        tdDist = ast.literal_eval(s)

    with open(coverageDistFile) as f:
        s = f.read()
        covDist = ast.literal_eval(s)

    return gcDist, tdDist, covDist


class Distributions(object):
    def __init__(self, gcDist, gcDistPer, tdDist, tdDistPer, covDist, covDistPer):
        self.logger = logging.getLogger()

        # determine distribution index values for each distribution
        # (this is a little messy due to floating point keys)
        if gcDist != None:
            self.gcDist = gcDist
            sampleMeanGC = gcDist.keys()[0]
            sampleSeqLen = gcDist[sampleMeanGC].keys()[0]
            d = gcDist[sampleMeanGC][sampleSeqLen]
            self.gcLowerBoundKey = self.__findNearest(d.keys(), (100 - gcDistPer) / 2.0)
            self.gcUpperBoundKey = self.__findNearest(d.keys(), (100 + gcDistPer) / 2.0)

        if tdDist != None:
            self.tdDist = tdDist
            self.tdKey = self.__findNearest(tdDist[tdDist.keys()[0]].keys(), tdDistPer)

        if covDist != None:
            self.covDist = covDist
            d = covDist[covDist.keys()[0]]
            self.covLowerBoundKey = self.__findNearest(d.keys(), (100 - covDistPer) / 2.0)
            self.covUpperBoundKey = self.__findNearest(d.keys(), (100 + covDistPer) / 2.0)

    def boundsDistGC(self, unbinnedContig, core, fixedSeqLen=None):
        closestMeanGC = self.__findNearest(self.gcDist.keys(), core.GC)

        if fixedSeqLen is None:
            closestSeqLen = self.__findNearest(self.gcDist[closestMeanGC].keys(), unbinnedContig.length)
        else:
            closestSeqLen = self.__findNearest(self.gcDist[closestMeanGC].keys(), fixedSeqLen)

        d = self.gcDist[closestMeanGC][closestSeqLen]
        lowerBound = d[self.gcLowerBoundKey]
        upperBound = d[self.gcUpperBoundKey]

        return lowerBound, upperBound

    def withinDistGC(self, unbinnedContig, core, fixedSeqLen=None):
        lowerBound, upperBound = self.boundsDistGC(unbinnedContig, core, fixedSeqLen)

        gcDiff = unbinnedContig.GC - core.GC
        return self.__inRange(gcDiff, lowerBound, upperBound)

    def boundDistTD(self, unbinnedContig, fixedSeqLen=None):
        if fixedSeqLen is None:
            closestSeqLen = self.__findNearest(self.tdDist.keys(), unbinnedContig.length)
        else:
            closestSeqLen = self.__findNearest(self.tdDist.keys(), fixedSeqLen)

        return self.tdDist[closestSeqLen][self.tdKey]

    def witinDistTD(self, tetraDist, unbinnedContig, fixedSeqLen=None):
        return tetraDist < self.boundDistTD(unbinnedContig, fixedSeqLen)

    def boundsDistCov(self, unbinnedContig, core, fixedSeqLen=None):
        if fixedSeqLen is None:
            closestSeqLen = self.__findNearest(self.covDist.keys(), unbinnedContig.length)
        else:
            closestSeqLen = self.__findNearest(self.covDist.keys(), fixedSeqLen)

        d = self.covDist[closestSeqLen]
        lowerBound = d[self.covLowerBoundKey]
        upperBound = d[self.covUpperBoundKey]

        return lowerBound, upperBound

    def withinDistCov(self, unbinnedContig, core, fixedSeqLen=None):
        lowerBound, upperBound = self.boundsDistCov(unbinnedContig, core, fixedSeqLen)

        if core.coverage > 0:
            covPerDiff = (unbinnedContig.coverage - core.coverage) * 100.0 / core.coverage
        else:
            # assume a coverage of 0.1
            covPerDiff = (unbinnedContig.coverage - 0.1) * 100.0 / 0.1

        return self.__inRange(covPerDiff, lowerBound, upperBound)

    def gcDistance(self, unbinnedContig, core):
        return abs(unbinnedContig.GC - core.GC)

    def covDistance(self, unbinnedContig, core):
        return abs(unbinnedContig.coverage - core.coverage)

    def __findNearest(self, array, value):
        '''Find nearest array element to a given value.'''
        idx = (np.abs(np.array(array) - value)).argmin()
        return array[idx]

    def __inRange(self, testValue, lowerBound, upperBound):
        '''Determine is test value falls in a given range.'''
        if testValue >= lowerBound and testValue <= upperBound:
            return True

        return False
