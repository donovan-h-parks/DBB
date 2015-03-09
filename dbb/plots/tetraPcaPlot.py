###############################################################################
#
# pcaPlot.py - create plots of PCA'ed data matrix
#
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
import random
from collections import defaultdict

from AbstractPlot import AbstractPlot
from gcCoveragePlot import GcCoveragePlot

from dbb.common import checkFileExists
from dbb.seqUtils import readFasta, readSeqStats
from dbb.greedy import Greedy

import numpy as np

from matplotlib.colors import rgb2hex


class TetraPcaPlot(AbstractPlot):
    def __init__(self, options):
        AbstractPlot.__init__(self, options)

    def plot(self, preprocessDir, binDir, seqIds, pc, variance, minSeqLen, minCoreLen, bBoundingBoxes, bLabels, highlightedClusterId=None):
        # read statistics for each scaffold
        scaffoldStatsFile = os.path.join(preprocessDir, 'scaffolds.seq_stats.tsv')
        checkFileExists(scaffoldStatsFile)
        scaffoldStats = readSeqStats(scaffoldStatsFile)

        # read scaffolds stats for binned sequences
        binFiles = os.listdir(binDir)
        binIdToSeqStats = defaultdict(list)
        binnedSeqs = set()
        for binFile in binFiles:
            if not binFile.endswith('.fna'):
                continue

            binId = binFile[0:binFile.find('.')]

            seqs = readFasta(os.path.join(binDir, binFile))
            for seqId in seqs.keys():
                binIdToSeqStats[binId].append([seqId] + scaffoldStats[seqId])
                binnedSeqs.add(seqId)

        # get scaffold stats for unbinned sequences
        for seqId, stats in scaffoldStats.iteritems():
            if seqId not in binnedSeqs:
                binIdToSeqStats[Greedy.UNBINNED].append([seqId] + stats)

        # ensure pc matrix has at least 3 dimensions
        if pc.shape[1] == 1:
            pc = np.append(pc, np.zeros((pc.shape[0], 2)), 1)
        elif pc.shape[1] == 2:
            pc = np.append(pc, np.zeros((pc.shape[0], 1)), 1)

        # Set size of figure
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)

        axesPC1vsPC2 = self.fig.add_subplot(221)
        axesPC3vsPC2 = self.fig.add_subplot(222)
        axesPC1vsPC3 = self.fig.add_subplot(223)
        axesGcCov = self.fig.add_subplot(224)

        # plot each bin
        clusterIdToColour = {}
        seqIds = list(seqIds)
        for binId in sorted(binIdToSeqStats.keys()):
            drPC1 = []
            drPC2 = []
            drPC3 = []

            corePC1 = []
            corePC2 = []
            corePC3 = []
            for seqStat in binIdToSeqStats[binId]:
                seqId, seqLen, GC, coverage = seqStat
                try:
                    seqIndex = seqIds.index(seqId)
                except:
                    continue

                if seqLen > minCoreLen:
                    corePC1.append(pc[seqIndex, 0])
                    corePC2.append(pc[seqIndex, 1])
                    corePC3.append(pc[seqIndex, 2])
                else:
                    drPC1.append(pc[seqIndex, 0])
                    drPC2.append(pc[seqIndex, 1])
                    drPC3.append(pc[seqIndex, 2])

            if highlightedClusterId == None:
                if binId != Greedy.UNBINNED:
                    color = (random.uniform(0.3, 1.0), random.uniform(0.3, 1.0), random.uniform(0.3, 1.0))
                    clusterIdToColour[binId] = color
                    alpha = 0.7
                    zorder = 3
                else:
                    color = (0.7, 0.7, 0.7)
                    alpha = 0.5
                    zorder = 1
            else:
                if binId == highlightedClusterId:
                    color = (1.0, 0.0, 0.0)
                    alpha = 0.7
                    zorder = 3
                else:
                    color = (0.7, 0.7, 0.7)
                    alpha = 0.5
                    zorder = 1

            colorStr = rgb2hex(color)
            if drPC1:
                axesPC1vsPC2.scatter(drPC1, drPC2, marker='o', s=10, lw=0.5, c=colorStr, alpha=alpha, zorder=zorder)
                axesPC3vsPC2.scatter(drPC3, drPC2, marker='o', s=10, lw=0.5, c=colorStr, alpha=alpha, zorder=zorder)
                axesPC1vsPC3.scatter(drPC1, drPC3, marker='o', s=10, lw=0.5, c=colorStr, alpha=alpha, zorder=zorder)

            if corePC1:
                axesPC1vsPC2.scatter(corePC1, corePC2, marker='s', s=10, lw=0.5, c=colorStr, alpha=alpha, zorder=zorder + 1)
                axesPC3vsPC2.scatter(corePC3, corePC2, marker='s', s=10, lw=0.5, c=colorStr, alpha=alpha, zorder=zorder + 1)
                axesPC1vsPC3.scatter(corePC1, corePC3, marker='s', s=10, lw=0.5, c=colorStr, alpha=alpha, zorder=zorder + 1)

            # plot bounding rect
            if highlightedClusterId == None or binId == highlightedClusterId:
                pc1 = drPC1 + corePC1
                pc2 = drPC2 + corePC2
                pc3 = drPC3 + corePC3
                if len(pc1) > 1 and binId != Greedy.UNBINNED:
                    self.boundingBox(zip(pc1, pc2), axesPC1vsPC2, str(binId), bBoundingBoxes, bLabels)
                    self.boundingBox(zip(pc3, pc2), axesPC3vsPC2, str(binId), bBoundingBoxes, bLabels)
                    self.boundingBox(zip(pc1, pc3), axesPC1vsPC3, str(binId), bBoundingBoxes, bLabels)

        # set labels
        axesPC1vsPC2.set_xlabel('PC1 (%.1f%%)' % (variance[0] * 100))
        axesPC1vsPC2.set_ylabel('PC2 (%.1f%%)' % (variance[1] * 100))
        axesPC3vsPC2.set_xlabel('PC3 (%.1f%%)' % (variance[2] * 100))
        axesPC3vsPC2.set_ylabel('PC2 (%.1f%%)' % (variance[1] * 100))
        axesPC1vsPC3.set_xlabel('PC1 (%.1f%%)' % (variance[0] * 100))
        axesPC1vsPC3.set_ylabel('PC3 (%.1f%%)' % (variance[2] * 100))

        # Prettify plot
        for axes in [axesPC1vsPC2, axesPC3vsPC2, axesPC1vsPC3]:
            for a in axes.yaxis.majorTicks:
                a.tick1On = True
                a.tick2On = False

            for a in axes.xaxis.majorTicks:
                a.tick1On = True
                a.tick2On = False

            for line in axes.yaxis.get_ticklines():
                line.set_color(self.axesColour)

            for line in axes.xaxis.get_ticklines():
                line.set_color(self.axesColour)

            for loc, spine in axes.spines.iteritems():
                if loc in ['right', 'top']:
                    spine.set_color('none')
                else:
                    spine.set_color(self.axesColour)

        # plot GC vs. coverage
        gcCoveragePlot = GcCoveragePlot(self.options)
        gcCoveragePlot.plotOnAxes(axesGcCov, binIdToSeqStats, minSeqLen, minCoreLen, False, bBoundingBoxes, bLabels, None, highlightedClusterId, clusterIdToColour, None)

        self.fig.tight_layout(pad=2, w_pad=2, h_pad=2)