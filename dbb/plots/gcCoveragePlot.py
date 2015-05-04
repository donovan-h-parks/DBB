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

from AbstractPlot import AbstractPlot

import os
import random
from collections import defaultdict

from dbb.greedy import Greedy
from dbb.common import checkFileExists
from dbb.seqUtils import readFasta, readSeqStats

from matplotlib.colors import rgb2hex

import numpy as np


class GcCoveragePlot(AbstractPlot):
    def __init__(self, options):
        AbstractPlot.__init__(self, options)

    def plot(self, preprocessDir, binDir, minSeqLen, minCoreLen, bCoverageLinear, bBoundingBoxes, bLabels, legendItems=None, highlightedBinId=None, binIdToColour=None, binIdToShape=None):
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

            binId = binFile[0:binFile.rfind('.')]

            seqs = readFasta(os.path.join(binDir, binFile))
            for seqId in seqs.keys():
                binIdToSeqStats[binId].append([seqId] + scaffoldStats[seqId])
                binnedSeqs.add(seqId)

        # get scaffold stats for unbinned sequences
        for seqId, stats in scaffoldStats.iteritems():
            if seqId not in binnedSeqs:
                binIdToSeqStats[Greedy.UNBINNED].append([seqId] + stats)

        # GC vs coverage plot
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)

        axes = self.fig.add_subplot(111)

        self.plotOnAxes(axes, binIdToSeqStats, minSeqLen, minCoreLen, bCoverageLinear, bBoundingBoxes, bLabels, legendItems, highlightedBinId, binIdToColour, binIdToShape)

        self.fig.tight_layout(pad=1)

    def plotOnAxes(self, axes, binIdToSeqStats, minSeqLen, minCoreLen, bCoverageLinear, bBoundingBoxes, bLabels, legendItems, highlightedBinId, binIdToColour, binIdToShape):
        # plot each bin
        highlightedBinSize = 0
        highlightedBinGC = []
        highlightedBinCoverage = []

        legendPlotItems = []
        legendLabels = []

        for binId in sorted(binIdToSeqStats.keys()):
            gcDR = []
            covDR = []

            gcCore = []
            covCore = []
            for seqStat in binIdToSeqStats[binId]:
                seqId, seqLen, GC, coverage = seqStat
                GC = GC * 100.0

                if seqLen < minSeqLen:
                    continue

                if seqLen > minCoreLen:
                    gcCore.append(GC)
                    covCore.append(coverage)
                else:
                    gcDR.append(GC)
                    covDR.append(coverage)

                if highlightedBinId == None:
                    if binId != Greedy.UNBINNED and (not legendItems or binId in legendItems):
                        if binIdToColour:
                            color = binIdToColour[binId]
                        else:
                            color = (random.uniform(0.3, 1.0), random.uniform(0.3, 1.0), random.uniform(0.3, 1.0))
                        alpha = 0.8
                        zorder = 3
                    else:
                        color = (0.7, 0.7, 0.7)
                        alpha = 0.5
                        zorder = 1
                else:
                    if binId == highlightedBinId:
                        color = (1.0, 0.0, 0.0)
                        alpha = 1
                        zorder = 3

                        highlightedBinSize += seqLen
                        highlightedBinGC.append(GC)
                        highlightedBinCoverage.append(coverage)
                    else:
                        color = (0.7, 0.7, 0.7)
                        alpha = 1
                        zorder = 1


            colorStr = rgb2hex(color)
            if binIdToShape:
                marker = binIdToShape.get(binId, 'o')
                t = axes.scatter(gcCore, covCore, marker=marker, s=10, lw=0.5, c=colorStr, alpha=alpha, zorder=zorder + 1)
            else:
                if gcDR:
                    t = axes.scatter(gcDR, covDR, marker='s', s=10, lw=0.5, c=colorStr, alpha=alpha, zorder=zorder)
                if gcCore:
                    t = axes.scatter(gcCore, covCore, marker='o', s=10, lw=0.5, c=colorStr, alpha=alpha, zorder=zorder + 1)

            if legendItems and binId in legendItems:
                legendPlotItems.append(t)
                legendLabels.append(legendItems[binId])

            # plot bounding rect
            if highlightedBinId == None or binId == highlightedBinId:
                gc = gcDR + gcCore
                cov = covDR + covCore
                if len(gc) > 1 and binId != Greedy.UNBINNED:
                    self.boundingBox(zip(gc, cov), axes, str(binId), bBoundingBoxes, bLabels)

        # legend
        if legendItems:
            sortIndices = np.argsort(legendLabels)
            legendPlotItems = np.array(legendPlotItems)
            legendLabels = np.array(legendLabels)
            legend = axes.legend(legendPlotItems[sortIndices], legendLabels[sortIndices], scatterpoints=1, ncol=1, markerscale=2)
            legend.draw_frame(False)

        # setup axes
        axes.set_xlabel('GC (%)')
        axes.set_ylabel('Scaffold coverage')

        if not bCoverageLinear:
            axes.set_yscale('log')

        if highlightedBinId != None:
            axes.set_title('%s, # scaffolds = %d, size = %.2fMbps\nGC = %.1f +/- %.2f, Coverage = %.1f +/- %.2f' %
                            (highlightedBinId,
                             len(binIdToSeqStats[highlightedBinId]),
                             float(highlightedBinSize) / 1e6,
                             np.mean(highlightedBinGC), np.std(highlightedBinGC),
                             np.mean(highlightedBinCoverage), np.std(highlightedBinCoverage)))

        # Prettify plot
        _, yMax = axes.get_ylim()
        axes.set_ylim([-0.5, yMax])

        axes.yaxis.tick_left()
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
