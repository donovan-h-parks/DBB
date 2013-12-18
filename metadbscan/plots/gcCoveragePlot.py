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

import random

from metadbscan.dbscan import DBSCAN, readSeqStatsForClusters

from matplotlib.colors import rgb2hex

class GcCoveragePlot(AbstractPlot):
    def __init__(self, options):
        AbstractPlot.__init__(self, options)

    def plot(self, binningFile, minSeqLen, minCoreLen, bCoverageLog, bBoundingBoxes, bLabels, highlightedClusterId = None, clusterIdToColour = None):
        # read GC and coverage for sequences
        seqStatsForClusters = readSeqStatsForClusters(binningFile)
        
        # GC vs coverage plot
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)
        
        axes = self.fig.add_subplot(111)
        
        self.plotOnAxes(axes, seqStatsForClusters, minSeqLen, minCoreLen, bCoverageLog, bBoundingBoxes, bLabels, highlightedClusterId, clusterIdToColour)
        
        self.fig.tight_layout(pad=5)
        
    def plotOnAxes(self, axes, seqStatsForClusters, minSeqLen, minCoreLen, bCoverageLog, bBoundingBoxes, bLabels, highlightedClusterId, clusterIdToColour):
        # plot each bin
        for clusterId in sorted(seqStatsForClusters.keys()):
            gcDR = []
            covDR = []
            
            gcCore = []
            covCore = []
            for seqStat in seqStatsForClusters[clusterId]:
                _, seqLen, GC, coverage = seqStat
                if seqLen > minCoreLen:
                    gcCore.append(GC)
                    covCore.append(coverage)
                else:
                    gcDR.append(GC)
                    covDR.append(coverage)
                
            if highlightedClusterId == None:
                if clusterId != DBSCAN.NOISE:
                    if clusterIdToColour != None:
                        color = clusterIdToColour[clusterId]
                    else:
                        color = (random.uniform(0.3, 1.0), random.uniform(0.3, 1.0), random.uniform(0.3, 1.0))
                    alpha = 0.7
                    zorder = 3
                else:
                    color = (0.7, 0.7, 0.7)
                    alpha = 0.5
                    zorder = 1
            else:
                if clusterId == highlightedClusterId:
                    color = (1.0, 0.0, 0.0)
                    alpha = 1
                    zorder = 3
                else:
                    color = (0.7, 0.7, 0.7)
                    alpha = 1
                    zorder = 1
                
            colorStr = rgb2hex(color)
            if gcDR:
                axes.scatter(gcDR, covDR, marker='o', s=10, lw=0.5, c=colorStr, alpha=alpha, zorder=zorder)

            if gcCore:
                axes.scatter(gcCore, covCore, marker='s', s=10, lw=0.5, c=colorStr, alpha=alpha, zorder=zorder+1)
                
            # plot bounding rect  
            if highlightedClusterId == None or clusterId == highlightedClusterId:
                gc = gcDR + gcCore
                cov = covDR + covCore
                if len(gc) > 1 and clusterId != DBSCAN.NOISE:
                    self.boundingBox(zip(gc, cov), axes, str(clusterId), bBoundingBoxes, bLabels)
                
        # setup axes 
        axes.set_xlabel('GC')
        axes.set_ylabel('Coverage')
        
        if bCoverageLog:
            axes.set_yscale('log')
            axes.set_ylabel('log(coverage)')
        
        # Prettify plot
        _, yMax = axes.get_ylim()
        axes.set_ylim([-0.5, yMax])
        
        axes.yaxis.tick_left()
        for a in axes.yaxis.majorTicks:
            a.tick1On=True
            a.tick2On=False
        
        for a in axes.xaxis.majorTicks:
            a.tick1On=True
            a.tick2On=False
            
        for line in axes.yaxis.get_ticklines():
            line.set_color(self.axesColour)
        
        for line in axes.xaxis.get_ticklines():
            line.set_color(self.axesColour)
        
        for loc, spine in axes.spines.iteritems():
            if loc in ['right','top']:
                spine.set_color('none')
            else:
                spine.set_color(self.axesColour)
