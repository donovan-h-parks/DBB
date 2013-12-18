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

__prog_desc__ = 'bin contigs into population genomes'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import logging

from timeKeeper import TimeKeeper

from common import isDirEmpty, checkFileExists

from preprocess import Preprocess
from dbscan import DBSCAN, readSeqStatsForClusters
from refineBins import RefineBins
from extractScaffolds import ExtractScaffolds
from splitScaffolds import SplitScaffolds
from pca import PCA

from plots.gcCoveragePlot import GcCoveragePlot
from plots.tetraPcaPlot import TetraPcaPlot

class MetaDBSCAN(object):
    def __init__(self):
        self.logger = logging.getLogger()
        
    def binFiles(self, options):
        targetFiles = []
        if options.bin_folder is not None:
            all_files = os.listdir(options.bin_folder)
            for f in all_files:
                if f.endswith(options.extension):
                    targetFiles.append(os.path.join(options.bin_folder, f))
                    
        if not targetFiles:
            self.logger.error("  [Error] No bins found. Check the extension used to identify bins.")
            sys.exit()
        
        return sorted(targetFiles)
    
    def preprocess(self, args):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [MetaDBSCAN - preprocess] Calculate sequence statistics required for binning.')
        self.logger.info('*******************************************************************************')
  
        if not isDirEmpty(args.output_dir):
            self.logger.warning('  [Warning] Specified output path must be empty: ' + args.output_dir + '\n')
            sys.exit()
  
        preprocess = Preprocess()
        preprocess.run(args.contig_file, args.bam_file, args.min_seq_len, args.percent, args.block_size, args.threads, args.output_dir)
   
        self.timeKeeper.printTimeStamp()  
        
    def bin(self, args):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [MetaDBSCAN - bin] Bin contig partitions.')
        self.logger.info('*******************************************************************************')
        
        dbscan = DBSCAN()
        dbscan.run(args.preprocess_dir, args.min_seq_len, args.min_bin_size, args.gc_dist_per, args.td_dist_per, args.cov_dist_per, args.min_pts, args.min_core_len, args.threads, args.binning_file)
            
        self.timeKeeper.printTimeStamp() 
        
    def refine(self, args):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [MetaDBSCAN - refine] Refine bins.')
        self.logger.info('*******************************************************************************')
        
        refineBins = RefineBins()
        refineBins.run(args.preprocess_dir, args.binning_file, args.min_seq_len, args.gc_dist_per, args.td_dist_per, args.cov_dist_per, args.refined_bin_file)
            
        self.timeKeeper.printTimeStamp() 
        
    def extract(self, args):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [MetaDBSCAN - extract] Extract scaffolds into population genome bins.')
        self.logger.info('*******************************************************************************')
        
        self.logger.info('  Extracting scaffolds assigned to each bin.')
        extractScaffolds = ExtractScaffolds()
        extractScaffolds.run(args.scaffold_file, args.binning_file, args.bin_dir)
        
        self.logger.info('')
        self.logger.info('  Bins written to: ' + args.bin_dir)
        
        self.timeKeeper.printTimeStamp() 
        
    def tetraPcaPlot(self, args):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [MetaDBSCAN - tetra_pca] PCA plot of tetranucleotide signatures.')
        self.logger.info('*******************************************************************************')
  
        if not os.path.exists(args.plot_folder):
            os.makedirs(args.plot_folder)
        
        # check if PCA of tetranuclotide signatures has been computed
        pca = PCA()
        
        pcaTetraFile = os.path.join(args.preprocess_dir, 'partitions.pca_tetra.npz')
        if not os.path.exists(pcaTetraFile):
            self.logger.info('  Computing PCA of tetranuclotide signatures.\n')
            tetraFile = os.path.join(args.preprocess_dir, 'partitions.tetra.tsv')
            seqIds, pc, variance = pca.pcaFile(tetraFile, fraction=1.0, bCenter=True, bScale=False)
            pca.save(pcaTetraFile, seqIds, pc, variance)
        else:
            seqIds, pc, variance = pca.loadFile(pcaTetraFile)
            
        self.logger.info('  Creating PCA plot of tetranucleotide signatures.')
        plot = TetraPcaPlot(args)      
        plot.plot(args.binning_file, seqIds, pc, variance, args.min_seq_len, args.min_core_len, not args.no_bounding_box, not args.no_labels)
        
        outputFile = os.path.join(args.plot_folder, 'tetra_pca_plot.' + args.image_type)
        plot.savePlot(outputFile, dpi=args.dpi)
        if args.individual:
            self.logger.info('    Plot written to: ' + outputFile)
        else:
            self.logger.info('')
            self.logger.info('  Plot written to: ' + outputFile)
        
        if args.individual:
            self.logger.info('')
            self.logger.info('  Creating PCA plot of tetranucleotide signatures for individual bins.')
            seqStatsForClusters = readSeqStatsForClusters(args.binning_file)
            for clusterId in sorted(seqStatsForClusters.keys()):
                if clusterId != DBSCAN.NOISE:
                    plot.plot(args.binning_file, seqIds, pc, variance, args.min_seq_len, args.min_core_len, not args.no_bounding_box, not args.no_labels, clusterId)
                    
                    outputFile = os.path.join(args.plot_folder, 'tetra_pca_plot.bin' + str(clusterId) + '.' + args.image_type)
                    plot.savePlot(outputFile, dpi=args.dpi)
                    self.logger.info('    Plot written to: ' + outputFile)
            
        self.timeKeeper.printTimeStamp()
        
    def gcCovPlot(self, args):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [MetaDBSCAN - gc_cov] Plot GC vs coverage.')
        self.logger.info('*******************************************************************************')
             
        if not os.path.exists(args.plot_folder):
            os.makedirs(args.plot_folder)
        
        self.logger.info('  Creating GC vs coverage plot.')
        plot = GcCoveragePlot(args)      
        plot.plot(args.binning_file, args.min_seq_len, args.min_core_len, args.coverage_log, not args.no_bounding_box, not args.no_labels)
        
        outputFile = os.path.join(args.plot_folder, 'gc_cov_plot.' + args.image_type)
        plot.savePlot(outputFile, dpi=args.dpi)
        if args.individual:
            self.logger.info('    Plot written to: ' + outputFile)
        else:
            self.logger.info('')
            self.logger.info('  Plot written to: ' + outputFile)
        
        if args.individual:
            self.logger.info('')
            self.logger.info('  Creating GC vs coverage plot for individual bins.')
            seqStatsForClusters = readSeqStatsForClusters(args.binning_file)
            for clusterId in sorted(seqStatsForClusters.keys()):
                if clusterId != DBSCAN.NOISE:
                    plot.plot(args.binning_file, args.min_seq_len, args.min_core_len, args.coverage_log, not args.no_bounding_box, not args.no_labels, clusterId)
                    
                    outputFile = os.path.join(args.plot_folder, 'gc_cov_plot.bin' + str(clusterId) + '.' + args.image_type)
                    plot.savePlot(outputFile, dpi=args.dpi)
                    self.logger.info('    Plot written to: ' + outputFile)
            
        self.timeKeeper.printTimeStamp()
        
    def splitScaffolds(self, args):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [MetaDBSCAN - split_scaf] Split scaffolds into contigs.')
        self.logger.info('*******************************************************************************')
            
        checkFileExists(args.scaffold_file)
        
        splitScaffolds = SplitScaffolds()
        splitScaffolds.run(args.scaffold_file, args.min_n, args.output_file)

        self.logger.info('')
        self.logger.info('  Plot written to: ' + args.output_file)
            
        self.timeKeeper.printTimeStamp()
     
    def run(self, args):
        """Parse user options and call the correct pipeline(s)"""
        self.timeKeeper = TimeKeeper()

        if args.quiet:
            logging.basicConfig(format='', level=logging.ERROR)
        else:
            #logging.basicConfig(format='', level=logging.INFO)
            logging.basicConfig(format='', level=logging.DEBUG) # prints extra information useful for 
                                                                # debugging or algorithm development
            
        # execute desired command
        if(args.subparser_name == 'preprocess'):
            self.preprocess(args)
        elif(args.subparser_name == 'bin'):
            self.bin(args)
        elif(args.subparser_name == 'refine'):
            self.refine(args)
        elif(args.subparser_name == 'extract'):
            self.extract(args)
        elif(args.subparser_name == 'tetra_pca'):
            self.tetraPcaPlot(args)
        elif(args.subparser_name == 'gc_cov'):
            self.gcCovPlot(args)
        elif(args.subparser_name == 'split_scaf'):
            self.splitScaffolds(args)
