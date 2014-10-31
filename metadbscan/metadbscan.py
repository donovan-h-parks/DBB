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
__version__ = '0.0.3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import logging

from timeKeeper import TimeKeeper

from common import isDirEmpty

from preprocess import Preprocess
from greedy import Greedy, readSeqStatsForBins
from refineBins import RefineBins
from extractScaffolds import ExtractScaffolds
from compare import Compare
from merge import Merge
from unbinned import Unbinned
from expand import Expand
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
  
        argsStr = '\n'.join(map(str.strip, str(args).replace('Namespace(', '').replace(')', '').split(',')))
  
        preprocess = Preprocess()
        preprocess.run(args.scaffold_file, args.bam_file, args.min_seq_len, args.percent, args.min_n, args.min_align, args.max_edit_dist, args.coverage_type, args.threads, args.output_dir, argsStr)
   
        self.timeKeeper.printTimeStamp()  
        
    def core(self, args):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [MetaDBSCAN - core] Bin contigs into cores using a greedy approach.')
        self.logger.info('*******************************************************************************')
        
        argsStr = ', '.join(map(str.strip, str(args).replace('Namespace(', '').replace(')', '').split(',')))
        
        greedy = Greedy()
        greedy.run(args.preprocess_dir, args.min_seq_len, args.min_bin_size, args.build_dist_per, args.merge_dist_per, args.threads, args.binning_file, argsStr)
            
        self.timeKeeper.printTimeStamp() 
        
    def refine(self, args):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [MetaDBSCAN - refine] Refine core bins.')
        self.logger.info('*******************************************************************************')
        
        argsStr = ', '.join(map(str.strip, str(args).replace('Namespace(', '').replace(')', '').split(',')))
        
        refineBins = RefineBins()
        refineBins.run(args.preprocess_dir, args.binning_file, args.min_seq_len, args.gc_dist_per, args.td_dist_per, args.cov_dist_per, args.threads, args.refined_bin_file, argsStr)
            
        self.timeKeeper.printTimeStamp() 
        
    def extract(self, args):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [MetaDBSCAN - extract] Extract scaffolds into population genome bins.')
        self.logger.info('*******************************************************************************')
        
        self.logger.info('  Extracting scaffolds assigned to bins.')
        extractScaffolds = ExtractScaffolds()
        extractScaffolds.run(args.preprocess_dir, args.scaffold_file, args.binning_file, args.bin_dir)
        
        self.logger.info('')
        self.logger.info('  Bins written to: ' + args.bin_dir)
        
        self.timeKeeper.printTimeStamp()
        
    def compare(self, args):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [MetaDBSCAN - compare] Assess if a scaffold might belong in a bin.')
        self.logger.info('*******************************************************************************')
        
        compare = Compare()
        compare.run(args.preprocess_dir, args.binning_file, args.bin_id, args.scaffold_id, args.gc_dist_per, args.td_dist_per, args.cov_dist_per)
        
        self.timeKeeper.printTimeStamp()     
        
    def merge_stats(self, args):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [MetaDBSCAN - merge_stats] Calculate merging statistics.')
        self.logger.info('*******************************************************************************')
        
        merge = Merge()
        merge.run(args.preprocess_dir, args.binning_file, args.gc_dist_per, args.td_dist_per, args.cov_dist_per, args.output_file)
        
        self.timeKeeper.printTimeStamp()    
        
    def unbinned(self, args):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [MetaDBSCAN - unbinned] Calculate statistics for unbinned scaffolds.')
        self.logger.info('*******************************************************************************')
        
        unbinned = Unbinned()
        unbinned.run(args.preprocess_dir, args.binning_file, args.gc_dist_per, args.td_dist_per, args.cov_dist_per, args.min_scaffold_len, args.min_contig_len, args.output_file)
        
        self.timeKeeper.printTimeStamp()    
        
    def merge(self, args):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [MetaDBSCAN - merge] Merge two bins.')
        self.logger.info('*******************************************************************************')
        
        # read bins file
        header = None
        rows = []
        for line in open(args.binning_file):
            if line[0] == '#':
                continue
            
            if not header:
                header = line
                continue
            
            rows.append(line)
            
        # write out modified binning file
        fout = open(args.output_file, 'w')
        fout.write(header)
        for row in rows:
            rowSplit = row.split('\t')
            
            if rowSplit[1] == 'unbinned':
                fout.write(row)
            else:   
                binId = int(rowSplit[1])
                if binId == args.bin1_id or binId == args.bin2_id:
                    rowSplit[1] = str(args.merge_id)
                    
                fout.write('\t'.join(rowSplit))
        fout.close()
        
        self.timeKeeper.printTimeStamp()    
        
    def expand(self, args):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [MetaDBSCAN - merge] Merge two bins.')
        self.logger.info('*******************************************************************************')
        
        expand = Expand()
        expand.run(args.unbinned_stats_file, args.binning_file, args.score_threshold, args.output_file)
        
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
            seqStatsForBins = readSeqStatsForBins(args.binning_file)
            for binId in sorted(seqStatsForBins.keys()):
                if binId != Greedy.UNBINNED:
                    plot.plot(args.binning_file, seqIds, pc, variance, args.min_seq_len, args.min_core_len, not args.no_bounding_box, not args.no_labels, binId)
                    
                    outputFile = os.path.join(args.plot_folder, 'tetra_pca_plot.bin' + str(binId) + '.' + args.image_type)
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
            
        legendItems = {}
        clusterIdToColour = {}
        clusterIdToShape = {}
        if args.legend:
            for line in open(args.legend):
                lineSplit = line.split('\t')
                legendItems[int(lineSplit[0])] = lineSplit[1].strip()
                if len(lineSplit) == 4:
                    color = []
                    for c in lineSplit[2].split(','):
                        color.append(float(c)/255)
                    clusterIdToColour[int(lineSplit[0])] = color
                    
                    shape = lineSplit[3].strip()
                    clusterIdToShape[int(lineSplit[0])] = shape
        
        self.logger.info('  Creating GC vs coverage plot.')
        plot = GcCoveragePlot(args)      
        plot.plot(args.binning_file, args.min_seq_len, args.min_core_len, args.coverage_linear, not args.hide_bounding_boxes, not args.hide_labels, legendItems, None, clusterIdToColour, clusterIdToShape)
        
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
            seqStatsForBins = readSeqStatsForBins(args.binning_file)
            for binId in sorted(seqStatsForBins.keys()):
                if binId != Greedy.UNBINNED:
                    plot.plot(args.binning_file, args.min_seq_len, args.min_core_len, args.coverage_linear, not args.hide_bounding_boxes, not args.hide_labels, None, binId)
                    
                    outputFile = os.path.join(args.plot_folder, 'gc_cov_plot.bin' + str(binId) + '.' + args.image_type)
                    plot.savePlot(outputFile, dpi=args.dpi)
                    self.logger.info('    Plot written to: ' + outputFile)
            
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
        elif(args.subparser_name == 'core'):
            self.core(args)
        elif(args.subparser_name == 'refine'):
            self.refine(args)
        elif(args.subparser_name == 'extract'):
            self.extract(args)
        elif(args.subparser_name == 'compare'):
            self.compare(args)
        elif(args.subparser_name == 'merge_stats'):
            self.merge_stats(args)
        elif(args.subparser_name == 'unbinned'):
            self.unbinned(args)
        elif(args.subparser_name == 'merge'):
            self.merge(args)
        elif(args.subparser_name == 'expand'):
            self.expand(args)
        elif(args.subparser_name == 'tetra_pca'):
            self.tetraPcaPlot(args)
        elif(args.subparser_name == 'gc_cov'):
            self.gcCovPlot(args)
