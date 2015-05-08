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
import re
import logging
from collections import defaultdict

from timeKeeper import TimeKeeper

from common import isDirEmpty

from preprocess import Preprocess
from greedy import Greedy
from refineBins import RefineBins
from extractScaffolds import ExtractScaffolds
from compare import Compare
from merge import Merge
from unbinned import Unbinned
from expand import Expand
from pca import PCA
from seqUtils import readFasta, readSeqStats

from plots.gcCoveragePlot import GcCoveragePlot
from plots.tetraPcaPlot import TetraPcaPlot


class DistributionBinner(object):
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
        self.logger.info(' [DBB - preprocess] Calculate sequence statistics required for binning.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        if not isDirEmpty(args.output_dir):
            self.logger.warning('  [Warning] Specified output path must be empty: ' + args.output_dir + '\n')
            sys.exit()

        argsStr = '\n'.join(map(str.strip, str(args).replace('Namespace(', '').replace(')', '').split(',')))

        preprocess = Preprocess()
        preprocess.run(args.scaffold_file,
                       args.bam_file,
                       args.min_seq_len,
                       args.percent,
                       args.min_n,
                       args.min_align,
                       args.max_edit_dist,
                       args.bAllowImproperPairs,
                       args.coverage_type,
                       args.threads,
                       args.output_dir,
                       argsStr)

        self.timeKeeper.printTimeStamp()

    def core(self, args):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [DBB - core] Bin contigs into cores using a greedy approach.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        argsStr = ', '.join(map(str.strip, str(args).replace('Namespace(', '').replace(')', '').split(',')))

        greedy = Greedy()
        greedy.run(args.preprocess_dir, args.min_core_len, args.min_bin_size, args.build_dist_per, args.merge_dist_per, args.threads, args.binning_file, argsStr)

        self.timeKeeper.printTimeStamp()

    def refine(self, args):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [DBB - refine] Refine core bins.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        argsStr = ', '.join(map(str.strip, str(args).replace('Namespace(', '').replace(')', '').split(',')))

        refineBins = RefineBins()
        refineBins.run(args.preprocess_dir, args.binning_file, args.min_refine_len, args.gc_dist_per, args.td_dist_per, args.cov_dist_per, args.threads, args.refined_bin_file, argsStr)

        self.timeKeeper.printTimeStamp()

    def extract(self, args):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [DBB - extract] Extract scaffolds into population genome bins.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        self.logger.info('  Extracting scaffolds assigned to bins.')
        extractScaffolds = ExtractScaffolds()
        extractScaffolds.run(args.preprocess_dir, args.scaffold_file, args.binning_file, args.bin_dir)

        self.logger.info('')
        self.logger.info('  Bins written to: ' + args.bin_dir)

        self.timeKeeper.printTimeStamp()

    def binFile(self, args):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [DBB - bin_file] Creating file indicating bin assignments.')
        self.logger.info('*******************************************************************************')

        # get bin for each scaffold
        binFiles = os.listdir(args.bin_dir)
        seqIdToBinId = {}
        for binFile in binFiles:
            if not binFile.endswith('.fna'):
                continue

            binId = binFile[0:binFile.rfind('.')]

            seqs = readFasta(os.path.join(args.bin_dir, binFile))
            for seqId in seqs.keys():
                seqIdToBinId[seqId] = binId

        # write out binning file
        fout = open(args.binning_file, 'w')
        fout.write('# Created with DBB binFile command\n')
        fout.write('Contig Id\tBin Id\tLength\tGC\tCoverage\n')
        with open(os.path.join(args.preprocess_dir, 'contigs.seq_stats.tsv')) as f:
            f.readline()
            for line in f:
                lineSplit = line.split('\t')

                scaffoldId = lineSplit[0]
                match = re.search('_c[0-9]*$', scaffoldId)
                if match:
                    scaffoldId = scaffoldId[0:match.start()]

                fout.write(lineSplit[0])
                fout.write('\t' + seqIdToBinId.get(scaffoldId, 'unbinned'))
                fout.write('\t' + '\t'.join(lineSplit[1:]))

        fout.close()

        self.logger.info('')
        self.logger.info('  Binning file written to: ' + args.binning_file)
        self.timeKeeper.printTimeStamp()

    def binStats(self, args):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [DBB - bin_stats] Creating file with statistics for binned scaffolds.')
        self.logger.info('*******************************************************************************')

        # get bin for each scaffold
        binFiles = os.listdir(args.bin_dir)
        binIdToSeqIds = defaultdict(set)
        for binFile in binFiles:
            if not binFile.endswith('fna'):
                continue

            binId = binFile[0:binFile.rfind('.')]

            seqs = readFasta(os.path.join(args.bin_dir, binFile))
            for seqId in seqs.keys():
                binIdToSeqIds[binId].add(seqId)

        # write out scaffold stats
        scaffoldStats = readSeqStats(os.path.join(args.preprocess_dir, 'scaffolds.seq_stats.tsv'))
        fout = open(args.stats_file, 'w')
        fout.write('Bin Id\tScaffold Id\tLength\tGC\tCoverage\n')
        for binId, seqIds in binIdToSeqIds.iteritems():
            for seqId in seqIds:
                seqLen, GC, coverage = scaffoldStats[seqId]
                fout.write('%s\t%s\t%d\t%.2f\t%.2f\n' % (binId, seqId, seqLen, GC, coverage))
        fout.close()

        self.logger.info('')
        self.logger.info('  Statistics written to: ' + args.stats_file)
        self.timeKeeper.printTimeStamp()

    def compare(self, args):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [DBB - compare] Assess if a scaffold might belong in a bin.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        compare = Compare()
        compare.run(args.preprocess_dir, args.binning_file, args.bin_id, args.scaffold_id, args.gc_dist_per, args.td_dist_per, args.cov_dist_per)

        self.timeKeeper.printTimeStamp()

    def merge_stats(self, args):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [DBB - merge_stats] Calculate merging statistics.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        merge = Merge()
        merge.run(args.preprocess_dir, args.binning_file, args.gc_dist_per, args.td_dist_per, args.cov_dist_per, args.output_file)

        self.timeKeeper.printTimeStamp()

    def unbinned(self, args):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [DBB - unbinned] Calculate statistics for unbinned scaffolds.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        unbinned = Unbinned()
        unbinned.run(args.preprocess_dir, args.binning_file,
                     args.gc_dist_per, args.td_dist_per, args.cov_dist_per,
                     args.min_scaffold_len, args.min_contig_len,
                     args.all_scaffolds, args.output_file)

        self.timeKeeper.printTimeStamp()

    def merge(self, args):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [DBB - merge] Merge two bins.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

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
        self.logger.info(' [DBB - merge] Merge two bins.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        expand = Expand()
        expand.run(args.unbinned_stats_file, args.binning_file, args.score_threshold, args.output_file)

        self.timeKeeper.printTimeStamp()

    def tetraPcaPlot(self, args):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [DBB - tetra_pca] PCA plot of tetranucleotide signatures.')
        self.logger.info('*******************************************************************************')

        if not os.path.exists(args.plot_folder):
            os.makedirs(args.plot_folder)

        # check if PCA of tetranuclotide signatures has been computed
        pca = PCA()

        self.logger.info('')
        pcaTetraFile = os.path.join(args.preprocess_dir, 'scaffolds.pca_tetra.npz')
        if not os.path.exists(pcaTetraFile):
            self.logger.info('  Computing PCA of tetranuclotide signatures.\n')
            tetraFile = os.path.join(args.preprocess_dir, 'scaffolds.tetra.tsv')
            seqIds, pc, variance = pca.pcaFile(tetraFile, fraction=1.0, bCenter=True, bScale=False)
            pca.save(pcaTetraFile, seqIds, pc, variance)
        else:
            seqIds, pc, variance = pca.loadFile(pcaTetraFile)

        self.logger.info('  Creating PCA plot of tetranucleotide signatures.')
        plot = TetraPcaPlot(args)
        plot.plot(args.preprocess_dir, args.bin_dir, seqIds, pc, variance, args.min_seq_len, args.min_core_len, not args.hide_bounding_boxes, not args.hide_labels)

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
            binFiles = os.listdir(args.bin_dir)
            for binFile in binFiles:
                if not binFile.endswith('.fna'):
                    continue

                binId = binFile[0:binFile.rfind('.')]
                plot.plot(args.preprocess_dir, args.bin_dir, seqIds, pc, variance, args.min_seq_len, args.min_core_len, not args.hide_bounding_boxes, not args.hide_labels, binId)

                outputFile = os.path.join(args.plot_folder, 'tetra_pca_plot.' + str(binId) + '.' + args.image_type)
                plot.savePlot(outputFile, dpi=args.dpi)
                self.logger.info('    Plot written to: ' + outputFile)

        self.timeKeeper.printTimeStamp()

    def gcCovPlot(self, args):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [DBB - gc_cov] Plot GC vs coverage.')
        self.logger.info('*******************************************************************************')

        if not os.path.exists(args.plot_folder):
            os.makedirs(args.plot_folder)

        legendItems = {}
        binIdToColour = {}
        binIdToShape = {}
        if args.legend:
            for line in open(args.legend):
                lineSplit = line.split('\t')
                legendItems[lineSplit[0]] = lineSplit[1].strip()
                if len(lineSplit) == 4:
                    color = []
                    for c in lineSplit[2].split(','):
                        color.append(float(c) / 255)
                    binIdToColour[lineSplit[0]] = color

                    shape = lineSplit[3].strip()
                    binIdToShape[lineSplit[0]] = shape

        self.logger.info('')
        self.logger.info('  Creating GC vs coverage plot.')
        plot = GcCoveragePlot(args)
        plot.plot(args.preprocess_dir, args.bin_dir, args.min_seq_len, args.min_core_len, args.coverage_linear, not args.hide_bounding_boxes, not args.hide_labels, legendItems, None, binIdToColour, binIdToShape)

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

            binFiles = os.listdir(args.bin_dir)
            for binFile in binFiles:
                if binFile.endswith('.fna'):
                    continue

                binId = binFile[0:binFile.rfind('.')]
                plot.plot(args.preprocess_dir, args.bin_dir, args.min_seq_len, args.min_core_len, args.coverage_linear, not args.hide_bounding_boxes, not args.hide_labels, None, binId)

                outputFile = os.path.join(args.plot_folder, 'gc_cov_plot.' + str(binId) + '.' + args.image_type)
                plot.savePlot(outputFile, dpi=args.dpi)
                self.logger.info('    Plot written to: ' + outputFile)

        self.timeKeeper.printTimeStamp()

    def run(self, args):
        """Parse user options and call the correct pipeline(s)"""
        self.timeKeeper = TimeKeeper()

        if args.quiet:
            logging.basicConfig(format='', level=logging.ERROR)
        else:
            # logging.basicConfig(format='', level=logging.INFO)
            logging.basicConfig(format='', level=logging.DEBUG)  # prints extra information useful for
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
        elif(args.subparser_name == 'bin_wf'):
            args.output_dir = args.preprocess_dir
            self.preprocess(args)

            args.min_seq_len = args.min_refine_len

            args.binning_file = args.binning_prefix + '.core.tsv'
            args.min_bin_size = args.min_core_bin_size
            self.core(args)

            args.refined_bin_file = args.binning_prefix + '.refine.tsv'
            self.refine(args)

            args.binning_file = args.binning_prefix + '.core.tsv'
            args.bin_dir = args.bin_dir_prefix + '_core'
            self.extract(args)

            args.binning_file = args.binning_prefix + '.refine.tsv'
            args.bin_dir = args.bin_dir_prefix + '_refine'
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
        elif(args.subparser_name == 'bin_file'):
            self.binFile(args)
        elif(args.subparser_name == 'bin_stats'):
            self.binStats(args)
        elif(args.subparser_name == 'tetra_pca'):
            self.tetraPcaPlot(args)
        elif(args.subparser_name == 'gc_cov'):
            self.gcCovPlot(args)
