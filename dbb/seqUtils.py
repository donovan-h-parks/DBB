###############################################################################
#
# seqUtils.py - Common functions for interacting with sequences
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
import sys
import gzip
import logging
import exceptions as exception


class InputFileException(exception.Exception):
    pass


def readSeqStats(seqStatsFile):
    '''Read file containing sequence statistics.'''
    seqStats = {}
    with open(seqStatsFile) as f:
        f.readline()
        for line in f:
            lineSplit = line.split('\t')
            seqStats[lineSplit[0]] = [int(lineSplit[1]), float(lineSplit[2]), float(lineSplit[3])]

    return seqStats


def readFasta(fastaFile):
    '''Read sequences from FASTA file.'''
    try:
        if fastaFile.endswith('.gz'):
            openFile = gzip.open
        else:
            openFile = open

        seqs = {}
        for line in openFile(fastaFile):
            if line[0] == '>':
                seqId = line[1:].partition(' ')[0].rstrip()
                seqs[seqId] = []
            else:
                seqs[seqId].append(line[0:-1])

        for seqId, seq in seqs.iteritems():
            seqs[seqId] = ''.join(seq)
    except:
        logger = logging.getLogger()
        logger.error("  [Error] Failed to process sequence file: " + fastaFile)
        sys.exit()

    return seqs


def readSeq(seq_file):
    """Generator function to read sequences from fasta/q file.

    This function is intended to be used as a generator
    in order to avoid having to have large sequence files
    in memory. Input file may be gzipped and in either
    fasta or fastq format. It is slightly more efficient
    to directly call read_fasta_seq() or read_fastq_seq()
    if the type of input file in known.

    Example:
    seq_io = SeqIO()
    for seq_id, seq in seq_io.read_seq(fasta_file):
        print seq_id
        print seq

    Parameters
    ----------
    seq_file : str
        Name of fasta/q file to read.

    Yields
    ------
    list : [seq_id, seq]
        Unique id of the sequence followed by the sequence itself.
    """

    if seq_file.endswith(('.fq.gz', '.fastq.gz', '.fq', '.fq.gz')):
        for rtn in readFastqSeq(seq_file):
            yield rtn
    else:
        for rtn in readFastaSeq(seq_file):
            yield rtn


def readFastaSeq(fasta_file):
    """Generator function to read sequences from fasta file.

    This function is intended to be used as a generator
    in order to avoid having to have large sequence files
    in memory. Input file may be gzipped.

    Example:
    seq_io = SeqIO()
    for seq_id, seq in seq_io.read_fasta_seq(fasta_file):
        print seq_id
        print seq

    Parameters
    ----------
    fasta_file : str
        Name of fasta file to read.

    Yields
    ------
    list : [seq_id, seq]
        Unique id of the sequence followed by the sequence itself.
    """

    if not os.path.exists(fasta_file) or os.stat(fasta_file).st_size == 0:
        raise InputFileException

    try:
        open_file = open
        if fasta_file.endswith('.gz'):
            open_file = gzip.open

        seq_id = None
        seq = None
        for line in open_file(fasta_file):
            # skip blank lines
            if not line.strip():
                continue

            if line[0] == '>':
                if seq_id != None:
                    yield seq_id, ''.join(seq)

                seq_id = line[1:].split(None, 1)[0]
                seq = []
            else:
                seq.append(line[0:-1])

        # report last sequence
        yield seq_id, ''.join(seq)
    except:
        print traceback.format_exc()
        print ''
        print  "[Error] Failed to process sequence file: " + fasta_file
        sys.exit()


def readFastqSeq(fastq_file):
    """Generator function to read sequences from fastq file.

    This function is intended to be used as a generator
    in order to avoid having to have large sequence files
    in memory. Input file may be gzipped.

    Example:
    seq_io = SeqIO()
    for seq_id, seq in seq_io.read_fastq_seq(fastq_file):
        print seq_id
        print seq

    Parameters
    ----------
    fastq_file : str
        Name of fastq file to read.

    Yields
    ------
    list : [seq_id, seq]
        Unique id of the sequence followed by the sequence itself.
    """

    if not os.path.exists(fastq_file) or os.stat(fastq_file).st_size == 0:
        raise InputFileException

    try:
        open_file = open
        if fastq_file.endswith('.gz'):
            open_file = gzip.open

        line_num = 0
        for line in open_file(fastq_file):
            line_num += 1

            if line_num == 1:
                seq_id = line[1:].split(None, 1)[0]
            elif line_num == 2:
                yield seq_id, line[0:-1]
            elif line_num == 4:
                line_num = 0
    except:
        print traceback.format_exc()
        print ''
        print  "[Error] Failed to process sequence file: " + fastq_file
        sys.exit()


def readFastaSeqIds(fastaFile):
    '''Read sequence ids from FASTA file.'''
    seqIds = []
    for line in open(fastaFile):
        if line[0] == '>':
            seqId = line[1:].partition(' ')[0].rstrip()
            seqIds.append(seqId)

    return seqIds


def readGenomicSeqsFromFasta(fastaFile, seqToIgnore=None):
    '''Read genomic sequences from FASTA file. Explicitly ignores sequences marked as plasmids.'''
    seqs = {}
    bRead = False
    for line in open(fastaFile):
        if line[0] == '>':
            if 'plasmid' in line.lower():
                bRead = False
            else:
                seqId = line[1:].partition(' ')[0]
                seqs[seqId] = []
                bRead = True
        elif bRead:
            seqs[seqId].append(line[0:-1])

    for seqId, seq in seqs.iteritems():
        seqs[seqId] = ''.join(seq)

    return seqs


def writeFasta(seqs, outputFile):
    '''Write sequences from FASTA file.'''
    if outputFile.endswith('.gz'):
        fout = gzip.open(outputFile, 'wb')
    else:
        fout = open(outputFile, 'w')

    for seqId, seq in seqs.iteritems():
        fout.write('>' + seqId + '\n')
        fout.write(seq + '\n')
    fout.close()


def baseCount(seq):
    testSeq = seq.upper()
    a = testSeq.count('A')
    c = testSeq.count('C')
    g = testSeq.count('G')
    t = testSeq.count('T') + testSeq.count('U')

    return a, c, g, t


def calculateN50(seqLens):
    thresholdN50 = sum(seqLens) / 2.0

    seqLens.sort(reverse=True)

    testSum = 0
    for seqLen in seqLens:
        testSum += seqLen
        if testSum >= thresholdN50:
            N50 = seqLen
            break

    return N50