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
import errno
import sys
import logging


def checkFileExists(inputFile):
    if not os.path.exists(inputFile):
        logger = logging.getLogger()
        logger.error('  [Error] Input file does not exists: ' + inputFile + '\n')
        sys.exit()


def isDirEmpty(path):
    if not os.path.exists(path):
        os.makedirs(path)

    return not os.listdir(path)


def makeSurePathExists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            logger = logging.getLogger()
            logger.error('  [Error] Specified path does not exist: ' + path + '\n')
            sys.exit()


def binIdFromFilename(filename):
    binId = os.path.basename(filename)
    binId = os.path.splitext(binId)[0]

    return binId
