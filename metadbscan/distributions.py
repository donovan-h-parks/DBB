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

__prog_desc__ = 'functions for working with GC and TD distributions'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import ast

import numpy as np

def readDistribution(distPrefix, distPer):
    """Read distribution file."""
    distFile = os.path.join(os.path.join('/srv/whitlam/home/users/uqdparks/git/MetaDBSCAN/data', distPrefix, 'distribution_' + str(distPer) + '.txt'))

    with open(distFile) as f:
        s = f.read()
        d = ast.literal_eval(s)
        
    return d

def findNearest(array, value):
    '''Find nearest array element to a given value.'''
    idx = (np.abs(np.array(array)-value)).argmin()
    return array[idx]

def inRange(testValue, lowerBound, upperBound):
    '''Determine is test value falls in a given range.'''
    if testValue >= lowerBound and testValue <= upperBound:
        return True
    
    return False
