##########################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
##########################################################################
"""Utils.py - Utilities for ruffus pipelines
============================================

Reference
---------

"""
import inspect
import sys


def isTest():
    '''return True if the pipeline is run in a "testing" mode
    (command line options --is-test has been given).'''

    # note: do not test GLOBAL_OPTIONS as this method might have
    # been called before main()
    return "--is-test" in sys.argv


def getCallerLocals(decorators=0):
    '''returns the locals of the calling function.

    from http://pylab.blogspot.com/2009/02/python-accessing-caller-locals-from.html

    Arguments
    ---------
    decorators : int
        Number of contexts to go up to reach calling function
        of interest.

    Returns
    -------
    locals : dict
        Dictionary of variable defined in the context of the
        calling function.
    '''
    f = sys._getframe(2 + decorators)
    args = inspect.getargvalues(f)
    return args[3]


def getCaller(decorators=0):
    """return the name of the calling module.

    Arguments
    ---------
    decorators : int
        Number of contexts to go up to reach calling function
        of interest.

    Returns
    -------
    mod : object
        The calling module
    """

    frm = inspect.stack()[2 + decorators]
    mod = inspect.getmodule(frm[0])
    return mod

