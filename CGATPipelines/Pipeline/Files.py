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
"""Files.py - Working with files in ruffus pipelines
====================================================

Reference
---------

"""
import os
import tempfile

import CGAT.IOTools as IOTools

# Set from Pipeline.py
PARAMS = {}


def getTempFile(dir=None, shared=False):
    '''get a temporary file.

    The file is created and the caller needs to close and delete
    the temporary file once it is not used any more.

    Arguments
    ---------
    dir : string
        Directory of the temporary file and if not given is set to the
        default temporary location in the global configuration dictionary.
    shared : bool
        If set, the tempory file will be in a shared temporary
        location (given by the global configuration directory).

    Returns
    -------
    file : File
        A file object of the temporary file.
    '''
    if dir is None:
        if shared:
            dir = PARAMS['shared_tmpdir']
        else:
            dir = PARAMS['tmpdir']

    return tempfile.NamedTemporaryFile(dir=dir, delete=False, prefix="ctmp")


def getTempFilename(dir=None, shared=False):
    '''return a temporary filename.

    The file is created and the caller needs to delete the temporary
    file once it is not used any more.

    Arguments
    ---------
    dir : string
        Directory of the temporary file and if not given is set to the
        default temporary location in the global configuration dictionary.
    shared : bool
        If set, the tempory file will be in a shared temporary
        location.

    Returns
    -------
    filename : string
        Absolute pathname of temporary file.

    '''
    tmpfile = getTempFile(dir=dir, shared=shared)
    tmpfile.close()
    return tmpfile.name


def getTempDir(dir=None, shared=False):
    '''get a temporary directory.

    The directory is created and the caller needs to delete the temporary
    directory once it is not used any more.

    Arguments
    ---------
    dir : string
        Directory of the temporary directory and if not given is set to the
        default temporary location in the global configuration dictionary.
    shared : bool
        If set, the tempory directory will be in a shared temporary
        location.

    Returns
    -------
    filename : string
        Absolute pathname of temporary file.

    '''
    if dir is None:
        if shared:
            dir = PARAMS['shared_tmpdir']
        else:
            dir = PARAMS['tmpdir']

    return tempfile.mkdtemp(dir=dir, prefix="ctmp")


def checkExecutables(filenames):
    """check for the presence/absence of executables"""

    missing = []

    for filename in filenames:
        if not IOTools.which(filename):
            missing.append(filename)

    if missing:
        raise ValueError("missing executables: %s" % ",".join(missing))


def checkScripts(filenames):
    """check for the presence/absence of scripts"""
    missing = []
    for filename in filenames:
        if not os.path.exists(filename):
            missing.append(filename)

    if missing:
        raise ValueError("missing scripts: %s" % ",".join(missing))


