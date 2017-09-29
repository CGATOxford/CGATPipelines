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


def getTempFile(dir=None, shared=False, suffix="", mode="w+", encoding="utf-8"):
    '''get a temporary file.

    The file is created and the caller needs to close and delete the
    temporary file once it is not used any more. By default, the file
    is opened as a text file (mode ``w+``) with encoding ``utf-8``
    instead of the default mode ``w+b`` used in
    :class:`tempfile.NamedTemporaryFile`

    Arguments
    ---------
    dir : string
        Directory of the temporary file and if not given is set to the
        default temporary location in the global configuration dictionary.
    shared : bool
        If set, the tempory file will be in a shared temporary
        location (given by the global configuration directory).
    suffix : string
        Filename suffix

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

    return tempfile.NamedTemporaryFile(dir=dir,
                                       delete=False,
                                       prefix="ctmp",
                                       mode=mode,
                                       encoding=encoding,
                                       suffix=suffix)


def getTempFilename(dir=None, shared=False, suffix=""):
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
    suffix : string
        Filename suffix

    Returns
    -------
    filename : string
        Absolute pathname of temporary file.

    '''
    tmpfile = getTempFile(dir=dir, shared=shared, suffix=suffix)
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
