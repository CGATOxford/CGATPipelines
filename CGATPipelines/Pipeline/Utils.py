"""Utils.py - Utilities for ruffus pipelines
============================================

Reference
---------

"""
import inspect
import sys


def isTest():
    """return True if the pipeline is run in a "testing" mode.

    This method checks if ``-is-test`` has been given as a
    command line option.
    """
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


def add_doc(value, replace=False):
    """add doc string of value to function that is decorated.

    The original doc-string is added as the first paragraph(s)
    inside the new doc-string.

    Parameter
    ---------

    replace : bool
       If True, replace documentation rather than appending
    """
    def _doc(func):
        if func.__doc__:
            lines = value.__doc__.split("\n")
            for x, line in enumerate(lines):
                if line.strip() == "":
                    break
            # insert appropriate indentiation
            # currently hard-coded, can be derived
            # from doc string?
            if not replace:
                lines.insert(x+1, " " * 4 +
                             func.__doc__)
                func.__doc__ = "\n".join(lines)
            else:
                func.__doc__ = value.__doc__
        else:
            func.__doc__ = value.__doc__
        return func
    return _doc
