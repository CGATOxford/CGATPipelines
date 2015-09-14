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
'''Local.py - CGAT project specific functions
=============================================

The :mod:`Local` module contains various utility functions for working
on CGAT projects and are very specific to the CGAT directory layout.

.. note::

   Methods in this module need to made to work with arbitrary project
   layouts.

CGAT project layout
-------------------

The method :func:`isCGAT` checks if the code is executed within the
CGAT systems. The functions :func:`getProjectDirectories`,
:func:`getPipelineName`, :func:`getProjectId`, :func:`getProjectName`
provide information about the pipeline executed and the project context.

Report publishing
-----------------

Once built, a report can be published by copying it to the publicly
visible directories on the CGAT systems. At the same time, references
to files on CGAT systems need to be replaced with links through the
public interfacet. The functions :func:`getPublishDestinations` and
:func:`publish_report` implement this functionality.

Reference
---------

'''
import os
import re
import shutil
import inspect

from CGAT import Experiment as E

PROJECT_ROOT = '/ifs/projects'

# Variables PARAMS and CONFIG will be set by Pipeline.py
# on import.
PARAMS = None
CONFIG = None


def isCGAT(curdir=None):
    '''return True if this is a CGAT project.

    This method works by checking if the current working directory
    is part of :var:`PROJECT_ROOT`.
    '''
    if curdir is None:
        curdir = os.path.abspath(os.getcwd())

    return curdir.startswith(PROJECT_ROOT)


def getProjectDirectories(sections=None):
    '''return directories relevant to this project.

    The entries of the dictionary are:

    webdir
       Directory for publishing information (without password access).
    exportdir
       Directory for storing files to be exported alongside
       the report.
    notebookdir
       Directory where project notebooks are located.

    Arguments
    ---------
    sections : list
        If given, only the named sections are returned.

    Returns
    -------
    directories : dict

    Raises
    ------
    ValueError
       If any of the directories does not exist

    '''

    if not isCGAT():
        raise ValueError(
            "getProjectDirectories called for a non-CGAT project")

    project_name = getProjectName()

    result = {
        'webdir': os.path.join(
            PROJECT_ROOT, PARAMS["web_dir"]),
        'exportdir': os.path.join(
            PARAMS["exportdir"]),
        'notebookdir': os.path.join(
            PROJECT_ROOT, project_name, "notebooks")
    }

    if sections:
        result = dict([(x, y) for x, y in result.items()
                       if x in sections])

    for x, y in result.items():
        if not os.path.exists(y):
            raise ValueError(
                "directory %s for %s does not exist" % (y, x))

    return result


def getPipelineName():
    '''return the name of the pipeline.

    The name of the pipeline is deduced by the name of the top-level
    python script. The pipeline name is the name of the script
    without any path information and the ``.py`` suffix.

    Returns
    -------
    string

    '''
    # use the file attribute of the caller
    for x in inspect.stack():
        if x[0].f_globals["__name__"] == "__main__":
            return os.path.basename(x[0].f_globals['__file__'])[:-3]


def getProjectId():
    '''get the (obfuscated) project id based on the current working
    directory.

    The project is located by finding the ``web_dir`` configuration
    variable and working backwards from that. ``web_dir`` should be
    link to the web directory in the project directory which then
    links to the web directory in the sftp directory which then links
    to the obfuscated directory::

        pipeline:web_dir
        -> /ifs/projects/.../web
        -> /ifs/sftp/.../web
        -> /ifs/sftp/.../aoeuCATAa (obfuscated directory)

    Returns
    =======
    string

    '''
    # return an id that has been explicitely set
    if "report_project_url" in PARAMS:
        return PARAMS["report_project_url"]

    curdir = os.path.abspath(os.getcwd())
    if not isCGAT(curdir):
        raise ValueError(
            "method getProjectId not called within %s" % PROJECT_ROOT)

    webdir = PARAMS['web_dir']
    if not os.path.islink(webdir):
        raise ValueError(
            "unknown configuration: webdir '%s' is not a link" % webdir)
    target = os.readlink(webdir)
    if not os.path.islink(target):
        raise ValueError(
            "unknown configuration: target '%s' is not a link" % target)

    return os.path.basename(os.readlink(target))


def getProjectName():
    '''cgat specific method: get the name of the project
    based on the current working directory.

    If called outside the Project hierarchy, the project name
    will be set to the name of the current directory.
    '''

    curdir = os.path.abspath(os.getcwd())
    if isCGAT(curdir):
        prefixes = len(PROJECT_ROOT.split("/"))
        return curdir.split("/")[prefixes]
    else:
        return os.path.basename(curdir)


def getPublishDestinations(prefix="", suffix=None):

    if not prefix:
        prefix = PARAMS.get("report_prefix", "default")

    if prefix == "default":
        prefix = getPipelineName() + "_"

    if not suffix:
        suffix = PARAMS.get("report_suffix", "")

    dest_report = prefix + "report"
    dest_export = prefix + "export"

    if suffix is not None:
        dest_report += suffix
        dest_export += suffix

    return dest_report, dest_export


def publish_report(prefix="",
                   patterns=[],
                   project_id=None,
                   prefix_project="/ifs/projects",
                   export_files=None,
                   suffix=None,
                   subdirs=False,
                   ):
    '''publish report into web directory.

    Links export directory into web directory.

    Copies html pages and fudges links to the pages in the
    export directory.

    If *prefix* is given, the directories will start with prefix,
    otherwise, it is looked up from the option ``report_prefix``.
    If report_prefix is "default", the prefix will be derived
    from the pipeline name. For example, pipeline_intervals will
    we copied to ``pipeline_intervals_report``.

    *patterns* is an optional list of two-element tuples (<pattern>,
    replacement_string).  Each substitutions will be applied on each
    file ending in .html.

    If *project_id* is not given, it will be looked up. This requires
    that this method is called within a subdirectory of PROJECT_ROOT.

    *export_files* is a dictionary of files to be exported. The key
    of the dictionary denotes the targetdirectory within the web
    directory. The values in the dictionary are the files to be
    linked to in the direcotry. For example::

        exportfiles = {
            "bamfiles" : glob.glob( "*/*.bam" ) + glob.glob( "*/*.bam.bai" ),
            "bigwigfiles" : glob.glob( "*/*.bw" ),
            }

    .. note::
       This function is CGAT specific.

    '''

    dest_report, dest_export = getPublishDestinations(prefix, suffix)

    web_dir = PARAMS["web_dir"]

    if project_id is None:
        project_id = getProjectId()

    src_export = os.path.abspath("export")
    curdir = os.path.abspath(os.getcwd())

    # substitute links to export and report
    base_url = "http://www.cgat.org/downloads/%s" % project_id
    _patterns = [
        # redirect export directory
        (re.compile(src_export),
         "%(base_url)s/%(dest_export)s" % locals()),
        # redirect report directory
        (re.compile(curdir),
         "%(base_url)s/%(dest_report)s" % locals()),
        (re.compile('(%s)/_static' %
                    curdir),
         "%(base_url)s/%(dest_report)s/_static" % locals())]

    _patterns.extend(patterns)

    # add intersphinx mapping - this requires that the name
    # for the interpshinx redirection (key) corresponds to the
    # export location with an appended "_report".
    if CONFIG.has_section("intersphinx"):
        for key, value in CONFIG.items("intersphinx"):
            _patterns.append((
                re.compile(os.path.abspath(value)),
                "%(base_url)s/%(key)s_report" % locals()))
            # check if the target exists in download location
            intersphinx_target = os.path.join(
                web_dir, key + "_report", "objects.inv")
            if not os.path.exists(intersphinx_target):
                E.warn("intersphinx mapping for '%s' does not exist at %s" %
                       (key, intersphinx_target))

    def _link(src, dest):
        '''create links.

        Only link to existing targets.
        '''
        if os.path.exists(dest):
            os.remove(dest)

        if not os.path.exists(src):
            E.warn("%s does not exist - skipped" % src)
            return

        # IMS: check if base path of dest exists. This allows for
        # prefix to be a nested path structure e.g. project_id/
        if not os.path.exists(os.path.dirname(os.path.abspath(dest))):
            E.info('creating directory %s' %
                   os.path.dirname(os.path.abspath(dest)))
            os.mkdir(os.path.dirname(os.path.abspath(dest)))

        os.symlink(os.path.abspath(src), dest)

    def _copy(src, dest):
        if os.path.exists(dest):
            shutil.rmtree(dest)
        if not os.path.exists(src):
            E.warn("%s does not exist - skipped" % src)
            return
        shutil.copytree(os.path.abspath(src), dest)

    # publish export dir via symlinking
    E.info("linking export directory in %s" % dest_export)
    _link(src_export,
          os.path.abspath(os.path.join(web_dir, dest_export)))

    # publish web pages by copying
    E.info("publishing web pages in %s" %
           os.path.abspath(os.path.join(web_dir, dest_report)))
    _copy(os.path.abspath("report/html"),
          os.path.abspath(os.path.join(web_dir, dest_report)))

    for root, dirs, files in os.walk(os.path.join(web_dir, dest_report)):
        for f in files:
            fn = os.path.join(root, f)
            if fn.endswith(".html"):
                with open(fn) as inf:
                    data = inf.read()
                for rx, repl in _patterns:
                    data = rx.sub(repl, data)
                outf = open(fn, "w")
                outf.write(data)
                outf.close()

    if export_files:
        bigwigs, bams, beds = [], [], []

        for targetdir, filenames in export_files.items():

            targetdir = os.path.join(web_dir, targetdir)
            if not os.path.exists(targetdir):
                os.makedirs(targetdir)

            for src in filenames:
                dest = os.path.join(targetdir, os.path.basename(src))
                if dest.endswith(".bam"):
                    bams.append((targetdir, dest))
                elif dest.endswith(".bw"):
                    bigwigs.append((targetdir, dest))
                elif dest.endswith(".bed.gz"):
                    beds.append((targetdir, dest))
                dest = os.path.abspath(dest)
                if not os.path.exists(dest):
                    try:
                        os.symlink(os.path.abspath(src), dest)
                    except OSError, msg:
                        E.warn("could not create symlink from %s to %s: %s" %
                               (os.path.abspath(src), dest, msg))

        # output ucsc links
        with open("urls.txt", "w") as outfile:
            for targetdir, fn in bams:
                filename = os.path.basename(fn)
                track = filename[:-len(".bam")]
                outfile.write(
                    """track type=bam name="%(track)s" bigDataUrl=http://www.cgat.org/downloads/%(project_id)s/%(targetdir)s/%(filename)s\n""" % locals())

            for targetdir, fn in bigwigs:
                filename = os.path.basename(fn)
                track = filename[:-len(".bw")]
                outfile.write(
                    """track type=bigWig name="%(track)s" bigDataUrl=http://www.cgat.org/downloads/%(project_id)s/%(targetdir)s/%(filename)s\n""" % locals())

            for targetdir, fn in beds:
                filename = os.path.basename(fn)
                track = filename[:-len(".bed.gz")]
                outfile.write(
                    """http://www.cgat.org/downloads/%(project_id)s/%(targetdir)s/%(filename)s\n""" % locals())

        E.info("UCSC urls are in urls.txt")

    E.info(
        "report has been published at http://www.cgat.org/downloads/%(project_id)s/%(dest_report)s" % locals())


