'''
CGAT.py - CGAT project specific functions
=========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

The :mod:`CGAT` module contains various utility functions
for working on CGAT projects.

API
----
'''

import os
import collections
import CGAT.Experiment as E
import CGATPipelines.Pipeline as P
import CGAT.IOTools as IOTools
import CGATPipelines.PipelineUCSC as PipelineUCSC


PROJECT_ROOT = '/ifs/projects'

# gets set by importing script
PARAMS = {}


def publish_tracks(export_files,
                   prefix="",
                   project_id=None,
                   project_name=None):
    '''publish a UCSC Track Hub.

    This method takes a dictionary of file types associated
    with files. For each file, a link will be created in
    the upload directory. The track will be stored under
    a project name, which will be derived from the location
    of the working directory.

    Information about the genome, the upload directory, etc. will be
    taken from the global configuration dictionary.

    For example, calling the following code in a pipeline executed
    in .../proj013/mapping::
        export_files = {
            "bamfiles": glob.glob("*/*.bam") + \
                         glob.glob("*/*.bam.bai"),
            "bigwigfiles": glob.glob("*/*.bw"),
        }
        publish_tracks(export_files)

    Will create a hub file at
    :file:`<uploaddir>/OBFUSID/mapping/ucsc.hub`, where
    OBFUSID is the obfuscated directory entry in the CGAT
    download directory for a particular project.

    Arguments
    =========
    export_files : dict
        Dictionary mapping filetypes to files.
    prefix : string
        will be added to each track.
    project_id : string
        The project identifier. If not given, it will be taken from
        the path of the project directory.
    project_name : string
        The project name, typically the project number. If not given,
        it will be taken from the current directory.

    '''

    if not prefix:
        prefix = PARAMS.get("report_prefix", "")

    web_dir = PARAMS["web_dir"]
    if project_id is None:
        project_id = P.getProjectId()
    if project_name is None:
        project_name = P.getProjectName()

    src_export = os.path.abspath("export")
    dest_report = prefix + "report"
    dest_export = prefix + "export"

    hubdir = os.path.join(PARAMS["web_dir"], "ucsc")

    if not os.path.exists(hubdir):
        E.info("creating %s" % hubdir)
        os.mkdir(hubdir)

    # write the UCSC hub file
    hubfile = os.path.join(hubdir, "hub.txt")
    genomesfile = os.path.join(hubdir, "genomes.txt")
    trackdir = os.path.join(hubdir, PARAMS["genome"])
    trackfile = os.path.join(hubdir, PARAMS["genome"], "trackDb.txt")
    trackrelpath = os.path.join(PARAMS["genome"], "trackDb.txt")

    if os.path.exists(hubfile):
        with IOTools.openFile(hubfile) as infile:
            hubdata = PipelineUCSC.readUCSCFile(infile)
    else:
        hubdata = [('hub', "CGAT-" + project_name),
                   ('shortLabel', "CGAT-" + project_name),
                   ('longLabel', "Data for CGAT project %s" % project_name),
                   ('genomesFile', "genomes.txt"),
                   ('email', 'andreas.heger@gmail.com')]

    E.info("writing to %s" % hubfile)
    with IOTools.openFile(hubfile, "w") as outfile:
        PipelineUCSC.writeUCSCFile(outfile, hubdata)

    # create the genomes.txt file - append to it if necessary.
    if os.path.exists(genomesfile):
        with IOTools.openFile(genomesfile) as infile:
            genomes = PipelineUCSC.readUCSCFile(infile)
    else:
        genomes = []

    if ("genome", PARAMS["genome"]) not in genomes:
        genomes.append(("genome", PARAMS["genome"]))
        genomes.append(("trackDb", trackrelpath))

    E.info("writing to %s" % genomesfile)
    with IOTools.openFile(genomesfile, "w") as outfile:
        PipelineUCSC.writeUCSCFile(outfile, genomes)

    # create the track data
    if not os.path.exists(trackdir):
        os.mkdir(trackdir)

    if os.path.exists(trackfile):
        E.debug('reading existing tracks from %s' % trackfile)
        with IOTools.openFile(trackfile) as infile:
            tracks = PipelineUCSC.readTrackFile(infile)
    else:
        tracks = []

    tracks = collections.OrderedDict(tracks)

    def getName(name):
        if name.endswith(".bam"):
            return "bam", name
        elif name.endswith(".bw") or name.endswith(".bigwig"):
            return "bigWig", name
        else:
            return None, None

    for targetdir, filenames in export_files.items():
        for src in filenames:
            dest = os.path.join(trackdir, prefix + os.path.basename(src))
            dest = os.path.abspath(dest)
            # create a symlink
            if not os.path.exists(dest):
                try:
                    os.symlink(os.path.abspath(src), dest)
                except OSError, msg:
                    E.warn("could not create symlink from %s to %s: %s" %
                           (os.path.abspath(src), dest, msg))
            ucsctype, trackname = getName(os.path.basename(dest))
            # ignore invalid types and other files (.bai files, ...)
            if ucsctype is None:
                continue
            tracks[trackname] = (("bigDataUrl", os.path.basename(dest)),
                                 ("shortLabel", trackname),
                                 ("longLabel", trackname),
                                 ("type", ucsctype))

    E.info("writing to %s" % trackfile)
    with IOTools.openFile(trackfile, "w") as outfile:
        PipelineUCSC.writeTrackFile(outfile, list(tracks.iteritems()))

    E.info(
        "data hub has been created at http://www.cgat.org/downloads/%(project_id)s/ucsc/hub.txt" % locals())
