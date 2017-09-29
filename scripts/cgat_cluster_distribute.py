'''cgat_cluster_distribute.py - distribute files to cluster nodes
=================================================================


Purpose
-------

This script distributes one or more files to the local
scratch directory on all nodes in the cluster.

The script uses ``rsync`` to only copy files that are newer.

Usage
-----

For example::

   python cgat_cluster_distribute.py -collection=blast /local/nrdb/uniref100.*.{pin,psd,psi,phr,psq}

The above command mirrors uniprot blast indexed databases in the
directory :file:`/local/nrdb` into the directory
:file:`/scratch/blast` on the nodes. The files will be put into a
subdirectory called ``blast`

Type::

   python cgat_cluster_distribute.py --help

for command line help.

.. note::

   The remote directories need to be cleaned up manually.

.. todo::

   Currently the list of nodes is hard-coded and files are copied to
   all nodes. This is potentially wasteful if jobs will only be
   executed on a few nodes.

Command line options
--------------------

'''

import os
import sys

import CGAT.Experiment as E


def getNodes(nodes=None):
    '''hack - allow ranges, ...'''
    if nodes is None or len(nodes) == 0:
        return ["cgat%03i" % x for x in list(range(1, 15)) + list(range(101, 117))] + \
            ["cgat150", "cgatsmp1", "cgatgpu1",
             "andromeda", "gandalf", "saruman"]
    return nodes


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option(
        "-s", "--scratch-dir", dest="scratchdir", type="string",
        help="the scratch directory on the nodes [default=%default].")

    parser.add_option(
        "-c", "--collection", dest="collection", type="string",
        help="files will be put into collection. This is a directory that "
        "will be created just below the scratch directory [default=%default].")

    parser.set_defaults(
        scratchdir="/scratch",
        collection="",
        nodes=[],
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if len(args) == 0:
        raise ValueError(
            "please specify a collection of files/directories "
            "that should be mirrored.")

    targetdir = os.path.join(options.scratchdir, options.collection)

    nodes = getNodes(options.nodes)

    E.info("copying to %s on nodes %s" % (targetdir, ",".join(nodes)))

    ninput, noutput, nskipped = 0, 0, 0

    filenames = " ".join(args)

    for node in nodes:
        E.info("copying to node %s" % node)
        ninput += 1
        statement = '''
               ssh %(node)s mkdir %(targetdir)s >& /dev/null;
               rsync --progress -az %(filenames)s %(node)s:%(targetdir)s
        ''' % locals()
        E.run(statement)
        noutput += 1

    E.info("ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput, nskipped))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
