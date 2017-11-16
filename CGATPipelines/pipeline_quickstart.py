'''pipeline_quickstart.py - setup a new pipeline
=============================================


Purpose
-------

This script creates a new pipeline according to the CGAT
pipeline layout. This is useful when starting a new pipeline
from scratch.

Usage
-----

To start a new project, use  :file:`pipeline_quickstart.py`::

   python <srcdir>pipeline_quickstart.py --set-name=chipseq

This will create a new directory called ``chipseq`` in the current directory
with the following layout::

  |-- [         55]  work
  |   |-- [         49]  conf.py -> ../src/pipeline_chipseq/conf.py
  |   `-- [         54]  pipeline.ini -> ../src/pipeline_chipseq/pipeline.ini
  `-- [        102]  src
      |-- [         55]  pipeline_chipseq
      |   |-- [      13142]  conf.py
      |   `-- [       1232]  pipeline.ini
      |-- [       6003]  pipeline_chipseq.py
      `-- [         58]  pipeline_docs
          |-- [        169]  pipeline_chipseq
          |   |-- [         24]  __init__.py
          |   |-- [         31]  _templates
          |   |   `-- [      21703]  cgat_logo.png
          |   |-- [       1052]  contents.rst
          |   |-- [         56]  pipeline
          |   |   |-- [        405]  Dummy.rst
          |   |   `-- [         45]  Methods.rst
          |   |-- [         79]  pipeline.rst
          |   `-- [         35]  trackers
          |       `-- [        301]  TemplateReport.py
          `-- [         22]  themes
              `-- [         81]  cgat
                  |-- [        319]  layout.html
                  |-- [         58]  static
                  |   |-- [       4362]  cgat.css_t
                  |   `-- [      16973]  sorttable.js
                  `-- [         69]  theme.conf


The layout has the following components::

work
   Directory for running the pipeline. Links to configuration files in
   the :file:`src` directory. This directory exists to separate code
   from data and results.
src
   Directory for pipeline code. This directory should be put under
   version control and backed-up.
src/pipeline_chipseq.py
   The main pipeline script
src/pipeline_chipseq
   Directory for the pipeline configuration files
src/pipeline_docs/pipeline_chipseq
   Directory for the pipeline report
src/pipeline_docs/themes
   The CGAT report theme.

Type::

   python pipeline_quickstart.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import re
import os
import shutil
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGATPipelines.Pipeline as P


def main(argv=sys.argv):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-d", "--dest", dest="destination", type="string",
                      help="destination directory.")

    parser.add_option(
        "-n", "--name", "--set-name", dest="name", type="string",
        help="name of this pipeline. 'pipeline_' will be prefixed.")

    parser.add_option(
        "-f", "--force-output", dest="force", action="store_true",
        help="overwrite existing files.")

    parser.add_option(
        "-t", "--pipeline-type", dest="pipeline_type", type="choice",
        choices=("full", "minimal"),
        help="type of pipeline to output. "
        "full=a complete pipeline for the CGAT environment "
        "minimum=minimum pipeline "
        "[%default]")

    parser.set_defaults(
        destination=".",
        name=None,
        force=False,
        pipeline_type="full",
    )

    (options, args) = E.Start(parser)

    if not options.name:
        raise ValueError("please provide a pipeline name")

    reportdir = os.path.abspath("src/pipeline_docs/pipeline_%s" % options.name)
    confdir = os.path.abspath("src/pipeline_%s" % (options.name))

    destination_dir = options.destination

    # create directories
    for d in ("", "src", "work",
              "src/pipeline_docs",
              "src/pipeline_%s" % options.name,
              reportdir,
              "%s/_templates" % reportdir,
              "%s/pipeline" % reportdir,
              "%s/trackers" % reportdir):

        dd = os.path.join(destination_dir, d)
        if not os.path.exists(dd):
            os.makedirs(dd)

    # copy files
    # replaces all instances of template with options.name within
    # filenames and inside files.
    rx_file = re.compile("template")
    rx_type = re.compile("_%s" % options.pipeline_type)
    rx_template = re.compile("@template@")
    rx_reportdir = re.compile("@reportdir@")

    srcdir = P.CGATPIPELINES_PIPELINE_DIR

    def copy(src, dst, name):

        # remove "template" and the pipeline type from file/directory
        # names.
        fn_dest = os.path.join(
            destination_dir,
            dst,
            rx_type.sub("", rx_file.sub(name, src)))

        fn_src = os.path.join(srcdir,
                              "pipeline_template_data", src)

        E.debug("fn_src=%s, fn_dest=%s, src=%s, dest=%s" %
                (fn_src, fn_dest, src, dst))

        if os.path.exists(fn_dest) and not options.force:
            raise OSError(
                "file %s already exists - not overwriting." % fn_dest)

        if fn_src.endswith(".png"):
            shutil.copyfile(fn_src, fn_dest)
        else:
            with IOTools.openFile(fn_dest, "w") as outfile:
                with IOTools.openFile(fn_src) as infile:
                    for line in infile:
                        outfile.write(rx_reportdir.sub(reportdir,
                                                       rx_template.sub(name, line)))

    def copytree(src, dst, name):

        fn_dest = os.path.join(destination_dir, dst, rx_file.sub(name, src))
        fn_src = os.path.join(srcdir, "pipeline_template_data", src)

        if os.path.exists(fn_dest) and not options.force:
            raise OSError(
                "file %s already exists - not overwriting." % fn_dest)

        shutil.copytree(fn_src, fn_dest)

    for f in ("conf.py",
              "pipeline.ini"):
        copy(f, 'src/pipeline_%s' % options.name, name=options.name)

    # copy the script
    copy("pipeline_template_%s.py" % options.pipeline_type, 'src',
         name=options.name)

    # create links
    for src, dest in (("conf.py", "conf.py"),
                      ("pipeline.ini", "pipeline.ini")):
        d = os.path.join("work", dest)
        if os.path.exists(d) and options.force:
            os.unlink(d)
        os.symlink(os.path.join(confdir, src), d)

    for f in ("cgat_logo.png",):
        copy(f, "%s/_templates" % reportdir,
             name=options.name)

    for f in ("themes",):
        copytree(f, "src/pipeline_docs",
                 name=options.name)

    for f in ("contents.rst",
              "pipeline.rst",
              "__init__.py"):
        copy(f, reportdir,
             name=options.name)

    for f in ("Dummy.rst",
              "Methods.rst"):
        copy(f, "%s/pipeline" % reportdir,
             name=options.name)

    for f in ("TemplateReport.py", ):
        copy(f, "%s/trackers" % reportdir,
             name=options.name)

    absdest = os.path.abspath(destination_dir)

    name = options.name

    print("""
Welcome to your new %(name)s CGAT pipeline.

All files have been successfully copied to `%(destination_dir)s`. In
order to start the pipeline, go to `%(destination_dir)s/work`

   cd %(destination_dir)s/work

You can start the pipeline by typing:

   cgatflow %(name)s -v 5 -p 5 make full

The source code for the pipeline is in %(destination_dir)s/src.

""" % locals())

    E.Stop()


if __name__ == "__main__":
    sys.exit(main())
