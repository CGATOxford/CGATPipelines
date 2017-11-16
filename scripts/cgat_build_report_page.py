'''
cgat_build_report_page.py - build a report summary across projects
=======================================================================


Purpose
-------

This script scans all of :file:`/ifs/projects/sftp` for
:file:`index.html` files and outputs an html formatted summary table
into :file:`/ifs/projects/overview`.

Usage
-----

To build the report, execute::

   python cgat_build_report_page.py

Type::

   python cgat_build_report_page.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import subprocess

import CGAT.Experiment as E
import CGAT.IOTools as IOTools


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id",
                            usage=globals()["__doc__"])

    parser.add_option("-p", "--path", dest="path", type="string",
                      help="path to scan for files [%default]")

    parser.add_option("-d", "--destination", dest="destination", type="string",
                      help="path to deposit files into [%defaul]")

    parser.set_defaults(path='/ifs/projects/sftp',
                        url='http://www.cgat.org/downloads/',
                        dest='/ifs/projects/overview')

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    statement = "find %s -name 'index.html'" % options.path

    process = subprocess.Popen(statement,
                               shell=True,
                               stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)

    stdout, stderr = process.communicate()

    files = stdout.split('\n')
    files.sort()

    outfile = IOTools.openFile(os.path.join(options.dest, "index.html"), "w")

    outfile.write('''
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
  "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">

<html>
<head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
    <title>CGAT project reports</title>
    <link rel="stylesheet" href="cgat.css" type="text/css" />
    <link rel="stylesheet" href="pygments.css" type="text/css" />
    <link rel="shortcut icon" href="http://cgatwiki.anat.ox.ac.uk/favicon.ico">
    <script type="text/javascript" src="sorttable.js"></script>
</head>

  <body>
    <div class="related">
      <h3>Navigation</h3>
      <ul>
        <li><a href="index.html">CGAT Projects Overview</a> &raquo;</li>
      </ul>
    </div>

    <div class="document">
      <div class="documentwrapper">
        <div class="bodywrapper">
          <div class="body">
 <div class="section" id="cgat-pipelines">
<H1>CGAT exported project pages</H1>

<p> 
This page is for internal use only. Do not distribute outside of
CGAT and do not make this page available on the world wide web.
</p>

<table class="sortable">\n''')

    outfile.write(
        '''<tr><th>Project</th><th>Report</th><th>Title</th></tr>\n''')

    for f in files:
        if f == '':
            continue

        proj = re.search('(proj\d+)', f).groups()[0]
        relpath = re.sub('.*proj\d+/', '', f)
        report = re.sub('^[^/]*/', '', os.path.dirname(relpath))

        lines = IOTools.openFile(f).readlines()
        titles = [x for x in lines if "<title>" in x]
        if titles:
            title = re.search("<title>(.*)</title>", titles[0]).groups()[0]
        else:
            title = "NA"

        if title.endswith("documentation"):
            title = title[:-len("documentation")]

        url = os.path.join(options.url, relpath)
        outfile.write(
            '<tr><td>%(proj)s</td><td><a HREF="%(url)s">%(report)s</td><td>%(title)s</td></tr>\n' % locals())

    outfile.write('''
</table>

</div>
</div>


          </div>
        </div>
      </div>
      <div class="sphinxsidebar">
        <div class="sphinxsidebarwrapper">
            <p class="logo"><a href="contents.html">
              <img class="logo" src="cgat_logo.png" alt="Logo"/>
            </a></p>





</body>
</html>\n''')

    outfile.close()

    E.info('created output file %s' % outfile.name)
    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
