from ruffus import *
from ruffus.combinatorics import *
import os
import sqlite3
import CGAT.Experiment as E
import CGATPipelines.Pipeline as P
import sys


# load options from the config file
PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])


@follows(mkdir("genome.dir"))
@originate("genome.dir/%s.fa" % (PARAMS['genome_name']))
def downloadGenomeFasta(outfile):
    genome_path = PARAMS['genome_path']
    filename = genome_path.split('/')[-1]
    filestem = PARAMS['genome_name']
    statement = '''wget %(genome_path)s && cgat index_fasta 
               %(filestem)s %(filename)s''' % locals()
    P.run()

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
