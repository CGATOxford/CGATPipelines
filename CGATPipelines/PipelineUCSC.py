'''PipelineUCSC.py - Tasks for accessing UCSC data
==================================================

This module provides methods for accessing ucsc_ data via the public
UCSC relational database server and building UCSC track hubs.

The function :func:`connectToUCSC` establishes a connection.

Data import
-----------

The functions :func:`getRepeatsFromUCSC`, :func:`getRefseqFromUCSC` and
:func:`getCpGIslandsFromUCSC` import various data sets and save them
in standard genomic file formats.

UCSC track hubs
---------------

The functions :func:`readUCSCFile`, :func:`writeUCSCFile`,
:func:`readTrackFile` and :func:`writeTrackFile` support building a
UCSC track hub.

Reference
---------

'''

# for UCSC import
import os
import collections
import MySQLdb
import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools


import CGATPipelines.Pipeline as P


def connectToUCSC(host="genome-mysql.cse.ucsc.edu",
                  user="genome",
                  database="hg19"):
    """connect to UCSC database.

    Arguments
    ---------
    host : string
        Host to connect to
    user : string
        Username to connect with
    Database : string
        database to use

    Returns
    -------
    Database handle

    """
    dbhandle = MySQLdb.Connect(host=host,
                               user=user)

    cc = dbhandle.cursor()
    cc.execute("USE %s " % database)

    return dbhandle


def getRepeatsFromUCSC(dbhandle,
                       repclasses,
                       outfile,
                       remove_contigs_regex=None):
    '''download repeats from UCSC database and write to `outfile` in
    :term:`gff` format.

    This method downloads repeats from the repeatmasker track at
    the UCSC.

    Arguments
    ---------
    dbhandle : object
       Database handle to UCSC mysql database
    repclasses : list
       List of repeat classes to select. If empty, all repeat classes
       will be collected.
    outfile : string
       Filename of output file in :term:`gff` format.
    remove_contigs_regex : string
       If given, remove repeats on contigs matching the regular
       expression given.

    '''

    # Repeats are either stored in a single ``rmsk`` table (hg19) or in
    # individual ``rmsk`` tables (mm9) like chr1_rmsk, chr2_rmsk, ....
    # In order to do a single statement, the ucsc mysql database is
    # queried for tables that end in rmsk.
    cc = dbhandle.cursor()
    cc.execute("SHOW TABLES LIKE '%rmsk'")
    tables = [x[0] for x in cc.fetchall()]
    if len(tables) == 0:
        raise ValueError("could not find any `rmsk` tables")

    # now collect repeats
    tmpfile = P.getTempFile(".")

    for table in tables:

        cc = dbhandle.cursor()
        sql = """SELECT genoName, 'repeat', 'exon', genoStart+1, genoEnd,
        '.', strand, '.',
        CONCAT('class \\"', repClass, '\\"; family \\"',
        repFamily, '\\"; repName \\"', repName, '\\";')
        FROM %(table)s"""

        if repclasses:
            repclasses_str = ",".join(
                ["'" + x.strip() + "'" for x in repclasses])
            sql += ''' WHERE repClass in (%(repclasses_str)s) ''' % locals()

        sql = sql % locals()

        E.debug("executing sql statement: %s" % sql)
        cc.execute(sql)
        for data in cc.fetchall():
            tmpfile.write("\t".join(map(str, data)) + "\n")

    tmpfile.close()

    # sort gff and make sure that names are correct
    tmpfilename = tmpfile.name

    statement = ['''cat %(tmpfilename)s
    | %(pipeline_scriptsdir)s/gff_sort pos
    | cgat gff2gff
    --method=sanitize
    --sanitize-method=genome
    --skip-missing
    --genome-file=%(genome_dir)s/%(genome)s
    --log=%(outfile)s.log ''']

    if remove_contigs_regex:
        statement.append(
            ''' --contig-pattern="%(remove_contigs_regex)s" ''')

    statement.append('''| gzip > %(outfile)s ''')

    statement = " ".join(statement)

    P.run()

    os.unlink(tmpfilename)


def getRefSeqFromUCSC(dbhandle, outfile, remove_duplicates=False):
    '''get refseq gene set from UCSC database and save as :term:`gtf`
    formatted file.

    Matches to ``chr_random`` are ignored (as does ENSEMBL).

    Note that this approach does not work as a gene set, as refseq
    maps are not real gene builds and unalignable parts cause
    differences that are not reconcilable.

    Arguments
    ---------
    dbhandle : object
       Database handle to UCSC mysql database
    outfile : string
       Filename of output file in :term:`gtf` format. The filename
       aims to be close to the ENSEMBL gtf format.
    remove_duplicate : bool
       If True, duplicate mappings are removed.

    '''

    duplicates = set()

    if remove_duplicates:
        cc = dbhandle.cursor()
        cc.execute("""SELECT name, COUNT(*) AS c FROM refGene
        WHERE chrom NOT LIKE '%_random'
        GROUP BY name HAVING c > 1""")
        duplicates = set([x[0] for x in cc.fetchall()])
        E.info("removing %i duplicates" % len(duplicates))

    # these are forward strand coordinates
    statement = '''
    SELECT gene.name, link.geneName, link.name, gene.name2, product,
    protAcc, chrom, strand, cdsStart, cdsEnd,
    exonCount, exonStarts, exonEnds, exonFrames
    FROM refGene as gene, refLink as link
    WHERE gene.name = link.mrnaAcc
    AND chrom NOT LIKE '%_random'
    ORDER by chrom, cdsStart
    '''

    outf = IOTools.openFile(outfile, "w")

    cc = dbhandle.cursor()
    cc.execute(statement)

    SQLResult = collections.namedtuple(
        'Result',
        '''transcript_id, gene_id, gene_name, gene_id2, description,
        protein_id, contig, strand, start, end,
        nexons, starts, ends, frames''')

    counts = E.Counter()
    counts.duplicates = len(duplicates)

    for r in map(SQLResult._make, cc.fetchall()):

        if r.transcript_id in duplicates:
            continue

        starts = list(map(int, r.starts.split(",")[:-1]))
        ends = list(map(int, r.ends.split(",")[:-1]))
        frames = list(map(int, r.frames.split(",")[:-1]))

        gtf = GTF.Entry()
        gtf.contig = r.contig
        gtf.source = "protein_coding"
        gtf.strand = r.strand
        gtf.gene_id = r.gene_id
        gtf.transcript_id = r.transcript_id
        gtf.addAttribute("protein_id", r.protein_id)
        gtf.addAttribute("transcript_name", r.transcript_id)
        gtf.addAttribute("gene_name", r.gene_name)

        assert len(starts) == len(ends) == len(frames)

        if gtf.strand == "-":
            starts.reverse()
            ends.reverse()
            frames.reverse()

        counts.transcripts += 1
        i = 0
        for start, end, frame in zip(starts, ends, frames):
            gtf.feature = "exon"
            counts.exons += 1
            i += 1
            gtf.addAttribute("exon_number", i)
            # frame of utr exons is set to -1 in UCSC
            gtf.start, gtf.end, gtf.frame = start, end, "."
            outf.write("%s\n" % str(gtf))

            cds_start, cds_end = max(r.start, start), min(r.end, end)
            if cds_start >= cds_end:
                # UTR exons have no CDS
                # do not expect any in UCSC
                continue
            gtf.feature = "CDS"
            # invert the frame
            frame = (3 - frame % 3) % 3
            gtf.start, gtf.end, gtf.frame = cds_start, cds_end, frame
            outf.write("%s\n" % str(gtf))

    outf.close()

    E.info("%s" % str(counts))


def getCpGIslandsFromUCSC(dbhandle, outfile):
    '''get CpG islands from UCSC database and save as a :term:`bed`
    formatted file.

    The name column in the bed file will be set to the UCSC name.

    Arguments
    ---------
    dbhandle : object
       Database handle to UCSC mysql database
    outfile : string
       Filename of output file in :term:`bed` format.
    '''

    cc = dbhandle.cursor()
    table = "cpgIslandExt"
    sql = """SELECT chrom, chromStart, chromEnd, name
    FROM %(table)s ORDER by chrom, chromStart"""
    sql = sql % locals()

    E.debug("executing sql statement: %s" % sql)
    try:
        cc.execute(sql)
        outfile = IOTools.openFile(outfile, "w")
        for data in cc.fetchall():
            outfile.write("\t".join(map(str, data)) + "\n")
        outfile.close()
    except Exception:
        E.warn("Failed to connect to table %s. %s is empty" % (table, outfile))
        P.touch(outfile)


def readUCSCFile(infile):
    '''read data within a UCSC formatted file.

    An example for a UCSC formatted file is :file:`hub.txt` with
    the following contents::

       hub hub_name
       shortLabel hub_short_label
       longLabel hub_long_label
       genomesFile genomes_filelist
       email email_address
       descriptionUrl descriptionUrl

    Arguments
    ---------
    infile : File
        Iterator over the contents of the file.

    Returns
    -------
    data : list
        List of tuples of key/value pairs in the file.
    '''
    result = []
    for line in infile:
        if line.startswith("#"):
            continue
        if line.strip() == "":
            continue
        data = line[:-1].split()
        result.append((data[0], " ".join(data[1:])))
    return result


def writeUCSCFile(outfile, data):
    '''write a UCSC formatted text file.

    See :func:`readUCSCFile`.

    Arguments
    ---------
    outfile : File
        Output file
    data : list
        List of tuples of key/value pairs to write.
    '''
    for key, value in data:
        outfile.write(" ".join((key, value)) + "\n")


def readTrackFile(infile):
    '''read a track file.

    An example of an UCSC track file is::

      track dnaseSignal
      bigDataUrl dnaseSignal.bigWig
      shortLabel DNAse Signal
      longLabel Depth of alignments of DNAse reads
      type bigWig

      track dnaseReads
      bigDataUrl dnaseReads.bam
      shortLabel DNAse Reads
      longLabel DNAse reads mapped with MAQ
      type bam

    Arguments
    ---------
    infile : File
        Iterator over the contents of the file.

    Returns
    -------
    tracks : list
        List of tracks, each being a dictionary of keywords.
    '''

    data = readUCSCFile(infile)

    def _yielder(data):
        track, block = None, []
        for key, value in data:
            if key == "track":
                if block:
                    yield track, block
                block = []
                track = value
                continue
            block.append((key, value))
        # T.S need to yield final block too
        yield track, block

    return list(_yielder(data))


def writeTrackFile(outfile, tracks):
    '''write list of tracks to file.

    See :func:`readTrackFile`.

    Arguments
    ---------
    outfile : File
        Output file
    data : list
        List of tracks, each being a dictionary of key/value pairs.
    '''

    for track, trackdata in tracks:
        outfile.write(" ".join(("track", track)) + "\n")
        for key, value in trackdata:
            outfile.write(" ".join((key, value)) + "\n")
        outfile.write("\n")
