"""
PipelineGtfsubset.py - Tasks for GTF subsetting
===============================================

Reference
---------

"""

import CGAT.Experiment as E
import os
import MySQLdb
import pysam
import CGAT.IOTools as IOTools
import CGAT.GTF as GTF
import CGATPipelines.Pipeline as P


class SubsetGTF():

    def __init__(self, infile, *args, **kwargs):

        self.gtf = GTF.iterator(IOTools.openFile(infile, "r"))

    def makeLineDict(self, line):
        D = line.asDict()
        D['chrom'] = line.contig
        D['source'] = line.source
        D['feature'] = line.feature
        D['start'] = line.start
        D['end'] = line.end
        D['score'] = line.score
        D['strand'] = line.strand
        D['frame'] = line.frame

        return D

    def filterGTF(self, outfile, filteroption, filteritem, operators):
        '''

        '''

        with IOTools.openFile(outfile, "w") as outf:
            for line in self.gtf:
                D = self.makeLineDict(line)
                if len(filteritem) == 1:
                    if D[filteroption] == filteritem[0]:
                        outf.write("%s\n" % str(line))

                elif len(filteritem) == 2:
                    filter1, filter2 = filteritem
                    filteroption1, filteroption2 = filteroption
                    if operators == "and":
                        if D[filteroption1] == filter1 and \
                           D[filteroption2] == filter2:
                            outf.write("%s\n" % str(line))
                    elif operators == "and not":
                        if D[filteroption1] == filter1 and not\
                           D[filteroption2] == filter2:
                            outf.write("%s\n" % str(line))

                else:
                    pass

        outf.close()


class SubsetGFF3():

    def __init__(self, infile, *args, **kwargs):
        self.gff = pysam.tabix_iterator(IOTools.openFile(infile),
                                        parser=pysam.asGFF3())

    def makeLineDict(self, line):
        D = line.attribute_string2dict(line.attributes)
        D['chrom'] = line.contig
        D['source'] = line.source
        D['feature'] = line.feature
        D['start'] = line.start
        D['end'] = line.end
        D['score'] = line.score
        D['strand'] = line.strand

        return D

    def filterGFF3(self, outfile, filteroption, filteritem):

        with IOTools.openFile(outfile, "w") as outf:
            for line in self.gff:
                D = self.makeLineDict(line)
                if D[filteroption] == filteritem[0]:
                    outf.write("%s\n" % str(line))
        outf.close()


def connectToUCSC(host="genome-mysql.cse.ucsc.edu",
                  user="genome",
                  database=None):
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


def getRepeatDataFromUCSC(dbhandle,
                          repclasses,
                          outfile,
                          remove_contigs_regex=None):
    '''download data from UCSC database and write to `outfile` in
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


def buildGenomicContext(infiles, outfile, distance=10):

    '''build a :term:`bed` formatted file with genomic context.
    The output is a bed formatted file, annotating genomic segments
    according to whether they are any of the ENSEMBL annotations.
    The function also adds the RNA and repeats annotations from the UCSC.
    The annotations can be partially or fully overlapping.
    The annotations can be partially or fully overlapping. Adjacent
    features (less than 10 bp apart) of the same type are merged.
    Arguments
    ---------
    infiles : list
       A list of input files to generate annotations from. The contents are
       1. ``repeats``, a :term:`gff` formatted file with repeat annotations
       2. ``rna``, a :term:`gff` formatted file with small, repetetive
          RNA annotations
       3. ``annotations``, a :term:`gtf` formatted file with genomic
            annotations, see :func:`annotateGenome`.
       4. ``geneset_flat``, a flattened gene set in :term:`gtf` format, see
            :func:`buildFlatGeneSet`.
    outfile : string
       Output filename in :term:`bed` format.
    distance : int
       Merge adajcent features of the same type within this distance.
    '''

    repeats_gff, rna_gff, annotations_gtf, utr_gtf, intron_gtf = infiles

    tmpfile = P.getTempFilename(shared=True)
    tmpfiles = ["%s_%i" % (tmpfile, x) for x in range(4)]

    # add ENSEMBL annotations
    statement = """
    zcat %(annotations_gtf)s
    | cgat gtf2gtf
    --method=sort --sort-order=gene
    | cgat gtf2gtf
    --method=merge-exons --log=%(outfile)s.log
    | cgat gff2bed
    --set-name=gene_biotype --is-gtf
    --log=%(outfile)s.log
    | sort -k 1,1 -k2,2n
    | cgat bed2bed --method=merge --merge-by-name
    --merge-distance=%(distance)i --log=%(outfile)s.log
    > %(tmpfile)s_0
    """
    P.run()

    # rna
    statement = '''
    zcat %(repeats_gff)s %(rna_gff)s
    | cgat gff2bed --set-name=family --is-gtf -v 0
    | sort -k1,1 -k2,2n
    | cgat bed2bed --method=merge --merge-by-name
    --merge-distance=%(distance)i --log=%(outfile)s.log
    > %(tmpfile)s_1'''
    P.run()

    # utr
    statement = '''zcat %(utr_gtf)s
    | cgat gff2bed --is-gtf --set-name=feature
    | sort -k1,1 -k2,2n
    | cgat bed2bed --method=merge --merge-by-name
    --merge-distance=%(distance)i --log=%(outfile)s.log
    > %(tmpfile)s_2'''
    P.run()

    # intron
    statement = '''zcat %(intron_gtf)s
    | cgat gff2bed --is-gtf --set-name=feature
    | sort -k1,1 -k2,2n
    | cgat bed2bed --method=merge --merge-by-name
    --merge-distance=%(distance)i --log=%(outfile)s.log
    > %(tmpfile)s_3'''
    P.run()

    # sort and merge
    # remove strand information as bedtools
    # complains if there are annotations with
    # different number of field
    files = " ".join(tmpfiles)
    statement = '''
    sort --merge -k1,1 -k2,2n %(files)s
    | cut -f 1-4
    | gzip
    > %(outfile)s
    '''
    P.run()

    for x in tmpfiles:
        os.unlink(x)


def buildFlatGeneSet(infile, outfile):
    '''build a flattened gene set.
    All transcripts in a gene are merged into a single transcript by
    combining overlapping exons.
    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       Output filename in :term:`gtf` format.
    '''
    # sort by contig+gene, as in refseq gene sets, genes on
    # chr_random might contain the same identifier as on chr
    # and hence merging will fail.
    # --permit-duplicates is set so that these cases will be
    # assigned new merged gene ids.

    statement = """gunzip
    < %(infile)s
    | awk '$3 == "exon"'
    | grep "transcript_id"
    | cgat gtf2gtf
    --method=sort
    --sort-order=contig+gene
    --log=%(outfile)s.log
    | cgat gtf2gtf
    --method=merge-exons
    --permit-duplicates
    --log=%(outfile)s.log
    | cgat gtf2gtf
    --method=set-transcript-to-gene
    --log=%(outfile)s.log
    | cgat gtf2gtf
    --method=sort
    --sort-order=position+gene
    --log=%(outfile)s.log
    | gzip
    > %(outfile)s
        """
    P.run()


def loadGeneInformation(infile, outfile, only_proteincoding=False):
    '''load gene-related attributes from :term:`gtf` file into database.
    This method takes transcript-associated features from an
    :term:`gtf` file and collects the gene-related attributes in the
    9th column of the gtf file, ignoring exon_id, transcript_id,
    transcript_name, protein_id and exon_number.
    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       Output filename, contains logging information. The
       table name is derived from the filename of outfile.
    only_proteincoding : bool
       If True, only consider protein coding genes.
    '''

    job_memory = "4G"
    table = P.toTable(outfile)

    if only_proteincoding:
        filter_cmd = """cgat gtf2gtf
        --method=filter --filter-method=proteincoding""" % PARAMS
    else:
        filter_cmd = "cat"

    load_statement = P.build_load_statement(
        table,
        options="--add-index=gene_id "
        "--add-index=gene_name "
        "--map=gene_name:str")

    statement = '''
    zcat %(infile)s
    | %(filter_cmd)s
    | grep "transcript_id"
    | cgat gtf2gtf
    --method=sort --sort-order=gene+transcript
    | cgat gtf2tsv
    --attributes-as-columns --output-only-attributes -v 0
    | python %(toolsdir)s/csv_cut.py
    --remove exon_id transcript_id transcript_name protein_id exon_number
    | %(pipeline_scriptsdir)s/hsort 1
    | uniq
    | %(load_statement)s
    > %(outfile)s'''

    P.run()
