"""
==============================
Pipeline somatic exome vardict annotate
==============================

Overview
========

This pipeline takes output from vardict and annotates, filters and converts to table.

Input files
-----------

vardict .vcf files

file naming:
input files from vardict, e.g.:
vardict_CM0403_1.vcf

Code
====

"""
from ruffus import *
import sys
import os
import gzip
import vcf
import collections
import CGAT.IOTools as IOTools
import CGAT.Experiment as E
import CGATPipelines.Pipeline as P

# load options from the config file
PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])


################# Filter vardict ################################
    
@follows(mkdir("vardict_filter"))
@transform("vardict/*.vcf", 
           regex(r"vardict/(.*).vcf"),
                r"vardict_filter/\1.PASS.vcf")
def vardict_filter_PASS(infile,outfile):
    '''filter vardict to include only PASS filter'''
    statement = '''bcftools view
                    -f PASS \
                    -Ov -o %(outfile)s
                    %(infile)s'''    
    P.run()
    
@transform(vardict_filter_PASS, 
           regex(r"vardict_filter/(.*).PASS.vcf"),
                r"vardict_filter/\1.PASS_filtered.vcf")
def vardict_filter2(infile,outfile):
    '''filter vardict include StrongSomatic and LikelySomatic status'''
    statement = '''bcftools view
                    -i "INFO/STATUS='LikelySomatic'||INFO/STATUS='StrongSomatic'"
                    -Ov -o %(outfile)s
                    %(infile)s'''
    P.run()
    
################## Variant annotation ###################

@follows(mkdir("vardict_annovar_annotation"))
@transform(vardict_filter2,
           regex(r"vardict_filter/(.*).PASS_filtered.vcf"),
           r"vardict_annovar_annotation/\1.bcf_normalised.vcf")
def bcftools(infile, outfile):
    '''Splits multiallelic sites into multiple lines'''
    statement = '''bcftools norm -m-both
                    -o %(outfile)s
                    %(infile)s''' 
    P.run()

    
@transform(bcftools, 
           regex(r"vardict_annovar_annotation/(.*).bcf_normalised.vcf"),
                r"vardict_annovar_annotation/\1.hg38_multianno.vcf")
def annovar_annotate(infile,outfile):
    '''annotate using Annovar vcf file input'''
    basename = P.snip(outfile, ".hg38_multianno.vcf")
    statement = '''module() {  eval `/usr/bin/modulecmd bash $*`; } &&
                   module load annovar/2018-03-06 &&
                    table_annovar.pl
                    %(infile)s
                    /databank/indices/annovar/humandb
                    --buildver hg38
                    --remove
                    --outfile %(basename)s
                    -protocol %(annovar_protocol)s
                    -operation %(annovar_operation)s
                    -vcfinput'''
    P.run()

################## Vcf to table ###################

@follows(mkdir("vardict_table_variants"))    
@transform(annovar_annotate, 
           regex(r"vardict_annovar_annotation/(.*).hg38_multianno.vcf"),
                r"vardict_table_variants/\1.tsv")
def VariantsToTable(infile,outfile):
    '''converts the vcf file into a tab-separated file while splitting the INFO and FORMAT field'''
    vcf_reader = vcf.Reader(open(infile))
    # requires installation of pyvcf

    list_IDs = []
    list_desc = []
    list_cstring = []
    list_cstring2 = []
        
    for k,v in vcf_reader.infos.items():
        if k != "Samples":
            list_IDs.append(k)
            list_cstring.append("-F %s" % k)
            
    for k,v in vcf_reader.formats.items():
        if k != "Samples":
            list_IDs.append(k)
            list_cstring2.append("-GF %s" % k)
    
    cstring = " ".join(list_cstring)
    cstring2 = " ".join(list_cstring2)
    
    statement = '''gatk VariantsToTable
                    -R %(bwa_index)s
                   -V %(infile)s
                   -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER 
                   %(cstring)s
                   %(cstring2)s
                   --show-filtered=true
                   -O %(outfile)s'''
    P.run()

    
@transform(VariantsToTable,
           regex("vardict_table_variants/(\S+).tsv"),
           r"vardict_table_variants/\1.table.tsv")
def Table(infile,outfile):
    '''replace \x3d with = and \x3b with ;'''
    
    statement = '''sed -e 's/\\\\x3d/=/g' %(infile)s |
                    sed -e 's/\\\\x3b/;/g' > %(outfile)s'''
    P.run() 

                    
################## Abbreviations #######################

@merge(annovar_annotate, 
       "vardict_table_variants/abbreviations.tsv")
def Abbreviations(infiles,outfile):
    '''get a list of abbreviations'''
    infile = infiles[0]
    
    with IOTools.openFile(outfile, "w") as outf:
        with IOTools.openFile(infile, "r") as inf:
            for line in inf.readlines():
                if line.startswith('##FILTER'):
                    outf.write(line)
                if line.startswith('##FORMAT'):
                    outf.write(line)
                if line.startswith('##INFO'):
                    outf.write(line)
    
######################################################################
   
@follows(vardict_filter_PASS, vardict_filter2)
def Filtering():
    pass
 
@follows(bcftools, annovar_annotate)
def Annotation():
    pass
 
@follows(VariantsToTable, Table, Abbreviations)
def Variant_Tables():
    pass
 
  
@follows(Filtering, Annotation, Variant_Tables)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))