"""
==============================
Pipeline Pindel
==============================

Overview
========

This pipeline calls somatic indels from matched tumour and normal 
whole exome sequencing data using the Pindel.
Pindel can detect breakpoints of large deletions, medium sized insertions,
inversions, tandem duplications and other structural variants at single-based
resolution from next-gen sequence data.

Input files
-----------

paired end bam files after mark duplicates

file naming:
input files from Pipeline somatic exome gatk4 - mark duplicates, e.g.:
CM0403_1-control.md.bam
CM0403_1-tumour.md.bam

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

################## Pindel ########## 

@follows(mkdir("pindel"))
@transform("mark_duplicates/*.md.bam",
           regex(r"mark_duplicates/(\S+)(-tumour).md.bam"),
                r"pindel/\1_config.txt")
def pindel_config(infile, outfile):
    
    basename = P.snip(os.path.basename(infile), ".md.bam")
    size_metrics = "bqsr_qc/" + basename + ".insert_size_metrics"
    
    with IOTools.openFile(size_metrics, "r") as inf:
        lines = inf.readlines()
        values = lines[7].split("\t")
        median_size = values[0]
        inf.close()
                
    with IOTools.openFile(outfile, "w") as outf:
        line = infile  + "\t" + median_size + "\t" + basename + "\n" 
        outf.write(line)
       
@transform(pindel_config,
           regex(r"pindel/(.*)_config.txt"),
           r"pindel/\1_pindel_BP")
def pindel(infile, outfile):
    '''run pindel'''
    job_threads = 2    
    basename = P.snip(os.path.basename(infile),"_config.txt")
    outfile2  = P.snip(outfile, "_BP")
    statement = '''pindel
                 -f %(bwa_index)s
                 -i %(infile)s
                 -c %(pindel_chromosome)s
                 -o %(outfile2)s'''
    P.run()   
    
@transform(pindel,
           regex(r"pindel/(.*)_pindel_BP"),
           r"pindel/\1_pindel.vcf")
def pindel2vcf(infile, outfile):
    basename = P.snip(outfile,".vcf")
    statement = '''pindel2vcf
                    -P %(basename)s
                    -G
                    -r %(bwa_index)s
                    -R %(pindel_genome_build)s
                    -d %(pindel_genome_date)s
                    --min_supporting_reads %(pindel_mut_reads)s
                    -v %(outfile)s'''
                    
    P.run()
    
 ################## Variant annotation Pindel ###################
    

@follows(mkdir("annovar_annotation_pindel"))     
@transform(pindel2vcf, 
           regex(r"pindel/(.*)_pindel.vcf"),
                r"annovar_annotation_pindel/\1.hg38_multianno.vcf")
def annovar_annotate_pindel(infile,outfile):
    '''annotate variants using Annovar vcf file input'''
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

################## Make tables for Pindel ###################
    
@follows(mkdir("table_variants_pindel"))
@transform(annovar_annotate_pindel,
           regex(r"annovar_annotation_pindel/(.*).hg38_multianno.vcf"),
                r"table_variants_pindel/\1_pindel.tsv")

def VariantsToTablePindel(infile, outfile):
    '''converts the vcf file into a tab-separated file while splitting
    the INFO and FORMAT fields'''
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
    
@transform(VariantsToTablePindel,
           regex("table_variants_pindel/(\S+)_pindel.tsv"),
           r"table_variants_pindel/\1_pindel.table.tsv")
def Table_pindel(infile,outfile):
    '''replace \x3d by = and \x3b by ;'''
    
    statement = '''sed -e 's/\\\\x3d/=/g' %(infile)s |
                    sed -e 's/\\\\x3b/;/g' > %(outfile)s'''
    P.run()   
    
################## Make tables for FLT3-ITD ###################
@transform(Table_pindel,
           regex("table_variants_pindel/(\S+)_pindel.table.tsv"),
           r"table_variants_pindel/\1_FLT3.tsv")
def FLT3_table(infile, outfile):
    '''limits pindel table to FLT3_ITD'''

    with IOTools.openFile(outfile, "w") as outf:
       with IOTools.openFile(infile, "r") as inf:
           for line in inf.readlines():
               if line.startswith('CHROM'):
                   outf.write(line)
                    
               if line.startswith('chr'):
                   values = line.split("\t")
                   
                   if 'FLT3' in values[16]:
                       if 'exonic' in values[15]:
                           if 'nonframeshift_insertion' in values[18]:
                               outf.write(line)

#################################################
        
@follows(pindel_config, pindel, pindel2vcf)
def Pindel():
    pass

@follows(annovar_annotate_pindel)
def Annotation_Pindel():
    pass

@follows(VariantsToTablePindel, Table_pindel, FLT3_table)
def Variant_Tables_Pindel():
    pass

@follows(Pindel, Annotation_Pindel, Variant_Tables_Pindel)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))