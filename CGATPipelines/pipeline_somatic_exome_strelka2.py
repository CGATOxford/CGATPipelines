"""
==============================
Pipeline somatic exome Strelka2
==============================

Overview
========

This pipeline calls somatic variants from matched tumour and normal 
whole exome sequencing data using the Strelka2.

Input files
-----------

paired end bam files after mark duplicates
bed file of exome enriched regions
bed file needs to be sorted, bgzipped and indexed, e.g.:
bedtools sort -i S07604514_Padded_full_hg38.bed > S07604514_Padded_full_hg38_sorted.bed
bgzip -c S07604514_Padded_full_hg38_sorted.bed > S07604514_Padded_full_hg38_sorted.bed.gz 
tabix S07604514_Padded_full_hg38_sorted.bed.gz

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

########## Pre-processing ########## 

@transform("mark_duplicates/*.md.bam",
           regex(r"mark_duplicates/(\S+).md.bam"),
           r"mark_duplicates/\1.md.bam.bai")
def index_md_bam(infile, outfile):
    '''index mark duplicates bam '''
    statement = '''samtools index %(infile)s'''
    P.run()
    
@follows(mkdir("strelka"))
@transform(index_md_bam,
           regex(r"mark_duplicates/(.*)(-tumour).md.bam.bai"),
           r"strelka/\1.pid")
def patientID(infiles, outfile):
    '''makes and empty file for patient ID based on the tumour samples'''
    to_cluster = False
    statement = '''touch %(outfile)s'''
    P.run()

############ Manta #####################################
    
@transform(patientID, 
           regex(r"strelka/(.*).pid"),
                r"strelka/\1_Manta/runWorkflow.py")
def manta_config(infile,outfile):
    '''configure Manta'''
    '''control must be called _1-control'''
    basename = P.snip(os.path.basename(infile),".pid")
    infile_tumour = "mark_duplicates/" + basename + "-tumour.md.bam"  
    basename2 = basename.split("_")
    infile_control = "mark_duplicates/" + basename2[0] + "_1-control.md.bam"
    runDir = "strelka/" + basename + "_Manta"
    statement = '''/t1-data/user/stoilova/Strelka2/manta-1.3.2.centos6_x86_64/bin/configManta.py 
                    --normalBam %(infile_control)s
                    --tumorBam %(infile_tumour)s
                    --referenceFasta %(bwa_index)s
                    --callRegions %(strelka_call_regions)s
                    --exome
                    --runDir %(runDir)s'''    
    P.run()
    
    
@transform(manta_config, 
           regex(r"strelka/(.*)_Manta/runWorkflow.py"),
                r"strelka/\1_Manta/results/variants/candidateSmallIndels.vcf.gz")
def run_manta(infile,outfile):
    '''run Manta to call structural variants and indels'''
    statement = '''%(infile)s -m local -j 8'''    
    P.run()
    

@follows(run_manta)
@transform(patientID, 
           regex(r"strelka/(.*).pid"),
                r"strelka/\1_Strelka/runWorkflow.py")
def strelka_config(infile,outfile):
    '''configure Strelka'''
    '''control must be called _1-control'''
    basename = P.snip(os.path.basename(infile),".pid")
    infile_tumour = "mark_duplicates/" + basename + "-tumour.md.bam"  
    basename2 = basename.split("_")
    infile_control = "mark_duplicates/" + basename2[0] + "_1-control.md.bam"
    manta_indels = "strelka/" + basename + "_Manta/results/variants/candidateSmallIndels.vcf.gz"
    runDir = "strelka/" + basename + "_Strelka"
    statement = '''/t1-data/user/stoilova/Strelka2/strelka-2.9.2/bin/configureStrelkaSomaticWorkflow.py 
                    --normalBam %(infile_control)s
                    --tumorBam %(infile_tumour)s
                    --referenceFasta %(bwa_index)s
                    --indelCandidates %(manta_indels)s
                    --callRegions %(strelka_call_regions)s
                    --exome
                    --runDir %(runDir)s'''    
    P.run()


@transform(strelka_config, 
           regex(r"strelka/(.*)_Strelka/runWorkflow.py"),
                r"strelka/\1_Strelka/results/variants/somatic.snvs.vcf")
def run_strelka(infile,outfile):
    '''run Strelka to call SNVs and indels'''
    outfile_2 = outfile + ".gz"
    outfile_indels = outfile.replace("snvs", "indels")
    outfile_indels_2 =outfile_2.replace("snvs", "indels")
    statement = '''%(infile)s -m local -j 8 &&
                    bgzip -d -c %(outfile_2)s > %(outfile)s &&
                    bgzip -d -c %(outfile_indels_2)s > %(outfile_indels)s'''    
    P.run()

################# Filter Strelka ################################
    
@follows(run_strelka, mkdir("strelka_filter"))
@transform(patientID, 
           regex(r"strelka/(.*).pid"),
                r"strelka_filter/\1.EVS.somatic.snvs.vcf")
def strelka_filter_snv(infile,outfile):
    '''filter Strelka SNVs include EVS >=5'''
    basename = P.snip(os.path.basename(infile),".pid")
    infile_snv = "strelka/" + basename + "_Strelka/results/variants/somatic.snvs.vcf"
    statement = '''bcftools view
                    -i 'INFO/SomaticEVS>=5'
                    -Ov -o %(outfile)s
                    %(infile_snv)s'''    
    P.run()
    
@transform(strelka_filter_snv, 
           regex(r"strelka_filter/(.*).EVS.somatic.snvs.vcf"),
                r"strelka_filter/\1.filtered.somatic.snvs.vcf")
def strelka_filter2_snv(infile,outfile):
    '''filter Strelka SNVs exclude LowDepth'''
    statement = '''bcftools view
                    -e 'FILTER="LowDepth"'
                    -Ov -o %(outfile)s
                    %(infile)s'''
    P.run()
    
@follows(run_strelka, mkdir("strelka_filter"))
@transform(patientID, 
           regex(r"strelka/(.*).pid"),
                r"strelka_filter/\1.PASS.somatic.indels.vcf")
def strelka_PASS_indel(infile,outfile):
    '''filter indels on Strelka PASS FILTER'''
    basename = P.snip(os.path.basename(infile),".pid")
    infile_snv = "strelka/" + basename + "_Strelka/results/variants/somatic.indels.vcf"
    statement = '''bcftools view
                    -f PASS
                    -Ov -o %(outfile)s
                    %(infile_snv)s'''    
    P.run()
    
################## Vcf processing ###################

@follows(mkdir("strelka_vcf_processing"))
@transform(strelka_filter2_snv,
           regex(r"strelka_filter/(.*).filtered.somatic.snvs.vcf"),
           r"strelka_vcf_processing/\1_snvs.bcf_normalised.vcf")
def bcftools_snv(infile, outfile):
    '''Splits multiallelic SNV sites into multiple lines'''
    statement = '''bcftools norm -m-both
                    -o %(outfile)s
                    %(infile)s''' 
    P.run()
    
@follows(mkdir("strelka_vcf_processing"))
@transform(strelka_PASS_indel,
           regex(r"strelka_filter/(.*).PASS.somatic.indels.vcf"),
           r"strelka_vcf_processing/\1_indels.bcf_normalised.vcf")
def bcftools_indels(infile, outfile):
    '''Splits multiallelic indel sites into multiple lines'''
    statement = '''bcftools norm -m-both
                    -o %(outfile)s
                    %(infile)s''' 
    P.run()

@transform(bcftools_snv,
           regex(r"strelka_vcf_processing/(.*).bcf_normalised.vcf"),
           r"strelka_vcf_processing/\1.gt.vcf")
def add_GT_snv(infile, outfile):
    '''add GT FORMAT required for annovar'''
    basename = P.snip(os.path.basename(infile),"_snvs.bcf_normalised.vcf")
    normal = basename + "_NORMAL"
    tumour = basename + "_TUMOUR"
    statement = '''bcftools view %(infile)s |
                awk -v NORMAL=%(normal)s -v TUMOR=%(tumour)s 'BEGIN {FS="\\t"; OFS=FS;} {if (NF < 10) print;
                else if ($1=="#CHROM") { $10=NORMAL; $11=TUMOR; print; }
                else {$9="GT:"$9; $10="0/0:"$10; $11="0/1:"$11; print;}}' |
                awk '!found && /##FORMAT/{print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\\"Genotype\\">"; found=1}1'
                > %(outfile)s''' 
    P.run()
    
@transform(bcftools_indels,
           regex(r"strelka_vcf_processing/(.*).bcf_normalised.vcf"),
           r"strelka_vcf_processing/\1.gt.vcf")
def add_GT_indels(infile, outfile):
    '''add GT FORMAT required for annovar'''
    basename = P.snip(os.path.basename(infile),"_indels.bcf_normalised.vcf")
    normal = basename + "_NORMAL"
    tumour = basename + "_TUMOUR"
    statement = '''bcftools view %(infile)s |
                awk -v NORMAL=%(normal)s -v TUMOR=%(tumour)s 'BEGIN {FS="\\t"; OFS=FS;} {if (NF < 10) print;
                else if ($1=="#CHROM") { $10=NORMAL; $11=TUMOR; print; }
                else {$9="GT:"$9; $10="0/0:"$10; $11="0/1:"$11; print;}}' |
                awk '!found && /##FORMAT/{print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\\"Genotype\\">"; found=1}1'
                > %(outfile)s''' 
    P.run()

################## Variant annotation ###################

@follows(mkdir("strelka_annovar_annotation"))     
@transform(add_GT_snv, 
           regex(r"strelka_vcf_processing/(.*).gt.vcf"),
                r"strelka_annovar_annotation/\1.hg38_multianno.vcf")
def annovar_annotate_snv(infile,outfile):
    '''annotate SNVs using Annovar vcf file input'''
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
    
@follows(mkdir("strelka_annovar_annotation"))     
@transform(add_GT_indels, 
           regex(r"strelka_vcf_processing/(.*).gt.vcf"),
                r"strelka_annovar_annotation/\1.hg38_multianno.vcf")
def annovar_annotate_indels(infile,outfile):
    '''annotate SNVs using Annovar vcf file input'''
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

@follows(mkdir("strelka_table_variants"))    
@transform(annovar_annotate_snv, 
           regex(r"strelka_annovar_annotation/(.*).hg38_multianno.vcf"),
                r"strelka_table_variants/\1.tsv")
def VariantsToTable_snv(infile,outfile):
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

@follows(mkdir("strelka_table_variants"))    
@transform(annovar_annotate_indels, 
           regex(r"strelka_annovar_annotation/(.*).hg38_multianno.vcf"),
                r"strelka_table_variants/\1.tsv")
def VariantsToTable_indels(infile,outfile):
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
    
@transform(VariantsToTable_snv,
           regex("strelka_table_variants/(\S+).tsv"),
           r"strelka_table_variants/\1.table.tsv")
def Table_snv(infile,outfile):
    '''replace \x3d with = and \x3b with ;'''
    
    statement = '''sed -e 's/\\\\x3d/=/g' %(infile)s |
                    sed -e 's/\\\\x3b/;/g' > %(outfile)s'''
    P.run()
    
@transform(VariantsToTable_indels,
           regex("strelka_table_variants/(\S+).tsv"),
           r"strelka_table_variants/\1.table.tsv")
def Table_indels(infile,outfile):
    '''replace \x3d with = and \x3b with ;'''
    
    statement = '''sed -e 's/\\\\x3d/=/g' %(infile)s |
                    sed -e 's/\\\\x3b/;/g' > %(outfile)s'''
    P.run()
    
    
@transform(Table_snv,
           regex("strelka_table_variants/(\S+).table.tsv"),
           r"strelka_table_variants/\1.vaf.table.tsv")
def vaf_snv(infile,outfile):
    '''add VAF to table''' 
    
    with IOTools.openFile(outfile, "w") as outf:    
        with IOTools.openFile(infile, "r") as inf:
            for line in inf.readlines():
                if line.startswith('CHROM'):
                    line = line.rstrip('\n') + "\tVAF_TUMOUR\tVAF_NORMAL\n"
                    outf.write(line)
    
                if line.startswith('chr'):
                    values = line.split("\t")
                    
                    if "C" in values[3]:
                        ref_T = values[-3].split(',')[0]
                        ref_N = values[-12].split(',')[0]
                    elif "T" in values[3]:
                        ref_T = values[-1].split(',')[0]
                        ref_N = values[-10].split(',')[0]
                    elif "A" in values[3]:
                        ref_T = values[-4].split(',')[0]
                        ref_N = values[-13].split(',')[0]
                    elif "G" in values[3]:
                        ref_T = values[-2].split(',')[0]
                        ref_N = values[-11].split(',')[0]
                    else:
                        raise ValueError("Invalid base")
                        
                    if "C" in values[4]:
                        alt_T = values[-3].split(',')[0]
                        alt_N = values[-12].split(',')[0]
                    elif "T" in values[4]:
                        alt_T = values[-1].split(',')[0]
                        alt_N = values[-10].split(',')[0]
                    elif "A" in values[4]:
                        alt_T = values[-4].split(',')[0]
                        alt_N = values[-13].split(',')[0]
                    elif "G" in values[4]:
                        alt_T = values[-2].split(',')[0]
                        alt_N = values[-11].split(',')[0]
                    else:
                        raise ValueError("Invalid base")
                        
                    vaf_T = int(alt_T) / (int(ref_T) + int(alt_T))
                    vaf_T = "{0:.3f}".format(vaf_T)
                    if int(ref_N) == 0:
                        if int(alt_N) == 0:
                            vaf_N = 0
                    elif int(ref_N) != 0:
                        vaf_N = int(alt_N) / (int(ref_N) + int(alt_N))
                    vaf_N = "{0:.3f}".format(vaf_N)
                    line = line.rstrip('\n') + "\t" + vaf_T + "\t" + vaf_N + "\n"
                    outf.write(line)


@transform(Table_indels,
           regex("strelka_table_variants/(\S+).table.tsv"),
           r"strelka_table_variants/\1.vaf.table.tsv")
def vaf_indels(infile,outfile):
    '''add VAF to table''' 
    
    with IOTools.openFile(outfile, "w") as outf:    
        with IOTools.openFile(infile, "r") as inf:
            for line in inf.readlines():
                if line.startswith('CHROM'):
                    line = line.rstrip('\n') + "\tVAF_TUMOUR\tVAF_NORMAL\n"
                    outf.write(line)
    
                if line.startswith('chr'):
                    values = line.split("\t")
                    ref_T = values[-7].split(',')[0]
                    ref_N = values[-17].split(',')[0]
                    alt_T = values[-6].split(',')[0] 
                    alt_N = values[-16].split(',')[0]
                        
                    vaf_T = int(alt_T) / (int(ref_T) + int(alt_T))
                    vaf_T = "{0:.3f}".format(vaf_T)
                    if int(ref_N) == 0:
                        if int(alt_N) == 0:
                            vaf_N = 0
                    else:
                        vaf_N = int(alt_N) / (int(ref_N) + int(alt_N))
                    vaf_N = "{0:.3f}".format(vaf_N)
                    line = line.rstrip('\n') + "\t" + vaf_T + "\t" + vaf_N + "\n"
                    outf.write(line)
                    
################## Abbreviations #######################

@merge(annovar_annotate_snv, 
       "strelka_table_variants/abbreviations_snv.tsv")
def Abbreviations_snv(infiles,outfile):
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
                    
@merge(annovar_annotate_indels,
       "strelka_table_variants/abbreviations_indels.tsv")
def Abbreviations_indels(infiles,outfile):
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

@follows(index_md_bam, patientID, manta_config, run_manta, strelka_config, run_strelka)
def Strelka2():
    pass
   
@follows(strelka_filter_snv, strelka_filter2_snv, strelka_PASS_indel)
def Filtering():
    pass

@follows(bcftools_snv, bcftools_indels, add_GT_snv, add_GT_indels)
def VCFprocessing():
    pass
 
@follows(annovar_annotate_snv, annovar_annotate_indels)
def Annotation():
    pass
 
@follows(VariantsToTable_snv, VariantsToTable_indels, Table_snv, Table_indels,
vaf_snv, vaf_indels, Abbreviations_snv, Abbreviations_indels)
def Variant_Tables():
    pass
 
  
@follows(Strelka2, Filtering, VCFprocessing, Annotation, Variant_Tables)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))