"""
==============================
Pipeline variant overlap
==============================

Overview
========

This pipeline takes output from Mutect2, Strelka2 and vardict and looks for variant overlap.

Input files
-----------

Mutect2 .filtered_table.tsv
Strelka2 _snvs.tsv _indels.tsv
vardict .table.tsv

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
import pandas as pd
import CGAT.IOTools as IOTools
import CGAT.Experiment as E
import CGATPipelines.Pipeline as P

# load options from the config file
PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])


################# Variant coordinate lists ################################
    
@follows(mkdir("variant_lists"))
@transform("table_variants/*.filtered_table.tsv", 
           regex(r"table_variants/(.*).filtered_table.tsv"),
                r"variant_lists/\1_mutect2.tsv")
def mutect_variants(infile,outfile):
    '''table with coordinates of Mutect2 variants'''
    statement = '''tail -n +2 %(infile)s | cut -f 1-5 | sort > %(outfile)s'''    
    P.run()
    
@follows(mkdir("variant_lists"))
@transform("strelka_table_variants/*_snvs.vaf.table.tsv", 
           regex(r"strelka_table_variants/(.*)_snvs.vaf.table.tsv"),
                r"variant_lists/\1_strelka.tsv")
def strelka_variants(infile,outfile):
    '''table with coordinates of Strelka2 variants'''
    infile_indels = infile.replace("snvs", "indels")
    statement = '''cat <(tail -n +2 %(infile)s | cut -f 1-5) <(tail -n +2 %(infile_indels)s | cut -f 1-5) | sort > %(outfile)s'''    
    P.run()
    
@follows(mkdir("variant_lists"))
@transform("vardict_table_variants/*.table.tsv", 
           regex(r"vardict_table_variants/vardict_(.*).table.tsv"),
                r"variant_lists/\1_vardict.tsv")
def vardict_variants(infile,outfile):
    '''table with coordinates of Vardict variants'''

    statement = '''tail -n +2 %(infile)s | cut -f 1-5 | sort > %(outfile)s'''    
    P.run()

@follows(strelka_variants, vardict_variants)
@transform(mutect_variants,
           regex(r"variant_lists/(.*)_mutect2.tsv"),
                r"variant_lists/\1_variants.tsv")
def combined_variants(infile,outfile):
    '''table with combined coordinates'''
    infile_strelka = infile.replace("mutect2", "strelka")
    infile_vardict = infile.replace("mutect2", "vardict")
    statement = '''cat %(infile)s %(infile_strelka)s %(infile_vardict)s | sort -k1,1V -k2,2n |
                    uniq -c | grep -v "1 " | sed 's/[^      ]*      //' | sed 's/[^ ]* //' > %(outfile)s'''    
    P.run()
    
################# File subsetting ################################

@follows(combined_variants, mkdir("variant_subset"))
@transform(mutect_variants, 
           regex(r"variant_lists/(.*)_mutect2.tsv"),
                r"variant_subset/\1_mutect2_selected.tsv")
def mutect_subset(infile,outfile):
    '''subset the Mutect2 output with the combined variant list'''
    basename = P.snip(os.path.basename(infile),"_mutect2.tsv")
    infile_mutect = "table_variants/" + basename + ".filtered_table.tsv"
    infile_combined_variants = "variant_lists/" + basename + "_variants.tsv"
    statement = '''cat <(head -n 1 %(infile_mutect)s) <(grep -Fwf %(infile_combined_variants)s %(infile_mutect)s) > %(outfile)s'''    
    P.run()
    
@follows(combined_variants, mkdir("variant_subset"))
@transform(strelka_variants, 
           regex(r"variant_lists/(.*)_strelka.tsv"),
                r"variant_subset/\1_strelka_snvs_selected.tsv")
def strelka_snv_subset(infile,outfile):
    '''subset the Strelka2 snv output with the combined variant list'''
    basename = P.snip(os.path.basename(infile),"_strelka.tsv")
    infile_strelka = "strelka_table_variants/" + basename + "_snvs.vaf.table.tsv"
    infile_combined_variants = "variant_lists/" + basename + "_variants.tsv"
    statement = '''cat <(head -n 1 %(infile_strelka)s) <(grep -Fwf %(infile_combined_variants)s %(infile_strelka)s) > %(outfile)s'''    
    P.run()
    
@follows(combined_variants, mkdir("variant_subset"))
@transform(strelka_variants, 
           regex(r"variant_lists/(.*)_strelka.tsv"),
                r"variant_subset/\1_strelka_indels_selected.tsv")
def strelka_indel_subset(infile,outfile):
    '''subset the Strelka2 indels output with the combined variant list'''
    basename = P.snip(os.path.basename(infile),"_strelka.tsv")
    infile_strelka = "strelka_table_variants/" + basename + "_indels.vaf.table.tsv"
    infile_combined_variants = "variant_lists/" + basename + "_variants.tsv"
    statement = '''cat <(head -n 1 %(infile_strelka)s) <(grep -Fwf %(infile_combined_variants)s %(infile_strelka)s) > %(outfile)s'''    
    P.run()
    
@follows(combined_variants, mkdir("variant_subset"))
@transform(vardict_variants, 
           regex(r"variant_lists/(.*)_vardict.tsv"),
                r"variant_subset/\1_vardict_selected.tsv")
def vardict_subset(infile,outfile):
    '''subset the Vardict output with the combined variant list'''
    basename = P.snip(os.path.basename(infile),"_vardict.tsv")
    infile_vardict = "vardict_table_variants/vardict_" + basename + ".table.tsv"
    infile_combined_variants = "variant_lists/" + basename + "_variants.tsv"
    statement = '''cat <(head -n 1 %(infile_vardict)s) <(grep -Fwf %(infile_combined_variants)s %(infile_vardict)s) > %(outfile)s'''    
    P.run()
    
####################### Additioanal FLT3 variant indels ##################
    
@follows(strelka_variants, vardict_variants, mkdir("variant_FLT3"))
@transform(mutect_variants, 
           regex(r"variant_lists/(.*)_mutect2.tsv"),
                r"variant_FLT3/\1_variants_chr13.tsv")
def chr13_variants(infile,outfile):
    '''table with chr13 (containing FLT3 gene) variant coordinates found by just one variant caller'''
    infile_strelka = infile.replace("mutect2", "strelka")
    infile_vardict = infile.replace("mutect2", "vardict")
    statement = '''cat %(infile)s %(infile_strelka)s %(infile_vardict)s | sort -k1,1V -k2,2n |
                    uniq -c | grep "1 " | sed 's/[^      ]*      //' | sed 's/[^ ]* //' | grep chr13 > %(outfile)s'''    
    P.run()
    
@follows(chr13_variants)
@transform(mutect_variants, 
           regex(r"variant_lists/(.*)_mutect2.tsv"),
                r"variant_FLT3/\1_mutect2_FLT3.tsv")
def mutect_FLT3(infile,outfile):
    '''subset the Mutect2 output with the chr13 FLT3 variant list'''
    basename = P.snip(os.path.basename(infile),"_mutect2.tsv")
    infile_mutect = "table_variants/" + basename + ".filtered_table.tsv"
    infile_chr13 = "variant_FLT3/" + basename + "_variants_chr13.tsv"
    statement = '''cat <(head -n 1 %(infile_mutect)s) <(grep -Fwf %(infile_chr13)s %(infile_mutect)s | grep FLT3) > %(outfile)s'''    
    P.run()
    
@follows(chr13_variants)
@transform(strelka_variants, 
           regex(r"variant_lists/(.*)_strelka.tsv"),
                r"variant_FLT3/\1_strelka_FLT3.tsv")
def strelka_FLT3(infile,outfile):
    '''subset the Strelka2 indels output with the chr13 FLT3 variant list'''
    basename = P.snip(os.path.basename(infile),"_strelka.tsv")
    infile_strelka = "strelka_table_variants/" + basename + "_indels.vaf.table.tsv"
    infile_chr13 = "variant_FLT3/" + basename + "_variants_chr13.tsv"
    statement = '''cat <(head -n 1 %(infile_strelka)s) <(grep -Fwf %(infile_chr13)s %(infile_strelka)s | grep FLT3) > %(outfile)s'''    
    P.run()
    
@follows(chr13_variants)
@transform(vardict_variants, 
           regex(r"variant_lists/(.*)_vardict.tsv"),
                r"variant_FLT3/\1_vardict_FLT3.tsv")
def vardict_FLT3(infile,outfile):
    '''subset the Vardict output with the chr13 FLT3 variant list'''
    basename = P.snip(os.path.basename(infile),"_vardict.tsv")
    infile_vardict = "vardict_table_variants/vardict_" + basename + ".table.tsv"
    infile_chr13 = "variant_FLT3/" + basename + "_variants_chr13.tsv"
    statement = '''cat <(head -n 1 %(infile_vardict)s) <(grep -Fwf %(infile_chr13)s %(infile_vardict)s | grep FLT3) > %(outfile)s'''    
    P.run()

######################## Final list #########################################

@follows(mkdir("variant_list_final"))   
@transform(mutect_subset,
           regex("variant_subset/(.*)_mutect2_selected.tsv"),
           r"variant_list_final/\1_mutect2_cols.tsv")
def mutect_cols(infile,outfile):
    '''select columns from Mutect list''' 
    
    with IOTools.openFile(outfile, "w") as outf:    
        with IOTools.openFile(infile, "r") as inf:
            for line in inf.readlines():
                if line.startswith('CHROM'):
                    line = line.rstrip('\n')
                    header = line.split("\t")
                    line = header[0] + "\t" + header[1] + "\t" + header[3] + "\t" + header[4] + "\t" + header[6] + "\t" + header[20] + "\t" + header[21] +"\t" + header[23] + "\t" + header[24] + "\t" + header[26] + "\t" + header[38] + "\t" + header[106] + "\t" + header[107] + "\t" + header[129] +"\t" + header[130] + "\tVariant_caller\n"
                    outf.write(line)
    
                if line.startswith('chr'):
                    line = line.rstrip('\n')
                    values = line.split("\t")
                    line = values[0] + "\t" + values[1] + "\t" + values[3] + "\t" + values[4] + "\t" + values[6] + "\t" + values[20] + "\t" + values[21] + "\t" + values[23] + "\t" + values[24] + "\t" + values[26] + "\t" + values[38] + "\t" + values[106] + "\t" + values[107] + "\t" + values[129] + "\t" + values[130] + "\tMutect2\n"
                    outf.write(line)
                    
@follows(mkdir("variant_list_final"))   
@transform(strelka_snv_subset,
           regex("variant_subset/(.*)_strelka_snvs_selected.tsv"),
           r"variant_list_final/\1_strelka_snvs_cols.tsv")
def strelka_snv_cols(infile,outfile):
    '''select columns from Strelka snv list''' 
    
    with IOTools.openFile(outfile, "w") as outf:    
        with IOTools.openFile(infile, "r") as inf:
            for line in inf.readlines():
                if line.startswith('CHROM'):
                    line = line.rstrip('\n')
                    header = line.split("\t")
                    line = header[0] + "\t" + header[1] + "\t" + header[3] + "\t" + header[4] + "\t" + header[6] + "\t" + header[23] + "\t" + header[24] +"\t" + header[26] + "\t" + header[27] + "\t" + header[29] + "\t" + header[41] + "\t" + header[110] + "\t" + header[-1] + "\t" + header[119] +"\t" + header[-2] + "\tVariant_caller\n"
                    outf.write(line)
    
                if line.startswith('chr'):
                    line = line.rstrip('\n')
                    values = line.split("\t")
                    line = values[0] + "\t" + values[1] + "\t" + values[3] + "\t" + values[4] + "\t" + values[6] + "\t" + values[23] + "\t" + values[24] + "\t" + values[26] + "\t" + values[27] + "\t" + values[29] + "\t" + values[41] + "\t" + values[110] + "\t" + values[-1] + "\t" + values[119] + "\t" + values[-2] + "\tStrelka2_snvs\n"
                    outf.write(line)
                    
@follows(mkdir("variant_list_final"))   
@transform(strelka_indel_subset,
           regex("variant_subset/(.*)_strelka_indels_selected.tsv"),
           r"variant_list_final/\1_strelka_indels_cols.tsv")
def strelka_indel_cols(infile,outfile):
    '''select columns from Strelka indel list''' 
    
    with IOTools.openFile(outfile, "w") as outf:    
        with IOTools.openFile(infile, "r") as inf:
            for line in inf.readlines():
                if line.startswith('CHROM'):
                    line = line.rstrip('\n')
                    header = line.split("\t")
                    line = header[0] + "\t" + header[1] + "\t" + header[3] + "\t" + header[4] + "\t" + header[6] + "\t" + header[23] + "\t" + header[24] +"\t" + header[26] + "\t" + header[27] + "\t" + header[29] + "\t" + header[41] + "\t" + header[110] + "\t" + header[-1] + "\t" + header[120] +"\t" + header[-2] + "\tVariant_caller\n"
                    outf.write(line)
    
                if line.startswith('chr'):
                    line = line.rstrip('\n')
                    values = line.split("\t")
                    line = values[0] + "\t" + values[1] + "\t" + values[3] + "\t" + values[4] + "\t" + values[6] + "\t" + values[23] + "\t" + values[24] + "\t" + values[26] + "\t" + values[27] + "\t" + values[29] + "\t" + values[41] + "\t" + values[110] + "\t" + values[-1] + "\t" + values[120] + "\t" + values[-2] + "\tStrelka2_indels\n"
                    outf.write(line)
                    
@follows(mkdir("variant_list_final"))   
@transform(vardict_subset,
           regex("variant_subset/(.*)_vardict_selected.tsv"),
           r"variant_list_final/\1_vardict_cols.tsv")
def vardict_cols(infile,outfile):
    '''select columns from Vardict list''' 
    
    with IOTools.openFile(outfile, "w") as outf:    
        with IOTools.openFile(infile, "r") as inf:
            for line in inf.readlines():
                if line.startswith('CHROM'):
                    line = line.rstrip('\n')
                    header = line.split("\t")
                    line = header[0] + "\t" + header[1] + "\t" + header[3] + "\t" + header[4] + "\t" + header[20] + "\t" + header[22] + "\t" + header[23] +"\t" + header[25] + "\t" + header[26] + "\t" + header[28] + "\t" + header[40] + "\t" + header[111] + "\t" + header[114] + "\t" + header[130] +"\t" + header[133] + "\tVariant_caller\n"
                    outf.write(line)
    
                if line.startswith('chr'):
                    line = line.rstrip('\n')
                    values = line.split("\t")
                    line = values[0] + "\t" + values[1] + "\t" + values[3] + "\t" + values[4] + "\t" + values[20] + "\t" + values[22] + "\t" + values[23] + "\t" + values[25] + "\t" + values[26] + "\t" + values[28] + "\t" + values[40] + "\t" + values[111] + "\t" + values[114] + "\t" + values[130] + "\t" + values[133] + "\tVardict\n"
                    outf.write(line)
                    
############################ Variant filtering #################################
                    
@follows(strelka_snv_cols, strelka_indel_cols, vardict_cols, mkdir("variant_filter_final"))
@transform(mutect_cols, 
           regex(r"variant_list_final/(.*)_mutect2_cols.tsv"),
                r"variant_filter_final/\1_table.tsv")
def combined_table(infile,outfile):
    '''table with combined variants'''
    infile_strelka_snvs = infile.replace("mutect2", "strelka_snvs")
    infile_strelka_indels = infile.replace("mutect2", "strelka_indels")
    infile_vardict = infile.replace("mutect2", "vardict")
    statement = '''cat %(infile)s <(tail -n +2 %(infile_strelka_snvs)s) 
                    <(tail -n +2 %(infile_strelka_indels)s) <(tail -n +2 %(infile_vardict)s) | 
                    sort -k1,1V -k2,2n -k16,16 | sed 's/nonsynonymous/missense/g' | grep -v 'synonymous\|ncRNA' | 
                    sed 's/missense/nonsynonymous/g' | grep 'CHROM\|exonic\|splicing' > %(outfile)s'''  
    P.run()
    
@transform(combined_table,
           regex("variant_filter_final/(.*)_table.tsv"),
           r"variant_filter_final/\1.vaf.filter.tsv")
def vaf_filter(infile,outfile):
    '''filter the table for VAF''' 
    
    #AML_table = pd.read_csv("/t1-data/user/stoilova/WES_GATK4_pipeline/AA_bed_files/AML_genes.txt", sep="\t")
    AML_table = pd.read_csv(PARAMS["filter_aml_genes"], sep="\t")
    AML_list = list(AML_table.AML_genes)
    
    with IOTools.openFile(outfile, "w") as outf:    
        with IOTools.openFile(infile, "r") as inf:
            for line in inf.readlines():
                if line.startswith('CHROM'):
                    outf.write(line)
    
                if line.startswith('chr'):
                    values = line.split("\t")
                    if values[6] in AML_list:
                        if float(values[12]) < 0.2:
                            if float(values[14])> 0.05:
                                outf.write(line)
                    if not values[6] in AML_list:
                        if float(values[12]) < 0.01:
                            if float(values[14])> 0.05:
                                outf.write(line)
                                
@follows(mkdir("variant_recurrent_AML_genes"))
@transform(vaf_filter,
           regex("variant_filter_final/(.*).vaf.filter.tsv"),
           r"variant_recurrent_AML_genes/\1.recurrent_AML.tsv")
def recurrent_aml_genes(infile,outfile):
    '''filter the table for recurant AML genes'''
    statement = '''cat <(head -n 1 %(infile)s) <(grep -Fwf %(filter_aml_recurrent_genes)s %(infile)s) > %(outfile)s'''    
    P.run()
  
#####################################################################################
@follows(mkdir("variants_rel_patients"))   
@transform(vaf_filter,
           regex("variant_filter_final/(.*)_2.vaf.filter.tsv"),
           r"variants_rel_patients/\1_1.sample.tsv")
def rel_patient_vars_dx(infile,outfile):
    '''adds Dx to Dx variants from relapse patients''' 
    infile_dx = infile.replace("_2.vaf", "_1.vaf")
    
    with IOTools.openFile(outfile, "w") as outf:    
        with IOTools.openFile(infile_dx, "r") as inf:
            for line in inf.readlines():
                if line.startswith('CHROM'):
                    line = line.rstrip('\n')
                    line = line + "\tSample\n"
                    outf.write(line)
    
                if line.startswith('chr'):
                    line = line.rstrip('\n')
                    line = line + "\tDiagnosis\n"
                    outf.write(line)
                    
@follows(mkdir("variants_rel_patients"))   
@transform(vaf_filter,
           regex("variant_filter_final/(.*)_2.vaf.filter.tsv"),
           r"variants_rel_patients/\1_2.sample.tsv")
def rel_patient_vars_rel(infile,outfile):
    '''adds Rel to Rel variants from relapse patients''' 
    
    with IOTools.openFile(outfile, "w") as outf:    
        with IOTools.openFile(infile, "r") as inf:
            for line in inf.readlines():
                if line.startswith('CHROM'):
                    line = line.rstrip('\n')
                    line = line + "\tSample\n"
                    outf.write(line)
    
                if line.startswith('chr'):
                    line = line.rstrip('\n')
                    line = line + "\tRelapse\n"
                    outf.write(line)
                    

@follows(rel_patient_vars_rel)
@transform(rel_patient_vars_dx, 
           regex(r"variants_rel_patients/(.*)_1.sample.tsv"),
                r"variants_rel_patients/\1.patient_vars.tsv")
def rel_patient_vars(infile,outfile):
    '''table with combined variants from Dx and Rel from relapse patients'''
    infile_rel = infile.replace("_1.sample", "_2.sample")
    statement = '''cat %(infile)s <(tail -n +2 %(infile_rel)s) | 
                    sort -k1,1V -k2,2n -k17,17 -k16,16 > %(outfile)s'''  
    P.run() 
    

######################################################################
 
@follows(mutect_variants, strelka_variants, vardict_variants, combined_variants)
def variants():
    pass

@follows(mutect_subset, strelka_snv_subset, strelka_indel_subset, vardict_subset)
def subsetting():
    pass

@follows(chr13_variants, mutect_FLT3, strelka_FLT3, vardict_FLT3)
def FLT3():
    pass

@follows(mutect_cols, strelka_snv_cols, strelka_indel_cols, vardict_cols)
def final_variant_list():
    pass

@follows(combined_table, vaf_filter, recurrent_aml_genes)
def final_variant_filter():
    pass

@follows(rel_patient_vars_dx, rel_patient_vars_rel, rel_patient_vars)
def rel_variants():
    pass

@follows(variants, subsetting, FLT3, final_variant_list, final_variant_filter)
def full():
    pass


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))