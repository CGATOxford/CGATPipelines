"""
PipelineExome.py - Tasks for variant calling
============================================

Reference
---------

"""
# Import modules
import os
import CGAT.IOTools as IOTools
import CGATPipelines.Pipeline as P
from CGATPipelines.Pipeline import cluster_runnable
import numpy as np
import pandas as pd
import CGAT.CSV as csv
import CGAT.VCF as VCF
import collections
import re
import urllib
import itertools
from bs4 import BeautifulSoup
from bs4 import NavigableString
from rpy2.robjects import pandas2ri
from rpy2.robjects import r as R


# Set PARAMS in calling module
PARAMS = {}


def getGATKOptions():
    return "-l mem_free=1.4G"


def makeSoup(address):
    sock = urllib.urlopen(address)
    htmlSource = sock.read()
    soup = BeautifulSoup(htmlSource)
    return soup

##############################################################################


def GATKReadGroups(infile, outfile, genome,
                   library="unknown", platform="Illumina",
                   platform_unit="1", track="unknown"):
    '''Reorders BAM according to reference fasta and adds read groups'''

    if track == 'unknown':
        track = P.snip(os.path.basename(infile), ".bam")
    tmpdir_gatk = P.getTempDir('.')
    job_options = getGATKOptions()
    job_threads = 3

    statement = '''ReorderSam
                    INPUT=%(infile)s
                    OUTPUT=%(tmpdir_gatk)s/%(track)s.reordered.bam
                    REFERENCE=%(genome)s
                    ALLOW_INCOMPLETE_DICT_CONCORDANCE=true
                    VALIDATION_STRINGENCY=SILENT ; checkpoint ;''' % locals()

    statement += '''samtools index %(tmpdir_gatk)s/%(track)s.reordered.bam ;
                    checkpoint ;''' % locals()

    statement += '''AddOrReplaceReadGroups
                    INPUT=%(tmpdir_gatk)s/%(track)s.reordered.bam
                    OUTPUT=%(outfile)s
                    RGLB=%(library)s
                    RGPL=%(platform)s
                    RGPU=%(platform_unit)s
                    RGSM=%(track)s
                    VALIDATION_STRINGENCY=SILENT ; checkpoint ;''' % locals()

    statement += '''samtools index %(outfile)s ;
                    checkpoint ;''' % locals()
    statement += '''rm -rf %(tmpdir_gatk)s ;''' % locals()

    P.run()

##############################################################################


def GATKIndelRealign(infile, outfile, genome, threads=4):
    '''Realigns BAMs around indels using GATK'''

    intervalfile = outfile.replace(".bam", ".intervals")
    job_options = getGATKOptions()
    job_threads = 3

    statement = '''GenomeAnalysisTK
                    -T RealignerTargetCreator
                    -o %(intervalfile)s
                    --num_threads %(threads)s
                    -R %(genome)s
                    -I %(infile)s; ''' % locals()

    statement += '''GenomeAnalysisTK
                    -T IndelRealigner
                    -o %(outfile)s
                    -R %(genome)s
                    -I %(infile)s
                    -targetIntervals %(intervalfile)s;''' % locals()
    P.run()

##############################################################################


def GATKBaseRecal(infile, outfile, genome, dbsnp, solid_options=""):
    '''Recalibrates base quality scores using GATK'''

    track = P.snip(os.path.basename(infile), ".bam")
    tmpdir_gatk = P.getTempDir('.')
    job_options = getGATKOptions()
    job_threads = 3

    statement = '''GenomeAnalysisTK
                    -T BaseRecalibrator
                    --out %(tmpdir_gatk)s/%(track)s.recal.grp
                    -R %(genome)s
                    -I %(infile)s
                    --knownSites %(dbsnp)s %(solid_options)s ;
                    checkpoint ;''' % locals()

    statement += '''GenomeAnalysisTK
                    -T PrintReads -o %(outfile)s
                    -BQSR %(tmpdir_gatk)s/%(track)s.recal.grp
                    -R %(genome)s
                    -I %(infile)s ;
                    checkpoint ;''' % locals()

    statement += '''rm -rf %(tmpdir_gatk)s ;''' % locals()
    P.run()

##############################################################################


def haplotypeCaller(infile, outfile, genome,
                    dbsnp, intervals, padding, options):
    '''Call SNVs and indels using GATK HaplotypeCaller in all members of a
    family together'''
    job_options = getGATKOptions()
    job_threads = 3

    statement = '''GenomeAnalysisTK
                    -T HaplotypeCaller
                    -o %(outfile)s
                    -R %(genome)s
                    -I %(infile)s
                    --dbsnp %(dbsnp)s
                    -L %(intervals)s
                    -ip %(padding)s
                    %(options)s''' % locals()
    P.run()

##############################################################################


def mutectSNPCaller(infile, outfile, mutect_log, genome, cosmic,
                    dbsnp, call_stats_out, job_memory, job_threads,
                    quality=20, max_alt_qual=150, max_alt=5,
                    max_fraction=0.05, tumor_LOD=6.3,
                    normal_panel=None, 
                    infile_matched=None,
                    gatk_key=None,
                    artifact=False):
    '''Call SNVs using Broad's muTect'''
    # TS. this is currently CGAT specific. How to generalise?

    job_memory_java = job_memory.lower()
    job_threads = 1

    statement = '''module load apps/java/jre1.6.0_26;
                   java -Xmx%(job_memory_java)s -jar
                   /ifs/apps/bio/muTect-1.1.4/muTect-1.1.4.jar
                   --analysis_type MuTect
                   --reference_sequence %(genome)s
                   --cosmic %(cosmic)s
                   --dbsnp %(dbsnp)s
                   --input_file:tumor %(infile)s
                   --out %(call_stats_out)s
                   --enable_extended_output
                   --vcf %(outfile)s
                ''' % locals()
    if artifact:
        statement += ''' --artifact_detection_mode '''

    if infile_matched:
        statement += '''--min_qscore %(quality)s
                        --gap_events_threshold 2
                        --max_alt_alleles_in_normal_qscore_sum %(max_alt_qual)s
                        --max_alt_alleles_in_normal_count %(max_alt)s
                        --max_alt_allele_in_normal_fraction %(max_fraction)s
                        --tumor_lod %(tumor_LOD)s
                        --input_file:normal %(infile_matched)s ''' % locals()
    if normal_panel:
        statement += ''' --normal_panel %(normal_panel)s ''' % locals()

    if gatk_key:
        statement += " -et NO_ET -K %(gatk_key)s " % locals()

    statement += " > %(mutect_log)s " % locals()

    P.run()

##############################################################################


def strelkaINDELCaller(infile_control, infile_tumor, outfile, genome, config,
                       outdir, job_memory):
    '''Call INDELs using Strelka'''

    statement = '''
    rm -rf %(outdir)s;
    /ifs/apps/bio/strelka-1.0.14/bin/configureStrelkaWorkflow.pl
    --normal=%(infile_control)s  --tumor=%(infile_tumor)s
    --ref=%(genome)s  --config=%(config)s  --output-dir=%(outdir)s;
    checkpoint ; make -j 12 -C %(outdir)s''' % locals()

    P.run()

##############################################################################


def variantAnnotator(vcffile, bamlist, outfile, genome,
                     dbsnp, annotations="", snpeff_file=""):
    '''Annotate variant file using GATK VariantAnnotator'''
    job_options = getGATKOptions()
    job_threads = 3

    if annotations != "":
        anno = annotations.split(",")
        anno = " -A " + " -A ".join(anno)
    else:
        anno = ""
    statement = '''GenomeAnalysisTK -T VariantAnnotator
                    -R %(genome)s
                    -I %(bamlist)s
                    -A SnpEff
                    --snpEffFile %(snpeff_file)s
                    -o %(outfile)s
                    --variant %(vcffile)s
                    -L %(vcffile)s
                    --dbsnp %(dbsnp)s
                    %(anno)s''' % locals()
    P.run()

##############################################################################


def variantRecalibrator(infile, outfile, genome, mode, dbsnp=None,
                        kgenomes=None, hapmap=None, omni=None, mills=None):
    '''Create variant recalibration file'''
    job_options = getGATKOptions()
    job_threads = 3

    track = P.snip(outfile, ".recal")
    if mode == 'SNP':
        statement = '''GenomeAnalysisTK -T VariantRecalibrator
        -R %(genome)s
        -input %(infile)s
        -resource:hapmap,known=false,training=true,truth=true,prior=15.0 %(hapmap)s
        -resource:omni,known=false,training=true,truth=true,prior=12.0 %(omni)s
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 %(dbsnp)s
        -resource:1000G,known=false,training=true,truth=false,prior=10.0 %(kgenomes)s
        -an QD -an SOR -an MQRankSum
        -an ReadPosRankSum -an FS -an MQ
        --maxGaussians 4
        -mode %(mode)s
        -recalFile %(outfile)s
        -tranchesFile %(track)s.tranches
        -rscriptFile %(track)s.plots.R ''' % locals()
        P.run()
    elif mode == 'INDEL':
        statement = '''GenomeAnalysisTK -T VariantRecalibrator
        -R %(genome)s
        -input %(infile)s
        -resource:mills,known=true,training=true,truth=true,prior=12.0 %(mills)s
        -an QD -an MQRankSum
        -an ReadPosRankSum -an FS -an MQ
        --maxGaussians 4
        --minNumBadVariants 5000
        -mode %(mode)s
        -recalFile %(outfile)s
        -tranchesFile %(track)s.tranches
        -rscriptFile %(track)s.plots.R ''' % locals()
        P.run()

##############################################################################


def applyVariantRecalibration(vcf, recal, tranches, outfile, genome, mode):
    '''Perform variant quality score recalibration using GATK '''
    job_options = getGATKOptions()
    job_threads = 3

    statement = '''GenomeAnalysisTK -T ApplyRecalibration
    -R %(genome)s
    -input:VCF %(vcf)s
    -recalFile %(recal)s
    -tranchesFile %(tranches)s
    --ts_filter_level 99.0
    -mode %(mode)s
    -o %(outfile)s ''' % locals()
    P.run()

##############################################################################


def vcfToTable(infile, outfile, genome, columns):
    '''Converts vcf to tab-delimited file'''
    job_options = getGATKOptions()
    job_threads = 3

    statement = '''GenomeAnalysisTK -T VariantsToTable
                   -R %(genome)s
                   -V %(infile)s
                   --showFiltered
                   --allowMissingData
                   %(columns)s
                   -o %(outfile)s''' % locals()
    P.run()

##############################################################################


def selectVariants(infile, outfile, genome, select):
    '''Filter de novo variants based on provided jexl expression'''
    statement = '''GenomeAnalysisTK -T SelectVariants
                    -R %(genome)s
                    --variant %(infile)s
                    -select '%(select)s'
                    -log %(outfile)s.log
                    -o %(outfile)s''' % locals()
    P.run()

##############################################################################


def buildSelectStatementfromPed(filter_type, pedfile, template):
    '''Build a select statement from a template and a pedigree file'''
    pedigree = csv.DictReader(
        IOTools.openFile(pedfile), delimiter='\t', fieldnames=[
            'family', 'sample', 'father', 'mother', 'sex', 'status'])
    affecteds = []
    unaffecteds = []
    parents = []
    select = None
    # loop over pedigree file and establish relationships
    for row in pedigree:
        if row['status'] == '2':
            if filter_type == "denovo":
                father = row['father']
                mother = row['mother']
                proband = row['sample']
            elif filter_type == "dominant" or filter_type == "recessive":
                affecteds += [row['sample']]
            if filter_type == "recessive":
                parents += [row['father'], row['mother']]
        if row['status'] == '1':
            if filter_type == "dominant":
                unaffecteds += [row['sample']]
            elif filter_type == "recessive":
                if row['sample'] not in parents:
                    unaffecteds += [row['sample']]

    # Build select statement from template
    if filter_type == "denovo":
        select = template.replace("father", father)
        select = select.replace("mother", mother)
        select = select.replace("proband", proband)
    elif filter_type == "dominant":
        affecteds_exp = '").getPL().1==0&&vc.getGenotype("'.join(affecteds)
        if len(unaffecteds) == 0:
            unaffecteds_exp = ''
        else:
            unaffecteds_exp = '&&vc.getGenotype("' + \
                ('").isHomRef()&&vc.getGenotype("'.join(unaffecteds)) + \
                '").isHomRef()'
        select = template.replace("affecteds_exp", affecteds_exp)
        select = select.replace("unaffecteds_exp", unaffecteds_exp)
    elif filter_type == "recessive":
        affecteds_exp = '").getPL().2==0&&vc.getGenotype("'.join(affecteds)
        unaffecteds_exp = '").getPL().2!=0&&vc.getGenotype("'.join(unaffecteds)
        if len(parents) == 0:
            parents_exp = ''
        else:
            parents_exp = '&&vc.getGenotype("' + \
                ('").getPL().1==0&&vc.getGenotype("'.join(parents)) + \
                '").getPL().1==0'
        select = template.replace("affecteds_exp", affecteds_exp)
        select = select.replace("unaffecteds_exp", unaffecteds_exp)
        select = select.replace("parents_exp", parents_exp)

    return select

##############################################################################


def guessSex(infile, outfile):
    '''Guess the sex of a sample based on ratio of reads
    per megabase of sequence on X and Y'''
    statement = '''calc `samtools idxstats %(infile)s
                    | grep 'X'
                    | awk '{print $3/($2/1000000)}'`
                    /`samtools idxstats %(infile)s | grep 'Y'
                    | awk '{print $3/($2/1000000)}'`
                    | tr -d " " | tr "=" "\\t" | tr "/" "\\t"
                    > %(outfile)s'''
    P.run()

##############################################################################


def filterMutect(infile, outfile,
                 control_id, tumour_id,
                 min_t_alt, min_n_depth,
                 max_n_alt_freq, min_t_alt_freq):
    ''' filter the mutect snps'''

    def comp(base):
        '''return complementary base'''
        comp_dict = {"C": "G", "G": "C", "A": "T", "T": "A"}
        return comp_dict[base]

    with IOTools.openFile(outfile, "w") as outf:
        with IOTools.openFile(infile, "r") as inf:
            for line in inf.readlines():
                # need to find location of control and tumor columns
                if line.startswith('#CHROM'):
                    columns = line.split("\t")
                    for x in range(0, len(columns)):
                        if control_id in columns[x]:
                            control_col = x
                        elif tumour_id in columns[x]:
                            tumor_col = x

                if not line.startswith('#'):
                    values = line.split("\t")
                    if values[6] == "PASS":
                        t_values = values[tumor_col].split(":")
                        t_ref, t_alt = map(float, (t_values[2].split(",")))
                        t_depth = t_alt + t_ref
                        n_values = values[control_col].split(":")
                        n_ref, n_alt = map(float, (n_values[2].split(",")))
                        n_depth = n_alt + n_ref
                        np.seterr(divide='ignore')

                        # filter
                        if not t_alt > min_t_alt:
                            continue

                        if not n_depth >= min_n_depth:
                            continue

                        if not np.divide(n_alt, n_depth) <= max_n_alt_freq:
                            continue

                        t_freq = np.divide(t_alt, t_depth)
                        n_freq = np.divide(n_alt, n_depth)
                        if (np.divide(t_freq, n_freq) >= min_t_alt_freq or
                            n_freq == 0):
                            outf.write(line)
                else:
                    # write out all comment lines
                    outf.write(line)                            


# the following two functions should be generalised
# currently they operate only on mutect output
def compileMutationalSignature(infiles, outfiles):
    '''takes a list of mutect output files and compiles per sample mutation
    signatures'''

    delim = ":"

    def lookup(b1, b2):
        '''return lookup key for a pair of bases'''
        return(b1 + delim + b2)

    def breakKey(key):
        '''take a lookup key and return the elements'''
        return key.split(delim)

    def comp(base):
        '''return complementary base'''
        comp_dict = {"C": "G", "G": "C", "A": "T", "T": "A"}
        return comp_dict[base]

    def getID(infile):
        return P.snip(os.path.basename(infile),
                      ".mutect.snp.annotated.filtered.vcf")

    outfile1 = IOTools.openFile(outfiles[0], "w")
    mutations = ["C:T", "C:A", "C:G", "A:C", "A:T", "A:G"]

    outfile1.write("%s\t%s\t%s\t%s\t%s\n" % ("patient_id", "base_change",
                                             "ref", "alt", "frequency"))
    patient_freq = {}

    for infile in infiles:
        patient_id = getID(infile)
        mut_dict = {}
        for comb in mutations:
            mut_dict[comb] = 0

        with IOTools.openFile(infile, "r") as f:
            for line in f.readlines():
                if line.startswith('#'):
                    continue
                values = line.split("\t")
                key = lookup(values[3], values[4])
                if key in mut_dict:
                    mut_dict[key] += 1
                else:
                    comp_key = lookup(
                        comp(values[3]), comp(values[4]))
                    mut_dict[comp_key] += 1

        patient_freq[patient_id] = mut_dict

    for mutation in mutations:
        base1, base2 = breakKey(mutation)
        for infile in infiles:
            patient_id = getID(infile)
            outfile1.write("%s\t%s\t%s\t%s\t%s\n" % (patient_id, mutation,
                                                     base1, base2,
                                                     patient_freq[patient_id]
                                                     [mutation]))
    outfile1.close()

    outfile2 = IOTools.openFile(outfiles[1], "w")
    outfile2.write("%s\t%s\n" % ("patient_id",
                                 "\t".join(mutations)))
    for infile in infiles:
        patient_id = getID(infile)
        frequencies = "\t".join(map(str, [patient_freq[patient_id][x]
                                          for x in mutations]))
        outfile2.write("%s\t%s\n" % (patient_id, frequencies))
    outfile2.close()


##############################################################################
# utility functions for pipeline_exome_cancer
##############################################################################

@cluster_runnable
def parseMutectCallStats(infile, outfile):
    '''take the call stats outfile from mutect and summarise the
    reasons for variant rejection'''

    single_dict = collections.defaultdict(int)
    combinations_dict = collections.defaultdict(int)

    with IOTools.openFile(infile, "rb") as infile:
        lines = infile.readlines()
        for i, line in enumerate(lines):
            if i < 2:
                continue
            values = line.strip().split("\t")
            judgement, justification = (values[-1], values[-2])
            if judgement == "REJECT":
                reasons = justification.split(",")
                if len(reasons) == 1:
                    single_dict[reasons[0]] += 1
                else:
                    for reason in reasons:
                        combinations_dict[reasons[0]] += 1

    df = pd.DataFrame([single_dict, combinations_dict])

    df = df.transpose()
    df.columns = ["single", "combination"]
    df = df.sort("single", ascending=False)
    df.index.name = "justification"
    df.to_csv(outfile, header=True, index=True, sep="\t")


@cluster_runnable
def defineEBioStudies(cancer_types, outfile):
    ''' For the cancer types specified in pipeline.ini, identify the
    relevent studies in eBio '''

    cancer_types = cancer_types.split(",")

    cancer_studies_url = "http://www.cbioportal.org/webservice.do?cmd=getCancerStudies"
    genetic_profiles_url = "http://www.cbioportal.org/webservice.do?cmd=getGeneticProfiles"

    type2study_dict = collections.defaultdict(list)
    study2table_dict = collections.defaultdict(list)

    soup = makeSoup(cancer_studies_url)
    for lines in soup.body:
        if isinstance(lines, NavigableString):
            for line in lines.split("\n"):
                if "cancer_study_id" not in line:
                    values = unicode(line).strip().split("\t")
                    if len(values) > 1:
                        cancer_type = values[1].split(" (")[0]
                        if cancer_type in cancer_types:
                            study = re.sub(".\n", "", values[0])
                            type2study_dict[cancer_type].append(study)

    for study_list in type2study_dict.values():
        for study in study_list:
            soup = makeSoup(genetic_profiles_url + "&cancer_study_id=" + study)
            lines = unicode(soup.body).split('\n')
            for line in lines:
                values = line.strip().split("\t")
                if len(values) > 1:
                    if values[1] == "Mutations":
                        genetic_profile = values[0]
                        study2table_dict[study] = genetic_profile

    outf = IOTools.openFile(outfile, "w")

    for cancer_type, study_id in type2study_dict.iteritems():
        for study in study_id:
            table_id = study2table_dict[study]
            outf.write("%s\t%s\t%s\n" % (cancer_type, study, table_id))

    outf.close()


@cluster_runnable
def extractEBioinfo(eBio_ids, vcfs, outfile):
    '''find the number of mutations identitified in previous studies (eBio_ids)
    for the mutated genes in the vcfs'''

    genes = set()

    for vcf in vcfs:
        infile = VCF.VCFFile(IOTools.openFile(vcf))
        for vcf_entry in infile:
            # assumes all vcf entries without "REJECT" are "PASS"
            if vcf_entry.filter != "REJECT":
                info_entries = vcf_entry.info.split(";")
                for entry in info_entries:
                    if "SNPEFF_GENE_NAME" in entry:
                        genes.update((entry.split("=")[1],))

    eBio_ids = IOTools.openFile(eBio_ids, "r")

    tissue_counts = collections.defaultdict(
        lambda: collections.defaultdict(
            lambda: collections.defaultdict(int)))

    def chunks(l, n):
        ''' Yield successive n-sized chunks from l '''
        for i in xrange(0, len(l), n):
            yield l[i:i+n]

    # delete me
    print "number of genes: ", len(list(genes))

    for line in eBio_ids:
        tissue, study, table = line.strip().split("\t")

        n = 0

        for i in xrange(0, len(list(genes)), 500):

            genes_chunk = list(genes)[i:i+500]

            # TS sporadic error when querying with a single gene at a time
            # "urllib2.URLError: <urlopen error [Errno 110] Connection timed out>"
            # max URL length appears to be 8200 characters,
            # try doing 500 genes at a time?

            gene_list = "+".join(list(genes_chunk))

            n += 500

            print "number of genes processed: ", n

            url = ("http://www.cbioportal.org/webservice.do?cmd=getProfileData&"
                   "case_set_id=%(study)s_all&genetic_profile_id=%(table)s&"
                   "gene_list=%(gene_list)s" % locals())

            #with IOTools.openFile(outfile+"_tmp", "w") as out:
            #    out.write("%s\n" % url)

            df = pd.io.parsers.read_csv(url, comment="#", sep="\t", index_col=0)

            # delete me
            print url

            for gene in genes_chunk:
                tmp_df = df[df['COMMON'] == gene]
                # check dataframe contains data!
                if tmp_df.shape[0] != 0:
                    # seem to be having issues with gene set containing duplicates!
                    # --> dataframe with repeated instances of gene after selection
                    # so splice to first row and recreate dataframe from series
                    if tmp_df.shape[0] > 1:
                        tmp_df = pd.DataFrame(tmp_df.iloc[0]).T

                    tissue_counts[tissue][gene]["total"] += tmp_df.shape[1]-2
                    tissue_counts[tissue][gene]["mutations"] += int(tmp_df.count(1))-1

    out = IOTools.openFile(outfile, "w")

    tissues = tissue_counts.keys()

    out.write("gene\t%s\n" % "\t".join([
        "%s_frequency" % x.replace(" ", "_") for x in tissues]))

    for gene in genes:
        freq_values = []
        for tissue in tissues:
            total = tissue_counts[tissue][gene]["total"]
            mutations = tissue_counts[tissue][gene]["mutations"]
            print "total: ", total, "mutations: ", mutations
            freq_values.append(round(np.divide(float(mutations), total), 4))

        out.write("%s\t%s\n" % (gene, "\t".join(map(str, freq_values))))

    out.close()


def intersectionHeatmap(infiles, outfile):
    ''' calculate the intersection between the infiles and plot'''

    pandas2ri.activate()

    name2genes = {}
    df = pd.DataFrame(columns=["id_1", "id_2", "intersection", "perc"])

    ix = 0
    for inf in infiles:

        name = P.snip(os.path.basename(inf)).split(".")[0]
        name = name.replace(".", "_")

        with IOTools.openFile(inf, "r") as f:
            genes = set()

            for line in f:
                if line[0] == "#":
                    continue

                values = line.strip().split("\t")
                info = values[7].split(";")

                for x in info:
                    if x.split("=")[0] == "SNPEFF_GENE_NAME":
                        gene_name = x.split("=")[1]
                        break

                # if no gene name found, line is skipped
                if gene_name:
                    genes.update((gene_name,))

        name2genes[name] = genes
        df.loc[ix] = [name, name, len(genes), 1.0]
        ix += 1

    for pair in itertools.permutations(name2genes.keys(), 2):
        id_1, id_2 = pair
        intersection = len(name2genes[id_1].intersection(name2genes[id_2]))
        not_intersecting = len(name2genes[id_1].symmetric_difference(name2genes[id_2]))
        intersection_perc = float(intersection) / (intersection +
                                                   not_intersecting)

        df.loc[ix] = [id_1, id_2, intersection, intersection_perc]
        ix += 1

    variant = os.path.basename(outfile).replace("overlap_", "").replace("_heatmap.png", "")

    plotIntersectionHeatmap = R('''
    function(df){
    library(ggplot2)
    m_txt = element_text(size=15)
    m_txt_90 = element_text(size=15, angle=90, vjust=0.5, hjust=1)
    l_txt = element_text(size=20)

    p = ggplot(df, aes(id_1, id_2, fill=100*perc)) +
    geom_tile() +
    geom_text(aes(label=intersection), size=3) +
    scale_fill_gradient(name="Intersection (%%)", limits=c(0,100),
                       low="yellow", high="dodgerblue4") +
    theme(axis.text.x = m_txt_90, axis.text.y = m_txt,
          legend.text = m_txt, legend.title = m_txt,
          aspect.ratio=1) +
    xlab("") + ylab("") +
    ggtitle("%(variant)s")

    ggsave("%(outfile)s", width=10, height=10)
    }''' % locals())

    plotIntersectionHeatmap(df)
