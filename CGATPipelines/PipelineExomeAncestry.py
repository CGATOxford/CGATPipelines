import pandas as pd
import numpy as np
import CGAT.IOTools as IOTools
import random
import itertools
from CGATPipelines.Pipeline import cluster_runnable
import json
import decimal
from sklearn.cluster import KMeans
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches


@cluster_runnable
def MakeSNPFreqDict(infiles, outfiles, rs):
    '''
    Generates a random set of SNPs to use to characterise the ancestry
    and relatedness of the samples.

    1. Finds all SNPs which have a known genotype frequency for all 11
       hapmap ancestries
    2. Picks a random 50000 of these SNPs
    3. Stores a list of these SNPs as randomsnps.tsv
    4. Builds a dictionary of dictionaries of dictionaries where:
           At Level 1 each key is a HapMap ancestry ID
               each of these keys leads to another dictionary - level 2.
           At Level 2 each key is a dbSNP SNP ID
               each of these keys leads to another dictionary - level 3.
           At Level 3 each key is a genotype:
               the values of these keys are the genotype frequency of this
               genotype at this SNP in this ancestry.
       e.g. snpdict['ASW']['rs000000001']['CT'] -->
            frequency of the CT genotype at the rs0000001 SNP in the
            ASW population
    5.  Stores this dictionary in json format as randomsnps.json
    '''

    # list all chromosomes
    chroms = np.array([f.split("/")[-1].split("_")[2] for f in infiles])
    uchroms = set(chroms)
    infiles = np.array(infiles)
    snpdict = dict()
    for chrom in uchroms:
        # get all the files for this chromosome
        c = np.where(chroms == chrom)[0]
        thischrom = infiles[c]
        snpidsets = []
        # make a set of all the snp ids of known frequency for
        # this chromsome for each ancestry
        for f in thischrom:
            snpids = []
            anc = f.split("/")[-1].split("_")[3]
            with IOTools.openFile(f) as inp:
                for line in inp:
                    line = line.strip().split(" ")
                    snpids.append(line[0])
            snpids = set(snpids)
            snpidsets.append(snpids)
        # find the snp ids where frequency is known in all ancestries
        thischromsnps = set.intersection(*snpidsets)
        snpdict[chrom] = thischromsnps

    vals = list(snpdict.values())
    # set of all the SNPs genotyped in all ancestries on all chroms
    pooled = set.union(*vals)
    # take a random sample of snps from this set
    random.seed(rs)
    sam = set(random.sample(pooled, 50000))

    # make a new dictionary of just the sampled snps for each chrom
    snpdictsam = dict()
    for key in snpdict:
        i = snpdict[key] & sam
        snpdictsam[key] = i

    # generate all possible genotypes as list of strings
    # (e.g. CC CT CG CA TT etc)
    res = ["%s%s" % (g[0], g[1])
           for g in itertools.permutations("CGAT", 2)]
    res += ['CC', 'GG', 'TT', 'AA']
    freqdict = dict()
    ancs = set([l.split("/")[-1].split("_")[3] for l in infiles])

    # build an empty dictionary of dictionaries of dictionaries
    # to store genotype freqs for each ancestry
    for anc in ancs:
        freqdict.setdefault(anc, dict())
        for s in sam:
            freqdict[anc].setdefault(s, dict())
            for r in res:
                freqdict[anc][s][r] = 0.0

    # read genotype freqs from the input file and store in the dictionary
    for f in infiles:
        bits = f.split("/")[-1].split("_")
        anc = bits[3]
        chrom = bits[2]
        snpset = snpdictsam[chrom]
        with IOTools.openFile(f) as input:
            for line in input:
                line = line.strip().split(" ")
                if line[0] in snpset:
                    G1 = "%s%s" % (line[10][0], line[10][-1])
                    try:
                        F1 = float(line[11])
                        G2 = "%s%s" % (line[13][0], line[13][-1])
                        F2 = float(line[14])
                        G3 = "%s%s" % (line[16][0], line[16][-1])
                        F3 = float(line[17])
                    except:
                        # because occasionally the GFs are not numeric
                        continue
                    freqdict[anc][line[0]][G1] = F1
                    freqdict[anc][line[0]][G2] = F2
                    freqdict[anc][line[0]][G3] = F3

    # dump dictionary in a json file to resurrect later
    outf = IOTools.openFile(outfiles[0], "w")
    json.dump(freqdict, outf)
    outf.close()

    # store the sampled list of snps
    out = IOTools.openFile(outfiles[1], "w")
    for snp in sam:
        out.write("%s\n" % snp)
    out.close()


@cluster_runnable
def GenotypeSNPs(infile, snplist, outfile):
    '''
    Fetches the genotype from the variant tables for all samples
    for SNPs in the hapmap sample from makeRandomSNPSet.

    Complex sites are ignored (as simple SNPs are sufficient for these
    calculations).
    These are:
        Sites which failed QC (column 3 in the variant table is not PASS)
        Sites with more than 2 alleles defined (column 6 in the variant table
        contains more than one alternative allele)
        SNPs with more than one ID
        Indels
    '''
    out = IOTools.openFile(outfile, "w")
    with IOTools.openFile(infile) as inf:
        for line in inf:
            line = line.strip().split()
            # if the variant passed QC
            if line[4] == "PASS":
                genotype = line[7]
                # if the genotype looks normal e.g. 1/1
                if len(genotype) == 3:
                    # get the actual genotype (rather than the index)
                    if genotype[0] != ".":
                        ind1 = int(genotype[0])
                    else:
                        ind1 = 0
                    if genotype[2] != ".":
                        ind2 = int(genotype[2])
                    else:
                        ind2 = 0
                    A1 = line[5]
                    A2 = line[6].split(",")
                    AS = [A1] + A2

                    if len(AS) <= 2:
                        GT = "%s%s" % (AS[ind1], AS[ind2])
                        refGT = "%s%s" % (A1, A1)
                        if len(GT) == 2:
                            if line[3][0:2] == "rs" and len(
                                    line[3].split(";")) == 1:
                                snpid = line[3]
                                chrom = line[0]
                                pos = line[1]
                                if snpid in snplist:
                                    out.write("%s\t%s\t%s\t%s\t%s\n"
                                              % (snpid, chrom, pos, GT,
                                                 refGT))
    out.close()


@cluster_runnable
def CalculateAncestry(infile, calledsnps, snpdict, outfiles):
    '''
    Takes the data stored in MakeRandomSNPSet and the genotype of each sample
    at each site in calledsnps.tsv and tabulates the frequency of this
    genotype in each of the HapMap ancestry categories.
    The overall probability of each ancestry is then calculated as the
    product of these frequencies. These can only be used in comparison to
    each other - to show which of the 11 ancestries is most probable.
    '''

    # List the SNPs in the sample where a variant has been called
    # Record the reference genotype at each of these SNPs
    called = set()
    refs = dict()
    for line in IOTools.openFile(calledsnps).readlines():
        line = line.strip().split("\t")
        refs[line[0]] = line[1]
        called.add(line[0])

    # Record the genotype of the current sample at each SNP
    currents = dict()
    for line in IOTools.openFile(infile).readlines():
        line = line.strip().split("\t")
        currents[line[0]] = line[3]

    # Regenerate the genotype frequency dictionary made earlier
    s = IOTools.openFile(snpdict)
    snpdict = json.load(s)
    s.close()
    ancs = list(snpdict.keys())

    # Build a table of the frequency of the genotype in this sample
    # in each of the hapmap ancestries
    out = IOTools.openFile(outfiles[0], "w")
    out.write("snp\tgenotype\t%s\n" % "\t".join(ancs))
    for snp in called:
        # where a variant hasn't been called in the current sample assume
        # the reference genotype
        if snp in currents:
            geno = currents[snp]
        else:
            geno = refs[snp]

        res = []
        for anc in ancs:
            D = snpdict[anc]
            sub = D[snp]
            freq = sub[geno]
            res.append(freq)

        res = [str(r) for r in res]
        out.write("%s\t%s\t%s\n" % (snp, geno, "\t".join(res)))

    out.close()
    out = IOTools.openFile(outfiles[1], "w")

    # read the output of the previous step into pandas
    res = pd.read_csv(outfiles[0], sep="\t")

    # calculate the probability of each ancestry as the product of all the
    # genotype frequencies
    # SNPs where a genotype not recorded in hapmap has been called are
    # skipped
    # decimal package is used because the results are very small floats
    res = res.replace(0, float('nan'))
    for anc in set(ancs):
        r = res[anc][np.invert(res[anc].isnull())]
        r = [decimal.Decimal(d) for d in r]
        out.write("%s\t%s\n" % (anc, np.product(r)))
    out.close()


@cluster_runnable
def MakePEDFile(infiles, outfiles):
    '''
    Generates the required input for the Plink and King software packages.
       - PED file - columns are SNPs and rows are samples, each cell is the
         genotype of the sample at this SNP
       - MAP file - rows are SNPs in the same order as the columns in the ped
         file, each row shows chromosome, snp id and position.
    '''
    # Record the chromosome and position of each SNP
    chromposdict = dict()
    for line in IOTools.openFile(infiles[-1]).readlines():
        line = line.strip().split("\t")
        chromposdict[line[0]] = ((line[2], line[3]))

    # Generate the PED file - takes the genotype column of the table from
    # the CalculateAncestry step and transposes it across a row

    out = IOTools.openFile(outfiles[0], "w")
    for inf in infiles[0:-1]:
        inf = inf[0]
        df = pd.read_csv(inf, sep="\t")
        geno = df['genotype']
        genolist = []
        for g in geno:
            genolist.append(" ".join(list(g)))
        out.write("%s\t%s\n" % (inf.split("/")[-1], "\t".join(genolist)))
    out.close()

    # Generate the MAP file - chromosome and position of each SNP in the same
    # order as the PED file
    mapf = IOTools.openFile(outfiles[1], "w")
    slist = df['snp']
    for snp in slist:
        mapf.write("%s\t%s\t%s\n" % (chromposdict[snp][0], snp,
                                     chromposdict[snp][1]))
    mapf.close()


def CalculateFamily(infile, outfile):
    '''
    Translates and filters the output from King.
    Pairs of related samples are written to the output file along with the
    degree of relatedness.  Degrees are decided using thresholds from
    the King documentation, here
    http://people.virginia.edu/~wc9c/KING/manual.html
    '''
    out = IOTools.openFile(outfile, "w")
    input = open(infile).readlines()[1:]
    for line in input:
        line = line.strip().split("\t")
        score = float(line[14])
        # if there is a significant relatedness score
        if score >= 0.0442:
            ID1 = line[0]
            ID2 = line[2]
            # 2nd degree relatives
            if score <= 0.177:
                degree = 2
            # 1st degree relatives
            elif score <= 0.354:
                degree = 1
            # 0th degree relatives (ID twins or duplicate)
            else:
                degree = 0
            out.write("%s\t%s\t%s\n" % (ID1, ID2, degree))
    out.close()


@cluster_runnable
def CalculateSex(infiles, outfile):
    '''
    Approximates the sex of the sample using the data in the variant table.
    Basic estimate based on heterozygosity on the X chromosome - genotypes
    are taken for all called variants on X passing QC and the percentage
    of heterozygotes is taken.
    This tends to produce two clear populations so Kmeans clustering
    is used to split the data into two - male and female.  Samples which are
    unclear are marked in the output.
    '''
    percs = []
    for f in infiles:
        homs = 0
        tot = 0
        # Count the homozygous and total variants called on the
        # X chromosome in each input file (after QC)
        with IOTools.openFile(f) as inf:
            for line in inf:
                line = line.strip().split("\t")
                if line[0] == "chrX" and line[4] == "PASS":
                    geno = line[7]
                    s = set(list(geno))
                    if len(s) == 2:
                        homs += 1
                    tot += 1
        perc = float(homs) / float(tot)
        percs.append(perc)

    # impose a grid on the data so clustering is possible
    ests = np.array(percs)
    ones = np.array([1] * len(ests))
    grid = np.column_stack((ones, ests))

    # divide data into two clusters
    means = KMeans(n_clusters=2)
    fit1 = means.fit(grid)
    res = fit1.predict(grid)

    # find the position of the cluster centres on the y axis
    ccs = fit1.cluster_centers_
    cc1 = ccs[0, 1]
    cc2 = ccs[1, 1]
    ccs = [cc1, cc2]

    # pick out the scores from each cluster and calculate the standard dev
    s1 = ests[res == 0]
    s2 = ests[res == 1]
    sds = [np.std(s1), np.std(s2)]
    # males will be more homozygous - the bigger value represents males
    if cc1 > cc2:
        sexes = ["male", "female"]
    else:
        sexes = ["female", "male"]

    x = 0
    out = IOTools.openFile(outfile, "w")
    for f in infiles:
        # for each input file take the group, the score, the standard deviation
        # of this group, the distance from the cluster centre and the sex
        group = res[x]
        est = ests[x]
        std = sds[group]
        dist = est - ccs[group]
        sex = sexes[group]

        # is the value with 3 SDs of the mean for this sex
        if sex == "male":
            if ((3 * std) - abs(dist) >= 0 or dist >= 0):
                sig = "*"
            else:
                sig = "-"
        elif sex == "unknown":
            sig = "-"
        else:
            if ((3 * std) - abs(dist) >= 0 or dist <= 0):
                sig = "*"
            else:
                sig = "-"
        out.write("%s\t%f\t%f\t%s\t%s\n" % (f, est, dist, sexes[group], sig))
        x += 1
    out.close()

    # plot the output - blue dots male, red dots female
    # lines represent 3 standard deviations from the cluster centre
    p = plt.figure()
    a = p.add_subplot('111')

    a.axis([0, len(infiles),
            (min([min(s1), min(s2)]) - 0.1), (max([max(s1), max(s2)]) + 0.1)])

    dist1 = np.linspace(0, len(infiles), len(s1))
    dist2 = np.linspace(0, len(infiles), len(s2))

    a.plot(dist1, s1, 'b.')
    a.plot(dist2, s2, 'r.')

    top = ccs[0] + (3 * sds[0])
    bottom = ccs[1] - (3 * sds[1])
    a.plot([0, len(infiles)], [top, top], 'b-')
    a.plot([0, len(infiles)], [bottom, bottom], 'r-')
    a.tick_params(axis='x', labelbottom='off', bottom='off')
    a.yaxis.set_label_text('percentage homozygous X variants')
    fig = outfile.replace(".tsv", ".png")
    p.savefig(fig)


def PlotAncestry(infile):
    '''
    Draws a plot showing the score (probabilty) for each individual
    for their assigned ancestry and the second closest match.
    x = individual
    y = score
    large diamonds represent the best match and small triangles the second
    best match
    '''
    ancestry = pd.read_csv(infile, sep="\t", header=None, dtype="str")

    # score using the index of the probability e.g. 5e-1000 as -1000
    ancestry[2] = [int(a[1]) * -1 for a in ancestry[2].str.split("-")]
    ancestry[4] = [int(a[1]) * -1 for a in ancestry[4].str.split("-")]

    # sort by assigned ancestry
    ancestry = ancestry.sort_values([1])
    ancs = list(set(ancestry[1].append(ancestry[3])))

    # assign colours to points by ancestry
    ancestry[5] = ancestry[1].apply(lambda x: ancs.index(x), 1)
    ancestry[6] = ancestry[3].apply(lambda x: ancs.index(x), 1)
    floats = np.linspace(0, 1, len(set(ancs)))
    colours = cm.jet(floats)

    f = plt.figure(figsize=(30, 15))
    a = f.add_subplot('111')
    ymin = min(ancestry[2]) - 100
    ymax = max(ancestry[4]) + 100

    a.axis([-5, len(ancestry) + 5, ymin, ymax])

    # plot best scoring ancestry
    a.scatter(list(range(len(ancestry))),
              ancestry[2], c=colours[ancestry[5].values],
              s=500, marker='D', edgecolors='None')

    # plot second best scoring ancestry
    a.scatter(list(range(len(ancestry))),
              ancestry[4], c=colours[ancestry[6].values],
              s=150, marker="^", edgecolors='None')
    sancs = set(ancestry[1])
    i = 20

    # draw lines between individuals assigned to different ancestries
    for anc in sancs:
        p = max(np.where(ancestry[1] == anc)[0]) + 0.5
        p2 = min(np.where(ancestry[1] == anc)[0]) + 0.5
        a.plot([p, p],
               [ymin, ymax], 'k-')
        a.text(p2, ymax - i, anc, ha='center')
        i += 10

    # legend
    ps = [mpatches.Patch(color=col) for col in colours]
    a.legend(handles=ps, labels=ancs, loc=3)

    f.savefig(infile.replace(".tsv", ".png"))
