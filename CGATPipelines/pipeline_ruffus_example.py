'''
This is a simple pipeline to demonstrate how ruffus is used to control
the flow of data through a series of functions.

There is a PowerPoint presentation with details of all the ruffus decorators
and which goes through the functions in this pipeline one by one at
https://www.cgat.org/downloads/public/training/CGATTrainingSessions/2017-02-14-KatyBrownRuffus.pptx

No input is required, to run the pipeline the steps are:

1. Create an empty folder in your "projects" directory

2. Log on to cgath1

3. Inside the empty folder, run:
    (replacing CGATPipelines with the path to your CGATPipelines directory)

    python CGATPipelines/pipeline_ruffus_example.py plot full
    to see an image of the structure of the pipeline

    python CGATPipelines/pipeline_ruffus_example.py show full -v 5
    to see what the input and output files will be

    python CGATPipelines/pipeline_ruffus_example.py make full -v 5
    to run the entire pipeline

    python CGATPipelines/pipeline_ruffus_example.py make FUNCTIONNAME -v 5
    where FUNCTIONNAME is the name of a function (e.g. exampleOriginate)
    to run a specific function.

full can also be replaced with the name of a function for the show and plot
examples above.

'''
from ruffus import *
from ruffus.combinatorics import *
import sys
import os
import CGATPipelines.Pipeline as P
import shutil


PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])


@originate(['a_originate.tsv', 'b_originate.tsv'])
def exampleOriginate(outfile):
    '''
    Example of ruffus originate decorator.

    @originate generates new files from scratch (0 to many operation)

    Here, this function generates two new files:
    a_originate.tsv
    b_originate.tsv

    '''
    out = open(outfile, "w")
    out.write("I am %s - generated using @originate" % outfile)
    out.close()


@transform(exampleOriginate,
           suffix("_originate.tsv"), "_transform.new")
def exampleTransform(infile, outfile):
    '''
    Example of the ruffus transform decorator.

    @transform runs a function on a single file and generates a single output
    (1 to 1 operation).

    Here this function transforms the output of exampleOriginate into two new
    files:
    a_originate.tsv
    -->
    a_transform.new
    b_originate.tsv
    -->
    b_transform.new
    '''
    content = open(infile).readline()
    content = content.replace("originate", "transform")
    content = content.replace("tsv", "new")
    out = open(outfile, "w")
    out.write(content)
    out.close()


@merge(exampleTransform, "summary_merge.tsv")
def exampleMerge(infiles, outfile):
    '''
    Example of the ruffus merge decorator.

    @merge operates on several files to give a single output (many to 1
    operation)

    Here, this function summarises the output from exampleTransform into a
    single file.

    [a_transform.new, b_transform.new] --> summary_merge.tsv

    '''
    out = open(outfile, "w")
    out.write("""I am summary_merge.tsv - generated using @merge\n""")
    out.write("""My input files contained this text:\n\n""")

    for infile in infiles:
        content = open(infile).readline()
        out.write("%s\n" % content)
    out.close()


@subdivide(exampleTransform,
           regex("(.*)_transform.new"),
           [r"\1_subdivide_1.tsv",
            r"\1_subdivide_2.tsv"])
def exampleSubdivide(infile, outfiles):
    '''
    Example of the ruffus subdivide decorator

    @subdivide takes a set of input files and generates mulitple output files
    for each (fewer to many operation).

    Here, this function splits each output from exampleTransform into
    two output files.

    a_transform.new
    -->
    a_subdivide_1.tsv
    a_subdivide_2.tsv

    b_transform.new
    -->
    b_subdivide_1.tsv
    b_subdivide_2.tsv

    '''
    content = open(infile).readline()
    content1, content2 = content.split("-")
    out1 = open(outfiles[0], "w")
    out1.write("I am %s generated using @subdivide\n\n" % outfiles[0])
    out2 = open(outfiles[1], "w")
    out2.write("I am %s generated using @subdivide\n\n" % outfiles[1])

    out1.write("The first half of %s was:   %s" % (infile, content1))
    out2.write("The second half of %s was:   %s" % (infile, content2))
    out1.close()
    out2.close()


@collate(exampleSubdivide,
         regex("(.*)_(.*)_(.*).tsv"),
         r"\1_collate.tsv")
def exampleCollate(infiles, outfile):
    '''
    Example of the ruffus collate decorator

    @collate uses a regular expression to subset the input files (many to fewer
    operation)

    Here, this function groups the output of exampleSubdivide based on the
    first letter of the filename (a or b)

    a_subdivide_1.tsv
    a_subdivide_2.tsv
    -->
    a_collate.tsv

    b_subdivide_1.tsv
    b_subdivide_2.tsv
    -->
    b_collate.tsv
    '''
    out = open(outfile, "w")
    out.write("I am %s generated using @collate\n\n" % outfile)
    for infile in infiles:
        content = open(infile).readline()
        out.write("The first line of input file %s was:  %s\n" % (infile,
                                                                  content))
    out.close()


@split(exampleMerge, "split*.oneline")
def exampleSplit(infile, outfiles):
    '''
    Example of the ruffus split decorator.

    @split splits a single file into an unspecified number of output files
    (1 to many operation)

    Here, this function splits the output of exampleMerge into one file for
    each line in the input named split_x.oneline
    summary_merge.tsv
    -->
    split_1.oneline
    split_2.oneline
    split_3.oneline
    split_4.oneline
    split_5.oneline
    '''
    content = open(infile).readlines()
    i = 1
    for line in content:
        out = open("split_%i.oneline" % i, "w")
        out.write("""I am split_%i.oneline, generated using @split.\n\n""" % i)
        out.write("""Line %i of the input file was:  %s""" % (i, line.strip()))
        out.close()
        i += 1


@follows(exampleSplit, exampleCollate, exampleMerge)
def basicRuffus():
    '''
    This is an example of a dummy task to run a subsection of a pipeline.
    Running the pipeline as make basicRuffus will update, if needed,
    exampleSplit, exampleCollate, exampleMerge, exampleTransform
    and exampleOriginate.

    exampleOriginate, exampleTransform and exampleSubdivide
    do not need to be specified after
    @follows because they will run automatically as the prerequisites to
    exampleSplit, exampleCollate and exampleMerge
    '''
    pass


@follows(mkdir("combinations"))
@combinations(exampleSubdivide,
              formatter(), 2,
              "combinations/{basename[0][0]}_{basename[1][0]}.tsv")
def exampleCombinations(infiles, outfile):
    '''
    Example of the ruffus combinations decorator.

    @combinations generates every possible combination of length n of
    a list of files, where order is not important.  n is specified after
    the formatter in the decorator (n here is 2)

    Here, every combination of the four subdivide output files is generated.

    a_subdivide_1.tsv
    a_subdivide_2.tsv
    b_subdivide_1.tsv
    b_subdivide_2.tsv
    -->
    combinations/a_subdivide_1_a_subdivide_2.tsv
    combinations/a_subdivide_1_b_subdivide_1.tsv
    combinations/a_subdivide_1_b_subdivide_2.tsv
    combinations/a_subdivide_2_b_subdivide_1.tsv
    combinations/a_subdivide_2_b_subdivide_2.tsv
    combinations/b_subdivide_1_b_subdivide_2.tsv
    '''
    out = open(outfile, "w")
    out.write("I am %s, generated from %s and %s" % (outfile, infiles[0],
                                                     infiles[1]))
    out.close()


@follows(mkdir("permutations"))
@permutations(exampleSubdivide,
              formatter(), 2,
              "permutations/{basename[0][0]}_{basename[1][0]}.tsv")
def examplePermutations(infiles, outfile):
    '''
    Example of the ruffus permutations decorator.

    @permutations generates every possible combination of length n of
    a list of files, where order is important.  n is specified after
    the formatter in the decorator (n here is 2)

    Here, every permutation of the four subdivide output files is generated.

    a_subdivide_1.tsv
    a_subdivide_2.tsv
    b_subdivide_1.tsv
    b_subdivide_2.tsv
    -->
    permutations/a_subdivide_1_a_subdivide_2.tsv
    permutations/a_subdivide_1_b_subdivide_1.tsv
    permutations/a_subdivide_1_b_subdivide_2.tsv
    permutations/a_subdivide_2_a_subdivide_1.tsv
    permutations/a_subdivide_2_b_subdivide_1.tsv
    permutations/a_subdivide_2_b_subdivide_2.tsv
    permutations/b_subdivide_1_a_subdivide_1.tsv
    permutations/b_subdivide_1_a_subdivide_2.tsv
    permutations/b_subdivide_1_b_subdivide_2.tsv
    permutations/ b_subdivide_2_a_subdivide_1.tsv
    permutations/b_subdivide_2_a_subdivide_2.tsv
    permutations/b_subdivide_2_b_subdivide_1.tsv
    '''
    out = open(outfile, "w")
    out.write("I am %s, generated from %s and %s" % (outfile, infiles[0],
                                                     infiles[1]))

    out.close()


@follows(mkdir("product"))
@product(exampleTransform,
         formatter(),
         exampleSubdivide,
         formatter(),
         "product/{basename[0][0]}_{basename[1][0]}.tsv")
def exampleProduct(infiles, outfile):
    '''
    Example of the ruffus product decorator.

    @product generates every possible pair from two lists of files, where
    order is not important.

    Here, every possible pairing of the two output files from exampleTransform
    and the four output files from exampleSubdivide is generated.

    a_subdivide_1.tsv
    a_subdivide_2.tsv
    b_subdivide_1.tsv
    b_subdivide_2.tsv

    and

    a_transform.new
    b_transform.new

    -->
    product/a_transform_a_subdivide_1.tsv
    product/a_transform_a_subdivide_2.tsv
    product/a_transform_b_subdivide_1.tsv
    product/a_transform_b_subdivide_2.tsv
    product/b_transform_a_subdivide_1.tsv
    product/b_transform_a_subdivide_2.tsv
    product/b_transform_b_subdivide_1.tsv
    product/b_transform_b_subdivide_2.tsv

    '''
    out = open(outfile, "w")
    out.write("I am %s, generated from %s and %s" % (outfile, infiles[0],
                                                     infiles[1]))

    out.close()


@follows(exampleCombinations, examplePermutations, exampleProduct)
def advancedRuffus():
    '''
    This is a dummy function to demonstrate the use of dummy functions to run
    subsections of the pipeline.
    Running the pipeline as make advancedRuffus will update, if needed,
    exampleCombinations, examplePermutations and exampleProduct,
    plus any prior steps they depend upon - these are exampleOriginate,
    exampleTransform and exampleSubdivide.
    exampleMerge, exampleSplit and exampleCollate will not be run.
    '''
    pass


@follows(basicRuffus, advancedRuffus)
def full():
    '''
    All CGAT pipelines should end with a full() function which updates,
    if needed, all branches of the pipeline.
    The @follows statement should ensure that all functions are covered,
    either directly or as prerequisites.
    '''
    pass


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
