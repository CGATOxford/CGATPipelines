import os
import subprocess
import CGAT.Experiment as E


def runSailfishIndex(fasta_file, outdir, threads,
                     kmer):
    '''
    Wrapper for sailfish index
    '''

    if fasta_file.endswith(".fa"):
        pass
    elif fasta_file.endswith(".fasta"):
        pass
    else:
        E.warn("are you sure this is a fasta file?")

    command = '''
    sailfish index --transcripts %s --out %s --threads %i --kmerSize %i
    ''' % (fasta_file, outdir, threads, kmer)

    os.system(command)


def runSailfishQuant(fasta_index, fastq_files, output_dir,
                     paired=False, library="ISF", threads=4,
                     gene_gtf=None):
    '''
    Wrapper for sailfish quant command
    '''

    decompress = False
    if len(fastq_files) > 1:
        if fastq_files[0].endswith(".gz"):
            decompress = True
        else:
            pass
    else:
        if fastq_files[0].endswith(".gz"):
            decompress = True
        else:
            pass

    # check output directory is an absolute path
    if os.path.isabs(output_dir):
        pass
    else:
        out_dir = os.path.abspath(output_dir)

    states = []
    command = " sailfish quant --index %s -l %s  -o %s " % (fasta_index,
                                                            library,
                                                            output_dir)

    states.append(command)

    if threads:
        states.append(" --threads %i " % threads)
    else:
        pass

    if gene_gtf:
        states.append(" --geneMap %s " % gene_gtf)
    else:
        pass

    # sailfish does not handle compress files natively,
    # need to decompress on the fly with advanced
    # bash syntax
    if decompress and paired:
        first_mates = [fq for fq in fastq_files if re.search("fastq.1.gz",
                                                             fq)]
        fstr_format = " ".join(["%s" for hq in first_mates])
        fdecomp_format = fstr_format % first_mates
        decomp_first = " -1 <(zcat %s)" % fdecomp_format

        states.append(decomp_first)

        second_mates = [sq for sq in fastq_files if re.search("fastq.2.gz",
                                                              sq)]
        sstr_format = " ".join(["%s" for aq in second_mates])
        sdecomp_format = sstr_format % second_mates
        decomp_second = " -2 <(zcat %s)" % sdecomp_format

        states.append(decomp_second)

    elif decompress and not paired:
        first_mates = [fq for fq in fastq_files if re.search("fastq.1.gz",
                                                             fq)]
        fstr_format = " ".join(["%s" for sq in first_mates])
        fdecomp_format = fstr_format % first_mates
        decomp_first = " -r <(zcat %s)" % fdecomp_format

        states.append(decomp_first)

    elif paired and not decompress:
        first_mates = [fq for fq in fastq_files if re.search("fastq.1",
                                                             fq)]
        fstr_format = " ".join(["%s" for sq in first_mates])
        fdecomp_format = fstr_format % first_mates
        decomp_first = " -1 %s " % fdecomp_format

        states.append(decomp_first)

        second_mates = [sq for sq in fastq_files if re.search("fastq.2",
                                                              sq)]
        sstr_format = " ".join(["%s" for aq in second_mates])
        sdecomp_format = sstr_format % second_mates
        decomp_second = " -2 %s " % sdecomp_format

        states.append(decomp_second)

    statement = " ".join(states)

    # subprocess cannot handle process substitution
    # therefore needs to be wrapped in /bin/bash -c '...'
    # for bash to interpret the substitution correctly

    process = subprocess.Popen(statement, shell=True,
                               executable="/bin/bash")

    stdout, stderr = process.communicate()

    if process.returncode != 0:
        raise OSError(
            "-------------------------------------------\n"
            "Child was terminated by signal %i: \n"
            "The stderr was \n%s\n%s\n"
            "-------------------------------------------" %
            (-process.returncode, stderr, statement))


def runKallistoIndex(fasta_file, outfile, kmer=31):
    '''
    Wrapper for kallisto index
    '''

    if fasta_file.endswith(".fa"):
        pass
    elif fast_file.endswith(".fasta"):
        pass
    else:
        E.warn("are you sure this is a fasta file?")

    command = "kallisto index --index=%s  %s" % (outfile,
                                                 fasta_file)

    os.system(command)


def runKallistoQuant(fasta_index, fastq_files, output_dir,
                     bias=False, bootstrap=None,
                     seed=1245, threads=None, plaintext=False):
    '''
    Wrapper for kallisto quant command
    '''

    if len(fastq_files) > 1:
        fastqs = " ".join(fastq_files)
    else:
        fastqs = fastq_files

    # check output directory is an absolute path
    if os.path.isabs(output_dir):
        pass
    else:
        out_dir = os.path.abspath(output_dir)

    states = []
    command = " kallisto quant --index=%s --output-dir=%s" % (fasta_index,
                                                              output_dir)
    states.append(command)

    if bias:
        states.append(" --use-bias ")
    else:
        pass

    if bootstrap:
        states.append(" --bootstrap=%i --seed=%i " % (bootstrap,
                                                      seed))
    else:
        pass

    if plaintext:
        states.append(" --plaintext ")
    else:
        pass

    if threads:
        states.append(" --threads=%i " % threads)
    else:
        pass

    states.append(" %s " % fastqs)

    statement = " ".join(states)

    # need to rename output files to conform to input/output
    # pattern as required.  Default name is abundance*.txt
    # when using plaintext output
    # kaliisto requires an output directory - create many small
    # directories, one for each file.
    # then extract the abundance.txt file and rename using the
    # input/output pattern

    os.system(statement)
