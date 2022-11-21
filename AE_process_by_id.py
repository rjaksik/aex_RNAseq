"""-------------------------------------------------------------------------------------------------
This script downloads and processes all samples belonging to an experiment with given ID
in ArrayExpress database. Steps are performed as follows:
1) all samples associated with given exp ID are searched in ArrayExpress database
2) TrimGalore is used to remove adapter sequences
3) reads are mapped to the genome given in -starGenome parameter
4) samtools quickcheck is performed on resulting BAM file
if BAM file if found to be fine following analysing tools are run:
5) samotools flagstat
6) bedtools coverage
7) FastQC
8) read_GC.py, tin.py, inner_distance.py, geneBody_coverage.py from RSeQC package

DEPENDENCIES:
Following programs must be avaiable on your PATH:
FastQC, TrimGalore-0.5.0, STAR, samtools_1.3, bedtools2
And following python packages must be installed:
cutadapt, RSeQC

Script might be run both using Python2 and Python3.
-------------------------------------------------------------------------------------------------"""

from __future__ import absolute_import
from __future__ import print_function
from optparse import OptionParser
try:
    from urllib.parse import urlparse, urlencode
    from urllib.request import urlopen, Request
    from urllib.error import HTTPError, URLError
except ImportError:
    from urlparse import urlparse
    from urllib import urlencode
    from urllib2 import urlopen, Request, HTTPError, URLError
import ssl
import os
import sys
import subprocess

try:
    import itertools.ifilter as filter
except ImportError:
    pass
from io import open

__author__ = "Pawel Kus & Roman Jaksik"
__email__ = "kpawel2210@gmail.com; roman.jaksik@polsl.pl"

COVERAGE = "coverage"
FASTQC = "fastqc"
FLAGSTAT = "flagstat"
GENEBODY_COVERAGE = "geneBody_coverage"
INNER_DISTANCE = "inner_distance"
READ_GC = "read_gc"
TIN = "tin"

files = {
    COVERAGE: ["coverage"],
    FASTQC: ["fastqc.zip"],
    FLAGSTAT: ["flagstat"],
    GENEBODY_COVERAGE: ["geneBodyCoverage.txt"],
    INNER_DISTANCE: ["inner_distance_freq.txt", "inner_distance.txt"],
    READ_GC: ["read_GC.GC.xls"],
    TIN: ["tin.summary", "tin.xls"]
}


def log(*arg):
    length = 100
    for a in arg:
        n = int(round((length - len('%s' % a)) / 2))
        a = '-' * n + '%s' % a
        a = a + '-' * (length - len(a))
        os.system("echo '%s'" % a)


def get_and_process(exp_id, genome_path, bed_file, simulate, simulate_all):

    def run_task(cmd):
        if simulate:
            os.system("echo 'command: %s'" % cmd)
        else:
            os.system(cmd)

    def check_if_results_exist(id):
        # get all non empty files
        ls = os.listdir("./")
        res = []
        for file in ls:
            if os.stat(file).st_size > 0:
                res.append(file)

        done = {}
        all_exist = True
        any_exist = False
        for (key, val) in files.items():
            if simulate_all:
                done[key] = False
                all_exist = False
                any_exist = False
            else:
                filenames = ["%s.%s" % (id, v) for v in val]
                exist = [file in res for file in filenames]
                done[key] = all(exist)
                if not all(exist):
                    all_exist = False
                if any(exist):
                    any_exist = True
        return all_exist, any_exist, done

    def adapters_trimmed(files):
        final = ["_val_" in fname for fname in files]
        temp = ["_trimmed_" in fname for fname in files]
        if simulate_all:
            return False
        elif any(final) and not any(temp):
            return True
        else:
            return False

    def verify_bam_file(bam_file):
        if simulate_all:
            os.system("echo 'command: samtools quickcheck %s'" % bam_file)
            bam_ok = False
        elif not os.path.exists(bam_file):
            bam_ok = False
        else:
            try:
                subprocess.check_output("samtools quickcheck %s" % bam_file, shell=True)
                bam_ok = True
            except subprocess.CalledProcessError as e:
                log("Error occured", e)
                bam_ok = False
        return bam_ok

    # create directories
    exp_dir = "./" + exp_id
    os.system("mkdir -p -v %s" % exp_dir)
    os.chdir(exp_dir)

    # download experiment info
    url = "https://www.ebi.ac.uk/arrayexpress/files/%s/%s.sdrf.txt" % (exp_id, exp_id)
    os.system("wget --no-check-certificate " + url)

    # find samples IDs
    file = open("%s.sdrf.txt" % exp_id)
    header = file.readline().split("\t")
    name = list(filter(lambda x: "[ENA_RUN]" in x, header))
    num = header.index(name[0])
    samples_ids = [line.split("\t")[num] for line in file]
    samples_ids = list(set(samples_ids))
    samples_ids.sort()
    os.system("rm %s.sdrf.txt" % exp_id)

    # save samples table
    f = open("%s.samples" % exp_id, "w")
    for id in samples_ids:
        f.write("%s\t%s\n" % (exp_id, id))

    # processing samples
    broken_bams = []
    for id in samples_ids:

        fastq_dir = id + "/fastq"
        trimmed_dir = id + "/trimmed"
        mapped_dir = id + "/mapped"
        bam_file = "%s/Aligned.sortedByCoord.out.bam" % mapped_dir

        all_done, any_done, done = check_if_results_exist(id)
        log("Results for sample %s already exist" % id if all_done else "Processing sample %s" % id)

        # getting and preprocessing data
        if (not all_done and not verify_bam_file(bam_file)) or simulate_all:

            os.system("mkdir -p -v %s" % fastq_dir)
            os.system("mkdir -p -v %s" % trimmed_dir)
            os.system("mkdir -p -v %s" % mapped_dir)

            # get samples info from ebi.ac.uk
            fields = ["run_accession", "fastq_ftp", "fastq_md5", "fastq_bytes", "library_layout"]
            url = "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=%s&result=read_run&fields=%s" \
                  % (id, ",".join(fields))
            try:
                res = urlopen(url).read().decode("utf-8").split("\n")
            except URLError as e:
                context = ssl._create_unverified_context()
                res = urlopen(url, context=context).read().decode("utf-8").split("\n")
            res = res[1].split("\t")

            layout = res[4]
            urls = res[1].split(";")
            # md5 = res[2].split(";")
            # filenames = res[1].split(";")
            # filenames = [x.split("/")[-1] for x in filenames]

            # download all files
            for url in urls:
                run_task("wget -cNP %s %s" % (fastq_dir, url))
                # run_task("md5sum " + f)

            # remove adapters and low quality reads if didnt before
            if not adapters_trimmed(os.listdir(trimmed_dir)):
                layout = "--paired" if layout.lower() == "paired" else ""
                run_task("trim_galore %s -o %s %s/*" % (layout, trimmed_dir, fastq_dir))

            # mapping
            log("Mapping to genome")
            fastq_files = filter(lambda filename: ".fq" in filename, os.listdir(trimmed_dir))
            fastq_files = [trimmed_dir + "/" + filename for filename in fastq_files]
            fastq_files = " ".join(fastq_files)
            run_task("STAR --runThreadN 24 --genomeDir %s --readFilesIn %s --sjdbOverhang 100 "
                     "--outFileNamePrefix %s/ --twopassMode Basic --outSAMtype BAM SortedByCoordinate "
                     "--outReadsUnmapped Fastx --readFilesCommand zcat" % (genome_path, fastq_files, mapped_dir))

        # running analysis scripts and programs
        bam_ok = verify_bam_file(bam_file)
        if (not all_done and bam_ok) or simulate_all:

            run_task("samtools index %s" % bam_file)

            if not done[COVERAGE]:
                log("Running bedtools coverage")
                # converting BAM to BED and sorting, due to differences in chromosomes order in BAM and BED
                # run_task("bamToBed -i %s > %s/Aligned.bed" % (bam_file, mapped_dir))
                # run_task("bedtools sort -i %s/Aligned.bed > %s/Sorted.bed" % (mapped_dir, mapped_dir))
                # run_task("bedtools coverage -a %s -b %s/Sorted.bed -sorted > %s.coverage"
                #          % (bed_file, mapped_dir, id))
                run_task("bedtools coverage -a %s -b %s -sorted > %s.coverage"
                         % (bed_file, bam_file, id))

            if not done[FLAGSTAT]:
                log("Running samtools flagstat")
                run_task("samtools flagstat %s > %s.flagstat" % (bam_file, id))

            if not done[FASTQC]:
                log("Running fastqc")
                run_task("fastqc -o ./ %s" % bam_file)
                run_task("mv Aligned.sortedByCoord.out_fastqc.zip %s.fastqc.zip" % id)
                run_task("rm Aligned.sortedByCoord.out_fastqc.html")

            if not done[READ_GC]:
                log("Running read_GC.py from RSeQC")
                run_task("read_GC.py -i %s -o %s.read_GC" % (bam_file, id))

            if not done[TIN]:
                log("Running tin.py from RSeQC")
                run_task("tin.py -i %s -r %s" % (bam_file, bed_file))
                run_task("mv Aligned.sortedByCoord.out.tin.xls %s.tin.xls" % id)
                run_task("mv Aligned.sortedByCoord.out.summary.txt %s.tin.summary" % id)

            if not done[INNER_DISTANCE]:
                log("Running inner_distance.py from RSeQC")
                run_task("inner_distance.py -i %s -r %s -o %s" % (bam_file, bed_file, id))

            if not done[GENEBODY_COVERAGE]:
                log("Running geneBody_coverage.py from RSeQC")
                run_task("geneBody_coverage.py -i %s -r %s -o %s" % (bam_file, bed_file, id))

        if not bam_ok and not all_done:
            log("Bam file for sample %s broken. Continuing with next sample" % id)
            broken_bams.append(id)
            continue

        all_done, any_done, done = check_if_results_exist(id)
        # log("Removing temporary files")
        if any_done and bam_ok:
            if os.path.exists(fastq_dir):
                os.system("rm -r %s" % fastq_dir)
            if os.path.exists(trimmed_dir):
                os.system("rm -r %s" % trimmed_dir)
        if all_done:
            if os.path.exists(id):
                os.system("rm -r %s" % id)

    fs = os.listdir('./')
    pdf = list(filter(lambda x: '.pdf' in x, fs))
    r = list(filter(lambda x: [-1] == 'r' and x[-2] == '.', fs))

    if pdf:
        log("Removing .pdf files")
        os.system("rm *.pdf")
    if r:
        log("Removing .r files")
        os.system("rm *.r")

    if broken_bams:
        log("List of samples that failed on BAM file verification:")
        # print(broken_bams)
        log(broken_bams)
    os.chdir('..')
    log("All samples downloaded and processed")


def main():
    usage = "%prog [options] [ EXP_IDs ]" + "\n" + __doc__ + "\n"
    parser = OptionParser(usage)
    # parser.add_option("-i", "--experiment-id", action="store", type="string", dest="exp_id",
    #                   help="ID of experiment in ArrayExpress database, eg. E-MTAB-5616. [required]")
    parser.add_option("-g", "--star-genome", action="store", type="string", dest="star_genome",
                      help="Indexed genome [required]")
    parser.add_option("-r", "--refgene", action="store", type="string", dest="ref_gene_model",
                      help="Reference gene model in BED format required by RSeQC scripts. "
                           "Must be strandard 12-column BED file. [required]")
    parser.add_option("-s", "--simulate", action="store_true", dest="simulate", default=False,
                      help="Prints commands instead of running them, checking if results or files already exist")
    parser.add_option("-S", "--simulate-hard", action="store_true", dest="simulate_all", default=False,
                      help="Prints all commands instead of running them, assuming that no results/files already exist")

    (options, args) = parser.parse_args()

    if not (args and options.ref_gene_model and options.star_genome):
        parser.print_help()
        sys.exit(0)

    if not os.path.exists(options.ref_gene_model):
        print("\n\n" + options.ref_gene_model + " does NOT exists" + "\n", file=sys.stderr)
        parser.print_help()
        sys.exit(0)

    if not os.path.exists(options.star_genome):
        print("\n\n" + options.star_genome + " does NOT exists" + "\n", file=sys.stderr)
        parser.print_help()
        sys.exit(0)

    simulate_all = options.simulate_all
    simulate = True if simulate_all else options.simulate
    genome_path = os.path.abspath(options.star_genome)
    ref_gene_path = os.path.abspath(options.ref_gene_model)

    for id in args:
        get_and_process(id, genome_path, ref_gene_path, simulate, simulate_all)


if __name__ == "__main__":
    main()
