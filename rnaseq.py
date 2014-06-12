import os
import pysam
import numpy

from contextlib import closing

from dish.pipeline import Pipeline


def has_two_fastqs(job):
    """Are there two .fastq files in the workdir?"""
    fastqs = [fname for fname in os.listdir(job["workdir"])
              if fname.endswith(".fastq")]
    if len(fastqs) > 2:
        raise RuntimeError("something is dramatically wrong")
    return len(fastqs) == 2


def innerdist_from_bowtie(job, logger):
    """Code from:
    https://github.com/chapmanb/bcbio-nextgen/blob/18dea9adb0850cb1b0d706f5d611a0b156ef775e/bcbio/ngsalign/tophat.py#L290-L304"""
    sam = os.path.join(job["workdir"], "bowtie_out.sam")
    with closing(pysam.Samfile(sam)) as work_sam:
        for read in work_sam:
            if read.is_proper_pair and read.is_read1:
                dists.append(abs(read.isize) - 2 * read.rlen)
        median = float(numpy.median(dists))
        deviations = []
        for d in dists:
            deviations.append(abs(d - median))
        # this is the median absolute deviation estimator of the
        # standard deviation
        mad = 1.4826 * float(numpy.median(deviations))
        job["inner_dist"] = int(median)
        job["std_dev"] = int(mad)


# first let's make some jobs
jobs = []
datadir = "/glusterfs/netapp/home2/data/rnaseq_data/READ/"
for dirname in os.listdir(datadir):
    contents = os.listdir(os.path.join(datadir, dirname))
    has_tarball = False
    for fname in contents:
        if fname.endswith(".tar.gz"):
            tarball = fname
            has_tarball = True
    if not has_tarball:
        continue  # there's no data. TODO would be nice to be able to log this?
    jobs.append({"description": dirname,
                 "tarball": os.path.join(datadir, dirname, tarball)})


p = Pipeline("/glusterfs/netapp/home2/PORTERJAMESJ/dishtest/work",
             jobs,
             60,
             "torque",
             "batch", retries=5)

p.start()

# transact on two fastqs being there. we could actually extract the
# manifest from each tarball and transact on those files, but that
# would be really slow because of the gzipping
with p.transaction(has_two_fastqs):
    p.run("tar xvzf {tarball}", max=25)  # running much more than this kills netapp

for job in p.jobs:
    fastqs = [fname for fname in os.listdir(job["workdir"])
              if fname.endswith(".fastq")]
    fastqs.sort()
    job["fastq1"] = fastqs[0]
    job["fastq2"] = fastqs[1]

# trim with trimmomatic
with p.transaction(["{fastq1}.trimmed", "{fastq2}.trimmed"]):
    p.run("java -jar /usr/local/java/Trimmomatic-0.32/trimmomatic-0.32.jar PE "
          " -phred33 -threads 8 {workdir}/{fastq1} {workdir}/{fastq2}"
          " {fastq1}.trimmed /dev/null"
          " {fastq2}.trimmed /dev/null"
          " LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36",
          cores=8,
          mem=8)

# run bowtie
with p.transaction("bowtie_out.sam"):
    p.run("bowtie2"
          " -p 8 -s 100000 -u 250000 -q"
          " -x /glusterfs/data/ICGC1/ref/bcbio-data/tcga/genomes/hg19/bowtie2/hg19"
          " -1 {workdir}/{fastq1}.trimmed -2 {workdir}/{fastq2}.trimmed"
          " -S bowtie_out.sam",
          cores=8,
          mem=8)


# run tophat
with p.transaction("tophat_out/accepted_hits.bam"):
    p.map(innerdist_from_bowtie)
    p.run("tophat"
          " --library-type=fr-unstranded"
          " --num-threads=8"
          " --GTF=/glusterfs/data/ICGC1/ref/bcbio-data/tcga/genomes/hg19/rnaseq/ref-transcripts.gtf"
          " --no-coverage-search"
          " --output-dir=tophat_out"
          " --mate-std-dev={std_dev} --mate-inner-dist={innder_dist}"
          " /glusterfs/data/ICGC1/ref/bcbio-data/tcga/genomes/hg19/bowtie2/hg19"
          " {workdir}/{fastq1}.trimmed"
          " {workdir}/{fastq2}.trimmed")

# run cufflinks
with p.transaction("cufflinks_out/transcripts.gtf"):
    p.run("cufflinks"
          " --num-threads 8"
          " --quiet"
          " --no-update-check"
          " --frag-bias-correct /glusterfs/data/ICGC1/ref/bcbio-data/tcga/genomes/hg19/seq/hg19.fa"
          " --multi-read-correct"
          " --upper-quartile-norm"
          " --GTF /usr/local/share/bcbio-nextgen/genomes/hg19/rnaseq/ref-transcripts.gtf"
          " --output-dir cufflinks_out"
          " {workdir}/tophat_out/accepted_hits.bam")
