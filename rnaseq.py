import os
import tarfile
from dish.pipeline import Pipeline


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
             30,
             "torque",
             "batch")
p.start()

# define pipeline functions with signature f(job, logger)
def get_tar_manifest(job, logger):
    tar = tarfile.open(job["tarball"])
    job["tar_manifest"] = [m.name for m in tar.getmembers()]

with p.group(max=20):
    # first we get a manifest from each fastq so we know what to
    # transact on. it's somewhat annoying to have to do this outside
    # the transaction because it's a fair amount of time wasted if the
    # untarring is already done. really this is the fault of using
    # .tar.gz archives because you can't quickly get a manifest out of
    # them.
    p.map(get_tar_manifest)
    # now transact on the existance of the fastqs from the tarball
    with p.transaction(["{tar_manifest[0]}","{tar_manifest[1]}"]):
        p.run("tar xvzf {tarball}")
    # at this point we can be sure that either all the tars have been
    # successfuly extracted OR an error has been thrown
    for job in p.jobs:
        # rename things for convenience
        job["fastq1"] = job["tar_manifest"][0]
        job["fastq2"] = job["tar_manifest"][1]

with p.transaction(["{fastq1}.trimmed", "{fastq2}.trimmed"]):
    p.run("java -jar /usr/local/java/Trimmomatic-0.32/trimmomatic-0.32.jar PE "
          " -phred33 -threads 8 {fastq1} {fastq2}"
          " {fastq1}.trimmed /dev/null"
          " {fastq2}.trimmed /dev/null"
          " LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36",
          cores=8,
          mem=8)

with p.transaction("bowtie_out.sam"):
    p.run("/usr/local/bin/bowtie2"
          " -p 8 -s 100000 -u 250000 -q"
          " -x /glusterfs/data/ICGC1/ref/bcbio-data/tcga/genomes/hg19/bowtie2/hg19"
          " -1 {fastq1}.trimmed -2 {fastq2}.trimmed"
          " -S bowtie_out.sam")
