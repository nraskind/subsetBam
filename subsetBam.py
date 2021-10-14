import os
import pysam
from collections import defaultdict
from itertools import chain
from sys import argv


# subset the reference fasta for realignment in a later step
def subset_fasta(fasta_in, region_one, region_two):
    if not os.path.exists("./out/"):
        os.mkdir("./out/")
    with pysam.FastaFile(fasta_in) as og_fasta:
        with open("./out/out.fa", "w") as new_fasta:
            new_fasta.write(">NewContig\n{}{}".format(og_fasta.fetch(region=region_one),
                                                      og_fasta.fetch(region=region_two)))
        os.system("bwa index ./out/out.fa")


def write_orphans(bamfile_path, read_dict):
    bam_name, _ = os.path.splitext(bamfile_path)
    with open(bam_name + "_orphans.fastq", mode="w") as orphans:
        for qname in read_dict:
            read = read_dict[qname][0] if read_dict[qname][0] else read_dict[qname][1]
            orphans.write("@{}\n{}\n+\n{}\n".format(read.query_name,
                                                    read.get_forward_sequence(),
                                                    pysam.qualities_to_qualitystring(read.get_forward_qualities())))


# cycle through each region and return read pairs
# writes the orphaned reads after completion
def read_pair_generator(bamfile_path, region_one, region_two):
    with pysam.AlignmentFile(bamfile_path, "rb") as bamfile:
        read_dict = defaultdict(lambda: [None, None])

        # remember this is 1-indexed whereas biobambam2 does 0 indexing
        for read in chain(bamfile.fetch(region=region_one, multiple_iterators=True),
                          bamfile.fetch(region=region_two, multiple_iterators=True)):

            if read.is_secondary:
               continue

            qname = read.query_name
            if qname not in read_dict:
                if read.is_read1:
                    read_dict[qname][0] = read
                else:
                    read_dict[qname][1] = read
            else:
                if read.is_read1:
                    yield read, read_dict[qname][1]
                else:
                    yield read_dict[qname][0], read
                del read_dict[qname]

        # write orphaned reads to a third fastq file
        write_orphans(bamfile_path, read_dict)


# writes paired end reads from a bam file to two separate fastq files in preparation for realignment
# additionally, writes orphaned reads to a third fastq file in preperation for realignment
def bamtofastq(bamfile_path, region_one, region_two):
    bam_name, _ = os.path.splitext(bamfile_path)
    with open(bam_name + "_1.fastq", mode="w") as read1s:
        with open(bam_name + "_2.fastq", mode="w") as read2s:
            for read1, read2 in read_pair_generator(bamfile_path, region_one, region_two):
                read1s.write("@{}/1\n{}\n+\n{}\n".format(read1.query_name,
                                                         read1.get_forward_sequence(),
                                                         pysam.qualities_to_qualitystring(read1.get_forward_qualities())))
                read2s.write("@{}/2\n{}\n+\n{}\n".format(read2.query_name,
                                                         read2.get_forward_sequence(),
                                                         pysam.qualities_to_qualitystring(read2.get_forward_qualities())))


# realigns paired end sequencing data, then orphaned reads
# finally, merges them together into a single file
def realign(bamfile_path):
    # need samtools and bwa, realistically samtools can be replaced with pysam
    bam_name, _ = os.path.splitext(bamfile_path)
    os.system(f"bwa mem ./out/out.fa {bam_name}_1.fastq {bam_name}_2.fastq | samtools view -b | samtools sort -o paired.bam")
    os.system(f"bwa mem ./out/out.fa {bam_name}_orphans.fastq | samtools view -b | samtools sort -o orphaned.bam")
    os.system("samtools merge -f ./out/out.bam paired.bam orphaned.bam")
    os.system("samtools index ./out/out.bam")


def clean_up(bamfile_path):
    bam_name, _ = os.path.splitext(bamfile_path)
    os.system(f"rm {bam_name}_1.fastq {bam_name}_2.fastq {bam_name}_orphans.fastq paired.bam orphaned.bam")
    os.system("rm ./out/out.fa.amb ./out/out.fa.ann ./out/out.fa.bwt ./out/out.fa.pac ./out/out.fa.sa")


def main(fasta_in, bamfile_path, region_one, region_two):
    subset_fasta(fasta_in, region_one, region_two)
    bamtofastq(bamfile_path, region_one, region_two)
    realign(bamfile_path)
    clean_up(bamfile_path)


if __name__ == '__main__':
    # usage: python3 subsetbam.py /path/to/reference/fasta /path/to/bam/file "region1 (1-indexed)" "region2 (1-indexed)"
    # output: out.fa* out.bam out.bai
    fasta_in, bamfile_path, region_one, region_two = argv[1], argv[2], argv[3], argv[4]
    main(fasta_in, bamfile_path, region_one, region_two)
