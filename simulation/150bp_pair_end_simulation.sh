#!/bin/bash

read -r -d '' PYTHON_SCRIPT << EOM
from Bio import SeqIO
from Bio.Seq import Seq

## Identified junction points of each chromosome end ##

points = {
    "chr1p": 2705,
    "chr1q": 248384173,
    "chr2p": 3615,
    "chr2q": 242694129,
    "chr3p": 2641,
    "chr3q": 201101334,
    "chr4p": 3269,
    "chr4q": 193572622,
    "chr5p": 2295,
    "chr5q": 182043926,
    "chr6p": 2898,
    "chr6q": 172123884,
    "chr7p": 3418,
    "chr7q": 160565207,
    "chr8p": 2512,
    "chr8q": 146256713,
    "chr9p": 3684,
    "chr9q": 150614275,
    "chr10p": 2635,
    "chr10q": 134754996,
    "chr11p": 1989,
    "chr11q": 135125178,
    "chr12p": 3101,
    "chr12q": 133322264,
    "chr13p": 2541,
    "chr13q": 113563248,
    "chr14p": 2073,
    "chr14q": 101159838,
    "chr15p": 3265,
    "chr15q": 99750312,
    "chr16p": 2322,
    "chr16q": 96327784,
    "chr17p": 1691,
    "chr17q": 84273930,
    "chr18p": 2016,
    "chr18q": 80539047,
    "chr19p": 2286,
    "chr19q": 61704428,
    "chr20p": 2729,
    "chr20q": 66207254,
    "chr21p": 3014,
    "chr21q": 45086114,
    "chr22p": 4579,
    "chr22q": 51321957,
    "chrXp": 1827,
    "chrXq": 154256621
}

### CHM13 reference genome fasta file ###

fasta_file = "./CHM13_reference_genome.fasta"

### Read1 and Read2 fastq generation ###

with open("./simulated_reads_R1.fastq", "w") as output_R1, open("./simulated_reads_R2.fastq", "w") as output_R2:
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = record.seq
        chrom = record.id

        for gap in range(-50, 201):  
            if chrom + "p" in points:
                point = points[chrom + "p"]
                for i in range(300):  
                    read1_start = point - 150 - 150 - gap + i
                    read2_start = point - 150 + i
                    if read2_start + 150 <= len(seq):
                        read1 = seq[read1_start:read1_start + 150].reverse_complement()
                        read2 = seq[read2_start:read2_start + 150]
                        quality = "I" * 150
                        output_R1.write(f"@{chrom}_p_read1_{i}_gap{gap}\n{read1}\n+\n{quality}\n")
                        output_R2.write(f"@{chrom}_p_read2_{i}_gap{gap}\n{read2}\n+\n{quality}\n")

            if chrom + "q" in points:
                point = points[chrom + "q"]
                for i in range(300):  
                    read1_start = point - 150 - 150 - gap + i
                    read2_start = point - 150 + i
                    if read2_start + 150 <= len(seq):
                        read1 = seq[read1_start:read1_start + 150]
                        read2 = seq[read2_start:read2_start + 150].reverse_complement()
                        quality = "I" * 150
                        output_R1.write(f"@{chrom}_q_read1_{i}_gap{gap}\n{read1}\n+\n{quality}\n")
                        output_R2.write(f"@{chrom}_q_read2_{i}_gap{gap}\n{read2}\n+\n{quality}\n")

print("Simulated reads generated successfully.")
EOM

python3 -c "$PYTHON_SCRIPT"

### Alignment of simulated reads ###
STAR \
    --runThreadN 14 \
    --readFilesCommand zcat \
    --genomeDir ./CHM13_STAR_index \
    --readFilesIn  ./simulated_reads_R1.fastq ./simulated_reads_R2.fastq\
    --outFileNamePrefix ./simulated \
    --outSAMtype BAM SortedByCoordinate

samtools index ./simulatedAligned.sortedByCoord.out.bam
samtools view -q 30 ./simulatedAligned.sortedByCoord.out.bam > ./simulatedAligned.sortedByCoord.q30.out.bam
samtools index ./simulatedAligned.sortedByCoord.q30.out.bam

read -r -d '' Ratio_caculate_SCRIPT << EOM
import pysam
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

bamfile = pysam.AlignmentFile(./simulatedAligned.sortedByCoord.out.bam", "rb")

chromosomes = ["chr" + str(i) for i in range(1, 23)] + ["chrX"]
arms = ['p', 'q']
index = pd.MultiIndex.from_product([chromosomes, arms], names=["Chromosome", "Arm"])
matrix = pd.DataFrame(0, index=index, columns=index)

correct_alignments = defaultdict(int)
total_alignments = defaultdict(int)

for read in bamfile.fetch():
    ref_chrom = bamfile.get_reference_name(read.reference_id)
    if ref_chrom in chromosomes:
        if 'p_read' in read.query_name:
            arm = 'p'
        elif 'q_read' in read.query_name:
            arm = 'q'
        else:
            continue
        
        chrom_part = read.query_name.split("_")[0]
        if chrom_part.startswith("chr"):
            chrom_part = chrom_part[3:]  

        simulated_chrom_arm = f"chr{chrom_part}_{arm}"
        print(f"Debug: read.query_name = {read.query_name}, chrom_part = {chrom_part}, simulated_chrom_arm = {simulated_chrom_arm}")

        if simulated_chrom_arm.split('_')[0] not in chromosomes:
            print(f"Warning: simulated chromosome {simulated_chrom_arm.split('_')[0]} not found in chromosomes list")
        if (simulated_chrom_arm.split('_')[0], simulated_chrom_arm.split('_')[1]) not in matrix.index:
            print(f"Warning: simulated chrom_arm {simulated_chrom_arm} not found in matrix index")

        matrix.at[(simulated_chrom_arm.split('_')[0], simulated_chrom_arm.split('_')[1]), (ref_chrom, arm)] += 1
        total_alignments[(simulated_chrom_arm.split('_')[0], simulated_chrom_arm.split('_')[1])] += 1
        if simulated_chrom_arm == f"{ref_chrom}_{arm}":
            correct_alignments[(ref_chrom, arm)] += 1

bamfile.close()

for chrom in chromosomes:
    for arm in arms:
        simulated_arm = f"{chrom}_{arm}"
        print(f"{simulated_arm}: Total alignments: {total_alignments[(chrom, arm)]}")

alignment_ratios = {k: correct_alignments[k] / total_alignments[k] if total_alignments[k] > 0 else 0 for k in total_alignments}

EOM

python3 -c "$Ratio_caculate_SCRIPT"

