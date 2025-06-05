#!/bin/bash

## Usage: sh 0.0_ps_Pipeline_v3.sh SRR3304509 c.0_RNAseq_QuantTERRA.cfg


expID=$1
cfgFile=$2
isDump=$3        
fqPath=$4        


if [ -z "$expID" ] || [ -z "$cfgFile" ]; then
  echo "Usage: sh 0.0_ps_Pipeline_v3.sh SRR_ID config.cfg [is-dump] [fqPath]"
  exit 1
fi


source $cfgFile


fastqDIR="$workdir/fastq"
if [ "$isDump" == "is-dump" ]; then
  if [ -z "$fqPath" ]; then
    echo "Error: is-dump flag was set but no fastq path provided."
    exit 1
  fi
  fastqDIR=$fqPath
fi


function mkFolder(){
  [ ! -d "$workdir/$1" ] && mkdir -p "$workdir/$1"
}
mkFolder reports

# ===================== #
# 3. fastq-dump + Trimgalore
# ===================== #
if [ "$isDump" != "is-dump" ]; then
  mkFolder fastq
  jid_fqdump=$(sbatch -A MST109178 -p ngs${partition}G -J fq_$expID -o reports/fq_$expID.o.txt -e reports/fq_$expID.e.txt -c $ncore --mem=${memory}g \
  0_ps_SRA_dump.sh $workdir/SRA/${expID}.sra $workdir/fastq)

  mkFolder trimed_fq
  jid_Trim=$(sbatch -A MST109178 -p ngs${partition}G -J trim_$expID -o reports/trim_$expID.o.txt -e reports/trim_$expID.e.txt -c $ncore --mem=${memory}g \
  --dependency=afterok:${jid_fqdump/"Submitted batch job "/} \
  1_ps_trimQC.sh $expID $workdir/trimed_fq $fastqDIR $ncore)
else
  mkFolder trimed_fq
  jid_Trim=$(sbatch -A MST109178 -p ngs${partition}G -J trim_$expID -o reports/trim_$expID.o.txt -e reports/trim_$expID.e.txt -c $ncore --mem=${memory}g \
  1_ps_trimQC.sh $expID $workdir/trimed_fq $fastqDIR $ncore)
fi

# ===================== #
# 4. STAR alignment
# ===================== #
mkFolder STAR_align
jid_STAR=$(sbatch -A MST109178 -p ngs${partition}G -J align_$expID -o reports/align_$expID.o.txt -e reports/align_$expID.e.txt -c $ncore --mem=${memory}g \
--dependency=afterok:${jid_Trim/"Submitted batch job "/} \
2_ps_STAR_align.sh $workdir/trimed_fq/$expID $ncore $workdir/STAR_align/$expID)

# ===================== #
# 5. bamCoverage
# ===================== #
mkFolder bamcoverage
jid_Bamcov=$(sbatch -A MST109178 -p ngs${partition}G -J bamcov_$expID -o reports/bamcov_$expID.o.txt -e reports/bamcov_$expID.e.txt -c $ncore --mem=${memory}g \
--dependency=afterok:${jid_STAR/"Submitted batch job "/} \
3_ps_bamcoverage.sh $workdir/STAR_align/${expID}*.bam $ncore $workdir/bamcoverage $expID)

# ===================== #
# 6. HTSeq-count
# ===================== #
mkFolder HTseq_count
jid_HTseq=$(sbatch -A MST109178 -p ngs${partition}G -J HTseq_$expID -o reports/HTseq_$expID.o.txt -e reports/HTseq_$expID.e.txt -c $ncore --mem=${memory}g \
--dependency=afterok:${jid_STAR/"Submitted batch job "/} \
4_ps_HTseq_count.sh $workdir/STAR_align/${expID}*.bam $qGTF $ncore $workdir/HTseq_count/$expID.count.txt)

# ===================== #
# 8. BBduk telomeric content
# ===================== #
mkFolder bbduk_telo
jid_BBduk=$(sbatch -A MST109178 -p ngs${partition}G -J bbduk_$expID -o reports/bbduk_$expID.o.txt -e reports/bbduk_$expID.e.txt -c $ncore --mem=${memory}g \
--dependency=afterok:${jid_Trim/"Submitted batch job "/} \
7_ps_bbduk_telo_contents.sh $expID $workdir/bbduk_telo $workdir/trimed_fq/$expID $teloSeq_ref $ncore)

# ===================== #
# Log submission time
# ===================== #
CYAN='\033[0;36m'
NC='\033[0m'
echo -e "${CYAN}Pipeline submitted at${NC} $(date '+%Y-%m-%d %H:%M:%S')"
