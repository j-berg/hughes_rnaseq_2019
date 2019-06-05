#!/bin/bash
#SBATCH --account=hci-kp
#SBATCH --partition=hci-kp
#SBATCH --nodes=1
#SBATCH -C "c20"
#SBATCH --time=177:02:12
#SBATCH -o /uufs/chpc.utah.edu/common/home/hcibcore/atlatl/job/u0672406_30649_218/stdout.txt
#SBATCH -e /uufs/chpc.utah.edu/common/home/hcibcore/atlatl/job/u0672406_30649_218/stderr.txt
#SBATCH --job-name=u0672406_30649_218
MEMTOTAL=`free | grep Mem | awk '{ print $2 }'`
MEMGB=`expr $MEMTOTAL / 1048576`
SMGB=`expr $MEMGB - 2`
NCPU=`nproc`
GCT=$NCPU
DISKAVAIL=`df -hP . | tail -n 1 | awk '{print $4}'`
if [ $NCPU -gt 8 ]
then
GCT=`expr $NCPU \* 5 / 8 + 3`
fi
HOST=`hostname`
echo "Hostname: " $HOST
echo "Total cpu: " $NCPU
echo "GC threads: " $GCT
echo "Total memory: " $MEMGB
echo "Java memory: " $SMGB
echo "Disk available: " $DISKAVAIL
cleanup()
{
cd /uufs/chpc.utah.edu/common/home/hcibcore/atlatl/job/u0672406_30649_218
cp -ur /scratch/local/u0672406_30649_218/* /uufs/chpc.utah.edu/common/home/hcibcore/atlatl/job/u0672406_30649_218
rm -fr /scratch/local/u0672406_30649_218
cat /uufs/chpc.utah.edu/common/home/hcibcore/atlatl/job/u0672406_30649_218/files_to_be_removed |xargs rm -f
}
failclean()
{    
cleanup
touch tmp.fail
echo failure
exit 1
}
trap 'failclean' ERR TERM
find /uufs/chpc.utah.edu/common/home/hcibcore/atlatl/job/u0672406_30649_218 -type f -print | egrep -v 'cmd.txt|pbs.sh|std....txt' > /uufs/chpc.utah.edu/common/home/hcibcore/atlatl/job/u0672406_30649_218/files_to_be_removed
rm -fr /scratch/local/* 2>/dev/null || true
umask 000
mkdir -p /scratch/local/u0672406_30649_218
cd /scratch/local/u0672406_30649_218
cp -r /uufs/chpc.utah.edu/common/home/hcibcore/atlatl/job/u0672406_30649_218/* .
#e chris.conley@hci.utah.edu -ab
#c kingspeak_20
#a A4301
DB=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/data/S_cerevisiae/Sc3
INDEX=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/data/S_cerevisiae/Sc3/star50
GTF=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/data/S_cerevisiae/Sc3/SacCer3_R64_all_genes_NoDubious_UTR.gtf
REFFLAT=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/data/S_cerevisiae/Sc3/SacCer3_R64_all_genes_NoDubious_UTR.refFlat
RIBOINT=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/data/S_cerevisiae/Sc3/Saccharomyces_cerevisiae.R64-1-1.rRNA.interval
CHROM=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/data/S_cerevisiae/Sc3/chrom.sizes
APP=/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/app
FASTQC=$APP/FastQC/0.11.5/fastqc
STAR=$APP/STAR/2.5.2b/STAR
FEATCOUNT=$APP/Subread/1.5.1/bin/featureCounts
SAMTOOLS=$APP/samtools/1.3/samtools
BIGWIG=$APP/UCSC/bedGraphToBigWig
RnaSeqMetrics=$APP/picard/1.87/CollectRnaSeqMetrics.jar
CUTADAPT=/uufs/chpc.utah.edu/common/home/u6008939/miniconda3/bin/cutadapt
GZ=`echo *.txt.gz`
OUT=`echo ${GZ%%_*}`
mv $GZ $OUT.fq.gz
$CUTADAPT -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -m 16 -o $OUT.cutadapt.fq.gz $OUT.fq.gz
mv $OUT.cutadapt.fq.gz $OUT.fq.gz
gunzip $OUT.fq.gz
$FASTQC -f fastq $OUT.fq
$STAR --genomeDir $INDEX \
--readFilesIn $OUT.fq \
--runThreadN $NCPU \
--twopassMode Basic \
--alignIntronMax 2000 \
--limitBAMsortRAM 1000000000 \
--outSAMtype BAM SortedByCoordinate \
--outWigType bedGraph \
--outWigStrand Unstranded \
--runMode alignReads \
--clip3pAdapterSeq AGATCGGAAGAGCACACGTCTGAACTCCAGTCA AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
mv Aligned.sortedByCoord.out.bam $OUT.bam
mv Log.final.out $OUT.Log.final.out
$FEATCOUNT -T 16 -s 2 -a $GTF -o $OUT.counts $OUT.bam
$SAMTOOLS index $OUT.bam
$SAMTOOLS idxstats $OUT.bam | sort -V > $OUT.idxstats
/uufs/chpc.utah.edu/common/home/hcibcore/atlatl/app/java/jre1.7.0_40/bin/java -Xmx20G -jar $RnaSeqMetrics REF_FLAT=$REFFLAT \
STRAND=SECOND_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=$RIBOINT \
I=$OUT.bam O=$OUT.rna_metrics
$BIGWIG Signal.Unique.str1.out.bg $CHROM $OUT.unique.bw
$BIGWIG Signal.UniqueMultiple.str1.out.bg $CHROM $OUT.multiple.bw
rm $OUT.fq
rm Signal.Unique.str1.out.bg
rm Signal.UniqueMultiple.str1.out.bg
rm -rf _STARgenome
rm -rf _STARpass1

cleanup
echo success
