#!/bin/bash
#PBS -N mapcontrol
#PBS -l walltime=02:30:00,select=1:ncpus=1:mem=120gb 
#PBS -A st-kdkortha-1 
#PBS -o /scratch/st-kdkortha-1/_pbs/MCC/out.txt
#PBS -e /scratch/st-kdkortha-1/_pbs/MCC/err.txt
#PBS -V 


TMPDIR=/scratch/st-kdkortha-1/MCC
mkdir -p $TMPDIR/PREPROCESS/RNA/reference/star
mkdir -p $TMPDIR/PREPROCESS/RNA/control_lines


## Script to align raw RNAseq reads from control lines
eval "$(conda shell.bash hook)"
conda activate cfmedip

# download reference sequence 
if [[ ! -e $TMPDIR/PREPROCESS/RNA/reference/gencode.v19.annotation.gtf ]]; then
  wget -O $TMPDIR/PREPROCESS/RNA/reference/gencode.v19.annotation.gtf "https://data.broadinstitute.org/snowman/hg19/star/gencode.v19.annotation.gtf"
fi

# also download: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz

# generate star index files
if [[ ! -e /scratch/st-kdkortha-1/MCC/PREPROCESS/RNA/reference/star//genomeParameters.txt ]]; then
  STAR --runThreadN 12 --runMode genomeGenerate --genomeDir $TMPDIR/PREPROCESS/RNA/reference/star/ \
  --genomeFastaFiles $TMPDIR/PREPROCESS/RNA/reference/GRCh37.p13.genome.fa \
  --outFileNamePrefix $TMPDIR/PREPROCESS/RNA/control_lines/$group \
  --outTmpDir $TMPDIR/PREPROCESS/RNA/reference/tmp
fi

FQFILES=$(find /arc/project/st-kdkortha-1/MCC/DATA/RNA/control_lines/ -type f -name '*.fastq.gz')

for i in $FQFILES; do
    echo $i

    samp=$(basename $i .fastq.gz)
    pref=$(dirname $i)
    group="$(basename "$(dirname "$i")")"
    mkdir -p $TMPDIR/PREPROCESS/RNA/control_lines/$group

    if [[ ! -e $TMPDIR/PREPROCESS/RNA/control_lines/$group/$samp\_Aligned.out.bam ]]; then
      STAR  --genomeDir $TMPDIR/PREPROCESS/RNA/reference/star/ \
      --runThreadN 1 \
      --readFilesIn $i \
      --twopassMode Basic \
      --outSAMstrandField intronMotif \
      --alignIntronMax 1000000 \
      --alignMatesGapMax 1000000 \
      --sjdbScore 2 \
      --outSAMtype BAM Unsorted \
      --outSAMattributes NH HI NM MD AS XS \
      --outFilterType BySJout \
      --outSAMunmapped Within \
      --genomeLoad NoSharedMemory \
      --outFilterScoreMinOverLread 0 \
      --outFilterMatchNminOverLread 0 \
      --outFilterMismatchNmax 999 \
      --outFilterMultimapNmax 20 \
      --outFileNamePrefix $TMPDIR/PREPROCESS/RNA/control_lines/$group/$samp\_ \
      --readFilesCommand zcat \
      --sjdbGTFfile $TMPDIR/PREPROCESS/RNA/reference/gencode.v19.annotation.gtf 
    fi


  if [[ ! -e $TMPDIR/PREPROCESS/RNA/control_lines/$group/$samp\_Aligned.out_sorted.bam ]]; then
    rm $TMPDIR/PREPROCESS/RNA/control_lines/$group/$samp\_Aligned.out_sorted.bam.tmp*
    samtools sort -o $TMPDIR/PREPROCESS/RNA/control_lines/$group/$samp\_Aligned.out_sorted.bam $TMPDIR/PREPROCESS/RNA/control_lines/$group/$samp\_Aligned.out.bam
    
  fi
    
  if [[ ! -e $TMPDIR/PREPROCESS/RNA/control_lines/$group/$samp\_Aligned.out_duplicates_marked.bam ]]; then
    picard MarkDuplicates \
      READ_ONE_BARCODE_TAG=RX \
      INPUT=$TMPDIR/PREPROCESS/RNA/control_lines/$group/$samp\_Aligned.out_sorted.bam \
      OUTPUT=$TMPDIR/PREPROCESS/RNA/control_lines/$group/$samp\_Aligned.out_duplicates_marked.bam \
      METRICS_FILE=$TMPDIR/PREPROCESS/RNA/control_lines/$group/$samp\_duplicate_metrics \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      TMP_DIR=$TMPDIR/PREPROCESS/RNA/control_lines/$group \
      VALIDATION_STRINGENCY=SILENT \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true \
      MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000 \
      MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 \
      SORTING_COLLECTION_SIZE_RATIO=0.25 \
      TAG_DUPLICATE_SET_MEMBERS=false \
      REMOVE_SEQUENCING_DUPLICATES=false \
      TAGGING_POLICY=DontTag \
      CLEAR_DT=true \
      DUPLEX_UMI=false \
      ADD_PG_TAG_TO_READS=true \
      REMOVE_DUPLICATES=false \
      ASSUME_SORTED=false \
      DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES \
      PROGRAM_RECORD_ID=MarkDuplicates \
      PROGRAM_GROUP_NAME=MarkDuplicates \
      MAX_OPTICAL_DUPLICATE_SET_SIZE=300000 \
      VERBOSITY=INFO \
      QUIET=false \
      COMPRESSION_LEVEL=5 \
      MAX_RECORDS_IN_RAM=500000 \
      USE_JDK_DEFLATER=false \
      USE_JDK_INFLATER=false 
  fi
done

