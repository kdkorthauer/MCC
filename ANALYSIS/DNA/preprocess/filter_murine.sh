# script to filter out murine reads to remove contamination from pdx samples

module load samtools

cd ../../../PREPROCESS/DNA/
export DATDIR=../../DATA/DNA
export FILTDIR=pdx_filter
mkdir -p $FILTDIR
mkdir -p $FILTDIR/fastq
mkdir -p $FILTDIR/cleanbam
mkdir -p $FILTDIR/xengsort

# 1. convert bam files to fastq

for id in 5350 5351 5367 5368 5369 5473 5474; do
  echo "Getting fastq for Tumor and Normal for ${id}"
  NORMAL_BAM=DFCI-${id}-N-01.bam
  TUMOR_BAM=DFCI-${id}-T-01.bam
  CELLLINE_BAM=DFCI-${id}-CL-01.bam
  CELLLINE_FQ=DFCI-${id}-CL-01_1.fastq
  CELLLINE_FQ2=DFCI-${id}-CL-01_2.fastq

  if [ $id = "5350" ] || [ $id = "5351" ] || [ $id = "5473" ] || [ $id = "5474" ]; then
    CELLLINE_BAM=DFCI-${id}-C-01.bam
    CELLLINE_FQ=DFCI-${id}-C-01_1.fastq
    CELLLINE_FQ2=DFCI-${id}-C-01_2.fastq
  fi
  
  if [ ! -f $FILTDIR/fastq/DFCI-${id}-N-01_1.fastq  ]; then
  	sbatch -o ../slurm/${id}-n-getfastq.out \
           -e ../slurm/${id}-n-getfastq.err \
           -J bam2fastq_n \
           -n 1 \
           -N 1 \
           --mem 20G \
           -t 0-20:00 \
           --wrap="java -jar /homes/keegan/picard-2.23.2/picard.jar SamToFastq \
      I=$DATDIR/$NORMAL_BAM \
      FASTQ=$FILTDIR/fastq/DFCI-${id}-N-01_1.fastq \
      SECOND_END_FASTQ=$FILTDIR/fastq/DFCI-${id}-N-01_2.fastq"
  fi

  if [ ! -f $FILTDIR/fastq/DFCI-${id}-T-01_1.fastq  ]; then
    sbatch -o ../slurm/${id}-t-getfastq.out \
           -e ../slurm/${id}-t-getfastq.err \
           -J bam2fastq_t \
           -n 1 \
           -N 1 \
           --mem 20G \
           -t 0-20:00 \
           --wrap="java -jar /homes/keegan/picard-2.23.2/picard.jar SamToFastq \
      I=$DATDIR/$TUMOR_BAM \
      FASTQ=$FILTDIR/fastq/DFCI-${id}-T-01_1.fastq \
      SECOND_END_FASTQ=$FILTDIR/fastq/DFCI-${id}-T-01_2.fastq"
  fi

  if [ ! -f $FILTDIR/fastq/$CELLLINE_FQ ]; then
    sbatch -o ../slurm/${id}-c-getfastq.out \
           -e ../slurm/${id}-c-getfastq.err \
           -J bam2fastq_c \
           -n 1 \
           -N 1 \
           --mem 20G \
           -t 0-20:00 \
           --wrap="java -jar /homes/keegan/picard-2.23.2/picard.jar SamToFastq \
      I=$DATDIR/$CELLLINE_BAM \
      FASTQ=$FILTDIR/fastq/$CELLLINE_FQ \
      SECOND_END_FASTQ=$FILTDIR/fastq/$CELLLINE_FQ2"
  fi
done

# 2. index hg19 and mm10 for xengsort

export TOPLEVEL_HUMAN=Homo_sapiens.GRCh38.dna.toplevel.fa.gz
export TOPLEVEL_HUMAN_DOWNLOAD=ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
export CDNA_HUMAN=Homo_sapiens.GRCh38.cdna.all.fa.gz
export CDNA_HUMAN_DOWNLOAD=ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
export TOPLEVEL_MOUSE=Mus_musculus.GRCm38.dna.toplevel.fa.gz
export TOPLEVEL_MOUSE_DOWNLOAD=ftp://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
export CDNA_MOUSE=Mus_musculus.GRCm38.cdna.all.fa.gz
export CDNA_MOUSE_DOWNLOAD=ftp://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
export REFPATH=./annotation/

if [ ! -f $REFPATH/$TOPLEVEL_HUMAN ]; then
  wget $TOPLEVEL_HUMAN_DOWNLOAD -P $REFPATH
fi

if [ ! -f $REFPATH/$TOPLEVEL_MOUSE ]; then
  wget $TOPLEVEL_MOUSE_DOWNLOAD -P $REFPATH
fi

if [ ! -f $REFPATH/$CDNA_HUMAN ]; then
  wget $CDNA_HUMAN_DOWNLOAD -P $REFPATH
fi 

if [ ! -f $REFPATH/$CDNA_MOUSE ]; then
  wget $CDNA_MOUSE_DOWNLOAD -P $REFPATH
fi 

export human_tl=$REFPATH/$TOPLEVEL_HUMAN
export mouse_tl=$REFPATH/$TOPLEVEL_MOUSE
export human_cdna=$REFPATH/$CDNA_HUMAN
export mouse_cdna=$REFPATH/$CDNA_MOUSE

echo -e '#!/bin/bash

xengsort index $REFPATH/xengsort_index-k25.h5 \
  -H <(zcat $human_tl) <(zcat $human_cdna) \
  -G <(zcat $mouse_tl) <(zcat $mouse_cdna) \
  -k 25 -P 4800000000 3FCVbb:u \
  --hashfunctions linear945:linear9123641:linear349341847 \
  -p 4 --fill 0.94 -T 8 \
  &> $REFPATH/index.log
' > xengsort-index.sh


conda activate xengsort

if [ ! -f $REFPATH/xengsort_index-k25.h5 ]; then
    sbatch -o ../slurm/index-xengsort.out \
           -e ../slurm/index-xengsort.err \
           -J index-xengsort \
           -n 8 \
           -N 1 \
           --mem 60G \
           -t 0-10:00 \
           xengsort-index.sh
    sleep 1 
fi


# 3. Use xengsort to classify each read 


for id in 5350 5351 5367 5368 5369 5473 5474; do
  echo "Xengsort for ${id}"

  export NORMAL=DFCI-${id}-N-01
  export TUMOR=DFCI-${id}-T-01
  export CELLLINE=DFCI-${id}-CL-01

  if [ $id = "5350" ] || [ $id = "5351" ] || [ $id = "5473" ] || [ $id = "5474" ]; then
    export CELLLINE=DFCI-${id}-C-01
  fi
  
  for TYPE in $NORMAL $TUMOR $CELLLINE; do
      if [ ! -f $FILTDIR/xengsort/$TYPE-host.1.fq ]; then
        sbatch -o ../slurm/${TYPE}.out \
               -e ../slurm/${TYPE}.err \
               -J classify-xengsort \
               -n 8 \
               -N 1 \
               --mem 40G \
               -t 0-20:00 \
               --wrap="xengsort classify \
                        --prefix $FILTDIR/xengsort/$TYPE \
                        --index $REFPATH/xengsort_index-k25.h5 \
                        --classification new \
                        -T 8 \
                        --fastq $FILTDIR/fastq/$TYPE\_1.fastq \
                        --pairs $FILTDIR/fastq/$TYPE\_2.fastq"
        sleep 1 
      fi
  done
done


# 4. get throw away read names

echo -e '#!/bin/bash
' > get-throw.sh

echo -e "awk 'NR%4==1 {print substr(\$1,2)}' \\
  \$FILTDIR/xengsort/\$TYPE-ambiguous.1.fq \\
  | sed 's/..$//' \\
  > \$FILTDIR/xengsort/\$TYPE-readnames_throw.txt \\

echo grabbed read names
" >> get-throw.sh

echo -e "awk 'NR%4==1 {print substr(\$1,2)}' \\
  \$FILTDIR/xengsort/\$TYPE-both.1.fq \\
  | sed 's/..$//' \\
  >> \$FILTDIR/xengsort/\$TYPE-readnames_throw.txt \\

echo grabbed read names
" >> get-throw.sh

echo -e "awk 'NR%4==1 {print substr(\$1,2)}' \\
  \$FILTDIR/xengsort/\$TYPE-graft.1.fq \\
  | sed 's/..$//' \\
  >> \$FILTDIR/xengsort/\$TYPE-readnames_throw.txt \\

echo grabbed read names
" >> get-throw.sh

echo -e "awk 'NR%4==1 {print substr(\$1,2)}' \\
  \$FILTDIR/xengsort/\$TYPE-neither.1.fq \\
  | sed 's/..$//' \\
  >> \$FILTDIR/xengsort/\$TYPE-readnames_throw.txt \\

echo grabbed read names
" >> get-throw.sh


for id in 5350 5351 5367 5368 5369 5473 5474; do
  echo "Getting fastq for Tumor and Normal for ${id}"
  export NORMAL=DFCI-${id}-N-01
  export TUMOR=DFCI-${id}-T-01
  export CELLLINE=DFCI-${id}-CL-01

  if [ $id = "5350" ] || [ $id = "5351" ] || [ $id = "5473" ] || [ $id = "5474" ]; then
    export CELLLINE=DFCI-${id}-C-01
  fi

  for TYPE in $NORMAL $TUMOR $CELLLINE; do
    export TYPE
    if [ ! -f $FILTDIR/xengsort/$TYPE-readnames_throw.txt ]; then
      sbatch -o ../slurm/filter-${TYPE}.out \
           -e ../slurm/filter-${TYPE}.err \
           -J filter-${TYPE} \
           -n 1 \
           -N 1 \
           --mem 20G \
           -t 0-6:00 \
           get-throw.sh
    fi 
  done
done

-rw-r--r-- 1 keegan 6.9M Aug 17 23:50 DFCI-5350-C-01-ambiguous.1.fq
-rw-r--r-- 1 keegan  23M Aug 17 23:50 DFCI-5350-C-01-both.1.fq
-rw-r--r-- 1 keegan 1.1M Aug 17 23:50 DFCI-5350-C-01-graft.1.fq
-rw-r--r-- 1 keegan  49M Aug 17 23:50 DFCI-5350-C-01-neither.1.fq

# clean up 
# rm $FILTDIR/xengsort/*fq

# 5. filter original bam file by read name (keep original alignment) 
#    to only keep reads that didn't align to mm10 in step 5

echo -e '#!/bin/bash
' > filter-bams.sh

echo -e "
java -jar /homes/keegan/picard-2.23.2/picard.jar \\
  FilterSamReads \\
  I=/rafalab/jill/MCC/DATA/DNA/\$TYPE.bam \\
  O=\$FILTDIR/cleanbam/\$TYPE.bam \\
  READ_LIST_FILE=\$FILTDIR/xengsort/\$TYPE-readnames_throw.txt \\
  MAX_RECORDS_IN_RAM=100000 \\
  FILTER=excludeReadList

echo filtering complete
" >> filter-bams.sh


for id in 5350 5351 5367 5368 5369 5473 5474; do
  echo "Getting fastq for Tumor and Normal for ${id}"
  export NORMAL=DFCI-${id}-N-01
  export TUMOR=DFCI-${id}-T-01
  export CELLLINE=DFCI-${id}-CL-01

  if [ $id = "5350" ] || [ $id = "5351" ] || [ $id = "5473" ] || [ $id = "5474" ]; then
    export CELLLINE=DFCI-${id}-C-01
  fi

  for TYPE in $NORMAL $TUMOR $CELLLINE; do
    export TYPE
    if [ ! -f $FILTDIR/cleanbam/$TYPE.bam ]; then
      sbatch -o ../slurm/filter-${TYPE}.out \
           -e ../slurm/filter-${TYPE}.err \
           -J filter-${TYPE} \
           -n 1 \
           -N 1 \
           --mem 40G \
           -t 0-10:00 \
           filter-bams.sh
    fi 
  done
done

# 7. clean up

# delete tmp fastq files
rm $FILTDIR/fastq/*.fastq

# 8. continue with Mutect2 pipeline using filtered bams
