# script to filter out murine reads to remove contamination from pdx samples

module load samtools
module load bwa

cd ../../../PREPROCESS/DNA/
export DATDIR=../../DATA/DNA
export FILTDIR=pdx_filter
mkdir -p $FILTDIR
mkdir -p $FILTDIR/fastq
mkdir -p $FILTDIR/bwa
mkdir -p $FILTDIR/disambiguate
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
           --mem 40G \
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
           --mem 40G \
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
           --mem 40G \
           -t 0-20:00 \
           --wrap="java -jar /homes/keegan/picard-2.23.2/picard.jar SamToFastq \
      I=$DATDIR/$CELLLINE_BAM \
      FASTQ=$FILTDIR/fastq/$CELLLINE_FQ \
      SECOND_END_FASTQ=$FILTDIR/fastq/$CELLLINE_FQ2"
  fi
done

# 2. index hg19 and mm10

  # index file from https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

bwa index ./annotation/hg19.fa.gz

   # index file from https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz

bwa index ./annotation/mm10.fa.gz

########

# xengsort index

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

echo -e '#!/bin/bash

xengsort index $REFPATH/xengsort_index-k25.h5 \
  -H <(zcat $human_tl) <(zcat $human_cdna) \
  -G <(zcat $mouse_tl) <(zcat $mouse_cdna) \
  -k 25 -P 4800000000 3FCVbb:u \
  --hashfunctions linear945:linear9123641:linear349341847 \
  -p 4 --fill 0.94 -T 8 \
  &> $REFPATH/index.log
' > xengsort-index.sh

########


# 3. align each bam file to both hg19 and mm10

for id in 5350 5351 5367 5368 5369 5473 5474; do
  echo "Mapping for ${id}"

  export NORMAL=DFCI-${id}-N-01
  export TUMOR=DFCI-${id}-T-01
  export CELLLINE=DFCI-${id}-CL-01
  export CELLLINE2=DFCI-${id}-CL-01

  if [ $id = "5350" ] || [ $id = "5351" ] || [ $id = "5473" ] || [ $id = "5474" ]; then
    export CELLLINE2=DFCI-${id}-C-01
  fi

  if [ ! -f $FILTDIR/bwa/$TUMOR\_human.bam ]; then
    sbatch -o ../slurm/${id}-tumor-human.out \
           -e ../slurm/${id}-tumor-human.err \
           -J map_tumor_human \
           -n 1 \
           -N 1 \
           --mem 40G \
           -t 0-40:00 \
           --wrap="bwa mem \
               ../../PREPROCESS/DNA/annotation/hg19.fa.gz \
               $FILTDIR/fastq/$TUMOR\_1.fastq \
               $FILTDIR/fastq/$TUMOR\_2.fastq > \
               $FILTDIR/bwa/$TUMOR\_human.sam && \
             samtools sort -o $FILTDIR/bwa/$TUMOR\_human.bam \
               $FILTDIR/bwa/$TUMOR\_human.sam"
    sleep 1 
  fi
  
  if [ ! -f $FILTDIR/bwa/$TUMOR\_mouse.bam ]; then
    sbatch -o ../slurm/${id}-tumor-mouse.out \
           -e ../slurm/${id}-tumor-mouse.err \
           -J map_tumor_mouse \
           -n 1 \
           -N 1 \
           --mem 40G \
           -t 0-40:00 \
           --wrap="bwa mem \
               ../../PREPROCESS/DNA/annotation/mm10.fa.gz \
               $FILTDIR/fastq/$TUMOR\_1.fastq \
               $FILTDIR/fastq/$TUMOR\_2.fastq > \
               $FILTDIR/bwa/$TUMOR\_mouse.sam && \
             samtools sort -o $FILTDIR/bwa/$TUMOR\_mouse.bam \
               $FILTDIR/bwa/$TUMOR\_mouse.sam"
    sleep 1 
  fi

  if [ ! -f $FILTDIR/bwa/$NORMAL\_human.bam ]; then
  	sbatch -o ../slurm/${id}-normal-human.out \
           -e ../slurm/${id}-normal-human.err \
           -J map_normal_human \
           -n 1 \
           -N 1 \
           --mem 40G \
           -t 0-40:00 \
           --wrap="bwa mem \
               ../../PREPROCESS/DNA/annotation/hg19.fa.gz \
               $FILTDIR/fastq/$NORMAL\_1.fastq \
               $FILTDIR/fastq/$NORMAL\_2.fastq > \
               $FILTDIR/bwa/$NORMAL\_human.sam && \
             samtools sort -o $FILTDIR/bwa/$NORMAL\_human.bam \
               $FILTDIR/bwa/$NORMAL\_human.sam"
  fi

  if [ ! -f $FILTDIR/bwa/$NORMAL\_mouse.bam ]; then
    sbatch -o ../slurm/${id}-normal-mouse.out \
           -e ../slurm/${id}-normal-mouse.err \
           -J map_normal_human \
           -n 1 \
           -N 1 \
           --mem 40G \
           -t 0-40:00 \
           --wrap="bwa mem \
               ./../PREPROCESS/DNA/annotation/mm10.fa.gz \
               $FILTDIR/fastq/$NORMAL\_1.fastq \
               $FILTDIR/fastq/$NORMAL\_2.fastq > \
               $FILTDIR/bwa/$NORMAL\_mouse.sam && \
             samtools sort -o $FILTDIR/bwa/$NORMAL\_mouse.bam \
               $FILTDIR/bwa/$NORMAL\_mouse.sam"
  fi

  if [ ! -f $FILTDIR/bwa/$CELLLINE\_human.bam ]; then
  	sbatch -o ../slurm/${id}-cl-human.out \
           -e ../slurm/${id}-cl-human.err \
           -J map_cl_human \
           -n 1 \
           -N 1 \
           --mem 40G \
           -t 0-40:00 \
           --wrap="bwa mem \
               ../../PREPROCESS/DNA/annotation/hg19.fa.gz \
               $FILTDIR/fastq/$CELLLINE2\_1.fastq \
               $FILTDIR/fastq/$CELLLINE2\_2.fastq > \
               $FILTDIR/bwa/$CELLLINE\_human.sam && \
             samtools sort -o $FILTDIR/bwa/$CELLLINE\_human.bam \
               $FILTDIR/bwa/$CELLLINE\_human.sam"
  fi

  if [ ! -f $FILTDIR/bwa/$CELLLINE\_mouse.bam ]; then
    sbatch -o ../slurm/${id}-cl-mouse.out \
           -e ../slurm/${id}-cl-mouse.err \
           -J map_cl_human \
           -n 1 \
           -N 1 \
           --mem 40G \
           -t 0-40:00 \
           --wrap="bwa mem \
               ../../PREPROCESS/DNA/annotation/mm10.fa.gz \
               $FILTDIR/fastq/$CELLLINE2\_1.fastq \
               $FILTDIR/fastq/$CELLLINE2\_2.fastq > \
               $FILTDIR/bwa/$CELLLINE\_mouse.sam && \
             samtools sort -o $FILTDIR/bwa/$CELLLINE\_mouse.bam \
               $FILTDIR/bwa/$CELLLINE\_mouse.sam"
  fi

done



# 5. use disambiguate to remove reads aligning to mm10

conda activate local

for id in 5350 5351 5367 5368 5369 5473 5474; do
  echo "Getting fastq for Tumor and Normal for ${id}"
  export NORMAL=DFCI-${id}-N-01
  export TUMOR=DFCI-${id}-T-01
  export CELLLINE=DFCI-${id}-CL-01

  for TYPE in $NORMAL $TUMOR $CELLLINE; do
    if [ ! -f $FILTDIR/disambiguate/$TYPE\_summary.txt ]; then
      sbatch -o ../slurm/disambiguate-${TYPE}.out \
           -e ../slurm/disambiguate-${TYPE}.err \
           -J disamb-${TYPE} \
           -n 1 \
           -N 1 \
           --mem 40G \
           -t 0-8:00 \
           --wrap="ngs_disambiguate -a bwa \
                   --prefix $TYPE \
                   --output-dir $FILTDIR/disambiguate \
                   $FILTDIR/bwa/$TYPE\_human.bam \
                   $FILTDIR/bwa/$TYPE\_mouse.bam"
    fi 
  done

done


for id in 5350 5351 5367 5368 5369 5473 5474; do
  echo "Getting fastq for Tumor and Normal for ${id}"
  export NORMAL=DFCI-${id}-N-01
  export TUMOR=DFCI-${id}-T-01
  export CELLLINE=DFCI-${id}-CL-01

  for TYPE in $NORMAL $TUMOR $CELLLINE; do
    if [ ! -f $FILTDIR/xenofilter/$TYPE\_summary.txt ]; then
      sbatch -o ../slurm/xenofilter-${TYPE}.out \
           -e ../slurm/xenofilter-${TYPE}.err \
           -J disamb-${TYPE} \
           -n 1 \
           -N 1 \
           --mem 40G \
           -t 0-8:00 \
           --wrap="perl /homes/keegan/XenofilteR/original/xenofilt_PE.pl \
             $FILTDIR/bwa/$TYPE\_human.bam \
             $FILTDIR/bwa/$TYPE\_mouse.bam"
    fi 
  done

done


###########



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


#########





# 6. filter original bam file by read name (keep original alignment) 
#    to only keep reads that didn't align to mm10 in step 5


for id in 5350 5351 5367 5368 5369 5473 5474; do
  echo "Getting fastq for Tumor and Normal for ${id}"
  export NORMAL=DFCI-${id}-N-01
  export TUMOR=DFCI-${id}-T-01
  export CELLLINE=DFCI-${id}-CL-01

  if [ $id = "5350" ] || [ $id = "5351" ]; then
    export CELLLINE=DFCI-${id}-C-01
  fi

  for TYPE in $NORMAL $TUMOR $CELLLINE; do
    if [ ! -f $FILTDIR/cleanbam/$TYPE.bam ]; then
      sbatch -o ../slurm/filter-${TYPE}.out \
           -e ../slurm/filter-${TYPE}.err \
           -J filter-${TYPE} \
           -n 1 \
           -N 1 \
           --mem 40G \
           -t 0-6:00 \
           --wrap="samtools view \
               $FILTDIR/disambiguate/$TYPE\.disambiguatedSpeciesA.bam | \
               cut -f 1 > $FILTDIR/disambiguate/$TYPE\.readnames_keep.txt && \
             java -jar /homes/keegan/picard-2.23.2/picard.jar \
               FilterSamReads \
               I=$DATDIR/$TYPE.bam\
               O=$FILTDIR/cleanbam/$TYPE.bam \
               READ_LIST_FILE=$FILTDIR/disambiguate/$TYPE\.readnames_keep.txt \
               FILTER=includeReadList"
    fi 
  done
done

# 7. clean up

# delete tmp fastq files & remapped bams
rm $FILTDIR/fastq/*.fastq
rm $FILTDIR/bwa/*.bam
rm $FILTDIR/bwa/*.sam

# 8. continue with Mutect2 pipeline using filtered bams
