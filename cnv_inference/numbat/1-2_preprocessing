# run preprocessing script ----------------------------------------------------------------------------------------------------
docker start numbat  # open docker container
docker exec -it numbat bash  # access running container

# install singularity
conda install -c conda-forge singularity

# singularity run \
#    --bind /mnt/isilon/ \
#    --cleanenv -H /mnt/isilon/cscb/xxx/code/04b_NumBat \
#    /mnt/isilon/cscb/xxx/code/04b_NumBat/numbat-rbase_latest.sif

# make allele count by donor
Rscript numbat/inst/bin/pileup_and_phase.R \
    --label rcc.1 \
    --samples rcc.1 \
    --bams /mnt/mydata/input/possorted_genome_bam.bam \
    --barcodes /mnt/mydata/input/rcc.1_barcodes.tsv \
    --outdir /mnt/mydata/output/rcc.1 \
    --gmap /Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
    --snpvcf /data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
    --paneldir /data/1000G_hg38 \
    --ncores 4

Rscript numbat/inst/bin/pileup_and_phase.R \
    --label rcc.2 \
    --samples rcc.2 \
    --bams /mnt/mydata/input/possorted_genome_bam.bam \
    --barcodes /mnt/mydata/input/rcc.2_barcodes.tsv \
    --outdir /mnt/mydata/output/rcc.2 \
    --gmap /Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
    --snpvcf /data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
    --paneldir /data/1000G_hg38 \
    --ncores 8

Rscript numbat/inst/bin/pileup_and_phase.R \
    --label rcc.3 \
    --samples rcc.3 \
    --bams /mnt/mydata/input/possorted_genome_bam.bam \
    --barcodes /mnt/mydata/input/rcc.3_barcodes.tsv \
    --outdir /mnt/mydata/output/rcc.3 \
    --gmap /Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
    --snpvcf /data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
    --paneldir /data/1000G_hg38 \
    --ncores 20

Rscript numbat/inst/bin/pileup_and_phase.R \
    --label rcc.4 \
    --samples rcc.4 \
    --bams /mnt/mydata/input/possorted_genome_bam.bam \
    --barcodes /mnt/mydata/input/rcc.4_barcodes.tsv \
    --outdir /mnt/mydata/output/rcc.4 \
    --gmap /Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
    --snpvcf /data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
    --paneldir /data/1000G_hg38 \
    --ncores 20

Rscript numbat/inst/bin/pileup_and_phase.R \
    --label rcc.5 \
    --samples rcc.5 \
    --bams /mnt/mydata/input/possorted_genome_bam.bam \
    --barcodes /mnt/mydata/input/rcc.5_barcodes.tsv \
    --outdir /mnt/mydata/output/rcc.5 \
    --gmap /Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
    --snpvcf /data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
    --paneldir /data/1000G_hg38 \
    --ncores 20

Rscript numbat/inst/bin/pileup_and_phase.R \
    --label rcc.6 \
    --samples rcc.6 \
    --bams /mnt/mydata/input/possorted_genome_bam.bam \
    --barcodes /mnt/mydata/input/rcc.6_barcodes.tsv \
    --outdir /mnt/mydata/output/rcc.6 \
    --gmap /Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
    --snpvcf /data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
    --paneldir /data/1000G_hg38 \
    --ncores 20

Rscript numbat/inst/bin/pileup_and_phase.R \
    --label rcc.7 \
    --samples rcc.7 \
    --bams /mnt/mydata/input/possorted_genome_bam.bam \
    --barcodes /mnt/mydata/input/rcc.7_barcodes.tsv \
    --outdir /mnt/mydata/output/rcc.7 \
    --gmap /Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
    --snpvcf /data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
    --paneldir /data/1000G_hg38 \
    --ncores 20

Rscript numbat/inst/bin/pileup_and_phase.R \
    --label rcc.8 \
    --samples rcc.8 \
    --bams /mnt/mydata/input/possorted_genome_bam.bam \
    --barcodes /mnt/mydata/input/rcc.8_barcodes.tsv \
    --outdir /mnt/mydata/output/rcc.8 \
    --gmap /Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
    --snpvcf /data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
    --paneldir /data/1000G_hg38 \
    --ncores 20

Rscript numbat/inst/bin/pileup_and_phase.R \
    --label norm_1 \
    --samples norm_1 \
    --bams /mnt/mydata/input/possorted_genome_bam.bam \
    --barcodes /mnt/mydata/input/norm_1_barcodes.tsv \
    --outdir /mnt/mydata/output/norm_1 \
    --gmap /Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
    --snpvcf /data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
    --paneldir /data/1000G_hg38 \
    --ncores 20

# make own references
# count_mat is a gene x cell raw count matrices
# cell_annot is a dataframe with columns "cell" and "group"
ref_internal = aggregate_counts(count_mat, cell_annot)  # didn't use
