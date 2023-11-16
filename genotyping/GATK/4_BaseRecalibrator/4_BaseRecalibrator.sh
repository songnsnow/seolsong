#BaseRecalibrator------------------------------------------------------------------
#download known indels
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

#Add read groups
picard AddOrReplaceReadGroups \
    I=/PATH/outs/Split.bam \
    O=/PATH/outs/Added.bam \
    RGID=5 \
    RGLB=lib1 \
    RGPL=illumina \
    RGPU=unit1 \
    RGSM=CSB

samtools view -H /PATH/outs/Added.bam | grep '^@RG'

gatk IndexFeatureFile \
    -F /reference/Mills_and_1000G_gold_standard.indels.hg38.vcf

gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
            -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
            -Xloggc:gc_log.log -Xms4000m" \
    BaseRecalibrator \
    -R /reference/Homo_sapiens_assembly38.fasta \
    -I /PATH/outs/Added.bam \
    --use-original-qualities \
    -known-sites /reference/Homo_sapiens_assembly38.dbsnp138.vcf \
    -known-sites /reference/Mills_and_1000G_gold_standard.indels.hg38.vcf \
    -O /PATH/outs/Recal.table
#----------------------------------------------------------------------------------
