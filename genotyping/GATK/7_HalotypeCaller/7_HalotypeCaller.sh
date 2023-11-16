
#HalotypeCaller--------------------------------------------------------------------
#download intervals
wget https://storage.googleapis.com/gatk-test-data/intervals/star.gencode.v19.transcripts.patched_contigs.gtf
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list

gatk --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
    HaplotypeCaller \
    -R /reference/Homo_sapiens_assembly38.fasta \
    -I /PATH/outs/bqsr.bam \
    -L /reference/wgs_calling_regions.hg38.interval_list \
    -O /PATH/outs/0.g.vcf.gz \
    -dont-use-soft-clipped-bases \
    -stand-call-conf 20 \
    -D /reference/Homo_sapiens_assembly38.dbsnp138.vcf

gatk --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
    HaplotypeCaller \
    -R /reference/Homo_sapiens_assembly38.fasta \
    -I /PATH/outs/bqsr.bam \
    -L /reference/38.interval_list \
    -O /PATH/outs/i.g.vcf.gz \
    -dont-use-soft-clipped-bases \
    -stand-call-conf 20 \
    -D /reference/Homo_sapiens_assembly38.dbsnp138.vcf
#----------------------------------------------------------------------------------
