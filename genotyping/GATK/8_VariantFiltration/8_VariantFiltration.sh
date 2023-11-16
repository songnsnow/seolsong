
#VariantFiltration-----------------------------------------------------------------
gatk VariantFiltration \
    -R /reference/Homo_sapiens_assembly38.fasta \
    -V /PATH/outs/i.g.vcf.gz \
    -O /PATH/outs/Name.vcf.gz \
    --window 35 \
    --cluster 3 \
    --filter-name "FS" \
    --filter "FS > 30.0" \
    --filter-name "QD" \
    --filter "QD < 2.0"
#----------------------------------------------------------------------------------
