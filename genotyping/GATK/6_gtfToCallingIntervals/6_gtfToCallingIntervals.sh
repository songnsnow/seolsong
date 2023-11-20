
#gtfToCallingIntervals--------------------------------------------------------------------
gtf <- read.table("/reference/star.gencode.v19.transcripts.patched_contigs.gtf", sep="\t")
gtf <- subset(gtf, V3 == "exon")
write.table(data.frame(chrom=gtf[,'V1'], start=gtf[,'V4'], end=gtf[,'V5']), "/reference/exome.bed", quote = F, sep="\t", col.names = F, row.names = F)

awk '{print $1 "\t" ($2 - 1) "\t" $3}' /reference/exome.bed > /reference/exome.fixed.bed

gatk \
    BedToIntervalList \
    -I /reference/exome.fixed.bed \
    -O /reference/interval_list \
    -SD /reference/Homo_sapiens_assembly38.dict
#didn't work.. just used wgs hg38 interval list. also, star.gencode.v19 was hg37.

#try again with gencode.v32.primary_assembly.annotation.gtf
gtf <- read.table("/reference/gencode.v32.primary_assembly.annotation.gtf", sep="\t")
gtf <- subset(gtf, V3 == "exon")
write.table(data.frame(chrom=gtf[,'V1'], start=gtf[,'V4'], end=gtf[,'V5']), "/reference/fexome.bed", quote = F, sep="\t", col.names = F, row.names = F)

awk '{print $1 "\t" ($2 - 1) "\t" $3}' /reference/fexome.bed > /reference/fexome.fixed.bed

gatk \
    BedToIntervalList \
    -I /reference/fexome.fixed.bed \
    -O /reference/38.interval_list \
    -SD /reference/Homo_sapiens_assembly38.dict

#remove 