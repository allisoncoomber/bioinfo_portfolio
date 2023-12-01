#!/bin/bash
  
/home/al/apps/gatk-4.2.6.1/gatk SelectVariants \
        -V chapter2_nuclear.vcf.gz \
        -select-type SNP \
        -O chapter2_nuclear.snps.vcf.gz

/home/al/apps/gatk-4.2.6.1/gatk SelectVariants \
        -V chapter2_nuclear.vcf.gz \
        -select-type INDEL \
        -O chapter2_nuclear.indels.vcf.gz

/home/al/apps/gatk-4.2.6.1/gatk VariantFiltration \
        -V chapter2_nuclear.snps.vcf.gz \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -O chapter2_nuclear.snps.filtered.vcf.gz

/home/al/apps/gatk-4.2.6.1/gatk VariantFiltration \
        -V chapter2_nuclear.indels.vcf.gz \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "FS > 200.0" --filter-name "FS200" \
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
        -O chapter2_nuclear.indels.filtered.vcf.gz