#!/bin/bash

CON="/home/varma/softwares/annovar/convert2annovar.pl"
VAR="/home/varma/softwares/annovar/annotate_variation.pl"
HUMANDB="/home/varma/softwares/annovar/humandb/"


######---- split files for CADD scores------
WDIR="/home/varma/proj/work_WGS/"
TMP="/home/varma/proj/work_WGS/TMP/"
INPUT="/media/varma/My_Passport_1/Backup_Varma/RenduOsler2/V6_WGS_VCF/"
OUTPUT="/media/varma/My_Passport_1/Backup_Varma/RenduOsler2/V6_WGS_VCF/"

#####----INPUT------

#./runwfall_annotation_norm_ANNOVAR_GATK_forallsamples_HHT2.bash V6_12345_1-22



##############----examples------
#(4) Run either score.sh or score_anno.sh on a gzip or block-gezipped VCF file:
#bin/score.sh path-to-variants/myVariants.vcf.gz path-to-output/output.tsv.gz
#bin/score.sh path-to-variants/myVariants.vcf.gz path-to-output/output.tsv.gz
#EX:
#bin/score.sh test.vcf.gz test_output.tsv.gz
#bin/score_anno.sh  ex3.vcf.gz  ex3.tsv.gz
#/home/varma/softwares/CADD_v1.2/bin/score_anno.sh  /home/varma/softwares/CADD_v1.2/ex3.vcf.gz  /home/varma/softwares/CADD_v1.2/ex3.tsv.gz

#######---dowload DB------

#/home/varma/softwares/annovar/annotate_variation.pl -buildver hg19 -webfrom annovar -downdb refGene /home/varma/softwares/annovar/humandb/
#/home/varma/softwares/annovar/annotate_variation.pl -buildver hg19 -webfrom annovar -downdb knownGene /home/varma/softwares/annovar/humandb/

#/home/varma/softwares/annovar/annotate_variation.pl -buildver hg19 -downdb cytoBand /home/varma/softwares/annovar/humandb/
#/home/varma/softwares/annovar/annotate_variation.pl -buildver hg19 -downdb genomicSuperDups /home/varma/softwares/annovar/humandb/

#/home/varma/softwares/annovar/annotate_variation.pl -buildver hg19 -webfrom annovar -downdb 1000g2015aug /home/varma/softwares/annovar/humandb/
#/home/varma/softwares/annovar/annotate_variation.pl -buildver hg19 -webfrom annovar -downdb avsnp138 /home/varma/softwares/annovar/humandb/
#/home/varma/softwares/annovar/annotate_variation.pl -buildver hg19 -webfrom annovar -downdb avsnp147 /home/varma/softwares/annovar/humandb/

#/home/varma/softwares/annovar/annotate_variation.pl -buildver hg19 -webfrom annovar -downdb clinvar_20170130 /home/varma/softwares/annovar/humandb/
##/home/varma/softwares/annovar/annotate_variation.pl -buildver hg19 -webfrom annovar -downdb exac03 /home/varma/softwares/annovar/humandb/
##/home/varma/softwares/annovar/annotate_variation.pl -buildver hg19 -webfrom annovar -downdb esp6500siv2_all /home/varma/softwares/annovar/humandb/
#/home/varma/softwares/annovar/annotate_variation.pl -buildver hg19 -webfrom annovar -downdb dbnsfp33a /home/varma/softwares/annovar/humandb/
#/home/varma/softwares/annovar/annotate_variation.pl -buildver hg19 -webfrom annovar -downdb gerp++gt2 /home/varma/softwares/annovar/humandb/
#/home/varma/softwares/annovar/annotate_variation.pl -buildver hg19 -webfrom annovar -downdb gnomad_exome /home/varma/softwares/annovar/humandb/
#/home/varma/softwares/annovar/annotate_variation.pl -buildver hg19 -webfrom annovar -downdb gnomad_genome /home/varma/softwares/annovar/humandb/

#/home/varma/softwares/annovar/annotate_variation.pl -buildver hg19 -webfrom annovar -downdb fathmm /home/varma/softwares/annovar/humandb/
#/home/varma/softwares/annovar/annotate_variation.pl -buildver hg19 -webfrom annovar -downdb gwava /home/varma/softwares/annovar/humandb/
#/home/varma/softwares/annovar/annotate_variation.pl -buildver hg19 -webfrom annovar -downdb eigen /home/varma/softwares/annovar/humandb/
#/home/varma/softwares/annovar/annotate_variation.pl -buildver hg19 -webfrom annovar -downdb intervar_20170202 /home/varma/softwares/annovar/humandb/

###############----run work flow-------
#/home/varma/softwares/CADD_v1.2/bin/score_anno.sh /media/varma/Maxtor/work_WGS/V6_12345_4_variants_norm_annovar.hg19_multianno.vcf.gz /media/varma/Maxtor/work_WGS/V6_12345_4_variants_norm_annovar.hg19_multianno.tsv.gz

#cat -n V6_12345_10_variants_step2.vcf | less
#except header
#tail -n+145 V6_12345_10_variants_step2.vcf | less
#V6_12345_22_variants_SNPnINDELcalls_PED_PhasebyT.vcf
#for ch in {1..22} #10 12 19 22 #{6..22}
#do
######-----preparing left-normalization vcf file-------

#mkdir -p $TMP

#/home/varma/softwares/bcftools-1.1/bcftools norm -m-both -o ''$TMP$1''_variants_GATK_step1.vcf ''$INPUT$1''_variants_GATK_step1_norm.vcf

#/home/varma/softwares/bcftools-1.1/bcftools norm -f /home/varma/softwares/human_g1k_v37.fasta -o ''$OUTPUT$1''_variants_GATK_step1_norm.vcf ''$TMP$1''_variants_GATK_step1.vcf

######-----without CADD scores ---------------
#/home/varma/softwares/annovar/table_annovar.pl ''$OUTPUT$1''_variants_GATK_step1_norm.vcf /home/varma/softwares/annovar/humandb/ -buildver hg19 -out ''$OUTPUT$1''_variants_GATK_step2_annovar -remove -protocol refGene,1000g2014oct_all,snp138,clinvar_20150629,exac03,esp6500siv2_all,ljb23_sift,ljb23_pp2hvar,ljb23_pp2hdiv,ljb23_gerp++ -operation g,f,f,f,f,f,f,f,f,f -nastring . -vcfinput


/home/varma/softwares/annovar/table_annovar.pl ''$OUTPUT$1''_variants_GATK_step1_norm.vcf /home/varma/softwares/annovar/humandb/ -buildver hg19 -out ''$OUTPUT$1''_variants_GATK_step2_annovar -thread 8 -remove -protocol refGene,1000g2015aug_all,avsnp147,clinvar_20170130,exac03,esp6500siv2_all,dbnsfp33a,gerp++gt2,gnomad_exome,gnomad_genome,intervar_20170202,fathmm,gwava,eigen -operation g,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput


######-----CADD--database--annotation-------
#for header
#rm ''$TMP''xa*
#rm ''$TMP''CADD_''$ch''_*

#rm ''$TMP$1''_variants_GATK_step1.vcf



