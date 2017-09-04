#!/bin/bash


###-----

#./runANNOTATIONS_all-samples_Example-DCM.bash > ./runANNOTATIONS_all-samples_Example-DCM.bash.log 2>&1

###-----Parse WGS all-----
########################
###---

python ./Pars_ANNOVAR-outfile_allsamp_include-zygo.py -fid HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT -sid B00H7EW,B00H7EX,B00H7EY -f HET,HET,HOM -d ./Example_WORKDIR/

echo "Done running script Parse from ANNOVAR output files...."

###---
########################
###---DCM--3 samples annotation----

python ./Annotation_withlocal-CADD-ReMM-LINSIGHT_sampall.py -i HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_Variants_Zyg_PASS_MAF5%.txt -o HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_Variants_Zyg_PASS_MAF5%_CRLD.txt -d ./Example_WORKDIR/

echo "Done running script... CADD-ReMM-LINSIGHT-DANN...."

python ./Annotation_checkfor_RE_RG_LDB-WES-WGS_sampall_CMH.py -i HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_Variants_Zyg_PASS_MAF5%_CRLD.txt -o HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_Variants_Zyg_PASS_MAF5%_CRLD_RE-RG-LDB.txt -d ./Example_WORKDIR/

echo "Done running script... RE_RG_LDB-WES-WGS...."

python ./Annotation_Gene-Epression-BodyMap2.0-n-GETx_GnC_REDIportal_UCSC-format_sampall.py -i HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_Variants_Zyg_PASS_MAF5%_CRLD_RE-RG-LDB.txt -o HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_Variants_Zyg_PASS_MAF5%_CRLD_RE-RG-LDB_GExp_GnC_REDIportal_UCSC.txt -d ./Example_WORKDIR/ -s 3

echo "Done running script... Gene-Epression-BodyMap2.0-n-GETx_UCSC-format...."

###---
#######################
###---Protein codon changes-----

python ./Annotation_Protein-codon-changes.py -i HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_variants_GATK_step2_annovar.hg19_multianno.vcf -o HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_variants_GATK_Protein-CODON_changes.txt -d ./Example_WORKDIR/

echo "Done running script... Protein codon changes...."
