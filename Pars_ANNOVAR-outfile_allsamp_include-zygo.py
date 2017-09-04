#!/usr/bin/python

"""
#Script to Pars vcf...
#####------Inputs-------
# script.py INPUT1=sampleIDs INPUT2=orderedsampleIDs_withcommasep(,) INPUT3=Filteringtype(HET or HOM)_withcommasep(,) INPUT4=PATH-TO-WORK

# please change the Het and Numbers of samples included 

Example: 
./Pars_ANNOVAR-outfile_allsamp_include-zygo.py -fid HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT -sid B00H7EW,B00H7EX,B00H7EY -f HET,HET,HOM -d ./Example_WORKDIR/

##----Based on these below ANNOVAR reference settings----
 
/home/varma/softwares/annovar/table_annovar.pl ''$OUTPUT$1''_variants_GATK_step1_norm.vcf /home/varma/softwares/annovar/humandb/ -buildver hg19 -out ''$OUTPUT$1''_variants_GATK_step2_annovar -thread 8 -remove -protocol refGene,1000g2015aug_all,avsnp147,clinvar_20170130,exac03,esp6500siv2_all,dbnsfp33a,gerp++gt2,gnomad_exome,gnomad_genome,intervar_20170202,fathmm,gwava,eigen -operation g,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput

"""

import sys
import re
import subprocess
from subprocess import *
from decimal import *

from argparse import ArgumentParser


class SearchDB:
	def __init__(self):
		self.data = []
	def ANNOVAR_pars(self,c,b1,HETHOM,writefl1):
		"""
		Parse ANNOVAR file method
		"""
		#------
		bb= list(map(str, b1.split(',')))
		def magic(numList):
			s = '-'.join(map(str, numList))
			return str(s)
		c2 = str(magic(bb))
		lst1=[];
		for xx in range(len(bb)):
			lst1.append("Sample_"+str(bb[xx])+"_genotype"+str("	")+"Sample_"+str(bb[xx])+"_AD:DP:GQ:PGT:PID:PL:TP")
		Joinall_gvcfs1 = '\t'.join(lst1)
		#------
		wfile1 = open(writefl1+c+"_Variants.txt",'w')
		wfile2 = open(writefl1+c+"_Variants_Zyg_MAF5%.txt",'w')
		wfile3 = open(writefl1+c+"_Variants_Zyg_PASS_MAF5%.txt",'w')
		
		with open(writefl1+c+"_variants_GATK_step2_annovar.hg19_multianno.txt",'r') as f1:
			first_line = f1.readline().strip()
			Header_title = str("annovar_Chr	annovar_Start-pos	annovar_End-pos	annovar_Ref('-' for insertion)	annovar_Alt('-' for deletion)	allele_freq	nb_homoz	nb_heteroz	POS	REF	ALT	QUAL	FILTER	DP(Approximate read depth)	MQ(Mean Mapping Quality)	QD(Qual By Depth)	VQSLOD	MLEAF	")+Joinall_gvcfs1+str("	annovar_loc(Func.refGene)	annovar_gene_name(Gene.refGene)	GeneDetail.refGene	annovar_effect(ExonicFunc.refGene)	annovar_exonic_effects(AAChange.refGene)	1000g2015aug_all	dbSNP147(avsnp147)	Clinvar-Significance(CLINSIG)	Clinvar-disease-name(CLNDBN)	Clinvar-Accession-and-Versions(CLNACC)	Clinvar-disease-database-name(CLNDSDB)	Clinvar-disease-database-ID(CLNDSDBID)	ExAC_ALL	gnomAD_exome_ALL	gnomAD_genome_ALL	ESP(esp6500siv2_all)	SIFT_score	SIFT_pred	Polyphen2_HVAR_score	Polyphen2_HVAR_pred	LRT_score	LRT_pred	MutationTaster_score	MutationTaster_pred	PROVEAN_score	PROVEAN_pred	VEST3_score	MetaSVM_score	MetaSVM_pred	DANN_score	integrated_fitCons_score	integrated_confidence_value	GERP++_annotation(Higher the score,the more conserved the site)gerp++gt2	InterVar(automated)	FATHMM_noncoding	FATHMM_coding	GWAVA_region_score	GWAVA_tss_score	GWAVA_unmatched_score	Eigen")
			wfile3.write(Header_title+"\n")
			wfile2.write(Header_title+"\n")
			wfile1.write(Header_title+"\n")
			count1 = 0; count2 = 0; count3 = 0;
			for line in f1:
				elm = line.strip()
				elm1 = elm.split("\t")
				hom = 0; het = 0;
				elm2 = elm1[-len(bb):]
				el1 =[]
				for zyg in elm2:
					el = zyg[0:3]
					el1.append(el)
					#----count hom and het
					el2 = re.compile("[|/]").split(el)
					if el2[0] == "0" and el2[1] == "0":
						False
					elif el2[0] == el2[1]:
						hom+=1
					elif el2[0] != el2[1]:
						het+=1
				lst2=[];lst3=[];lst4=[];lst22=[];lst33=[]
				for xx in range(len(bb)):
					####---Heterozygous-check----
					lst2.append(re.compile("[|/]").split(el1[xx])[0] != re.compile("[|/]").split(el1[xx])[1])
					lst3.append(elm2[xx][0:3]+str("	")+elm2[xx][4:])
					####---Homozygous-check----
					lst22.append((re.compile("[|/]").split(el1[xx])[0] == "1") and (re.compile("[|/]").split(el1[xx])[1] == "1"))
					lst33.append((re.compile("[|/]").split(el1[xx])[0] == "0") and (re.compile("[|/]").split(el1[xx])[1] == "0"))
					Joinall_gvcfs3 = '\t'.join(lst3)

				#-----file2-info-------
				dit1={}
				glm2 = elm1[-int(len(bb)+int(12)):]
				glm3 = '\t'.join(glm2[0:1])
				glm4 = '\t'.join(glm2[4:5])
				glm5 = '\t'.join(glm2[6:10])
				glm6 = elm1[-int(len(bb)+int(2))]
				for fll in glm6.split(";"):
					fl = fll.split("=")
					if len(fl) == 1: #fl[0] == "NEGATIVE_TRAIN_SITE" or fl[0] == "ALLELE_END":
						False
					else:
						dit1[fl[0]] = fl[1]
				if "QD" not in dit1:
					dit1['QD'] = "NA"
				if "VQSLOD" not in dit1:
					dit1['VQSLOD'] = "NA"
					
				#-----file3-info-------
				glm7 = '\t'.join(elm1[5:18])
				glm8 = '\t'.join(elm1[139:145])
				glm9 = '\t'.join(glm2[9:10])

				#----check het-hom in all-------------------
				writeout = str(str(elm1[0])+"\t"+str(elm1[1])+"\t"+str(elm1[2])+"\t"+str(elm1[3])+"\t"+str(elm1[4])+"\t"+str(glm3)+"\t"+str(hom)+"\t"+str(het)+"\t"+str(glm4)+"\t"+str(glm5)+"\t"+str(dit1['DP'])+"\t"+str(dit1['MQ'])+"\t"+str(dit1['QD'])+"\t"+str(dit1['VQSLOD'])+"\t"+str(dit1['MLEAF'])+"\t"+Joinall_gvcfs3.rstrip('\n')+"\t"+str(glm7)+"\t"+str(elm1[93])+"\t"+str(elm1[102])+"\t"+str(elm1[25])+"\t"+str(elm1[26])+"\t"+str(elm1[28])+"\t"+str(elm1[32])+"\t"+str(elm1[34])+"\t"+str(elm1[35])+"\t"+str(elm1[37])+"\t"+str(elm1[38])+"\t"+str(elm1[40])+"\t"+str(elm1[47])+"\t"+str(elm1[49])+"\t"+str(elm1[50])+"\t"+str(elm1[52])+"\t"+str(elm1[54])+"\t"+str(elm1[64])+"\t"+str(elm1[74])+"\t"+str(elm1[76])+"\t"+str(elm1[92])+"\t"+str(elm1[110])+"\t"+str(glm8)+"\n")
				
				def convint_repdot(List):
					True
					List1 = []
					for i in range(0,len(List)):
						if List[i] == '.':
							List1.append(List[i].replace(".","0"))
						else:
							List1.append(List[i])
					List2 = [float(n) for n in List1]
					return List2

				mafvar = []
				mafvar = elm1[10],elm1[17],elm1[93],elm1[102]
				mafcon = convint_repdot(mafvar)

				##-Take Heterozyous (HET) and Homozygous (HOM) input and process; ex:HET,HOM-------
				hlst = []; ht=0
				H1 = HETHOM.split(",")
				for hh in H1:
					if hh == "HET":
						hlst.append(lst2[ht])
					if hh == "HOM":
						hlst.append(lst33[ht])
					ht+=1
				#--------
				if ( all([x == True for x in hlst]) and (glm9 == "PASS") and ((Decimal(mafcon[0]) <= Decimal(0.05)) or (Decimal(mafcon[0]) == Decimal(0))) and ((Decimal(mafcon[1]) <= Decimal(0.05)) or (Decimal(mafcon[1]) == Decimal(0))) and ((Decimal(mafcon[2]) <= Decimal(0.05)) or (Decimal(mafcon[2]) == Decimal(0))) and ((Decimal(mafcon[3]) <= Decimal(0.05)) or (Decimal(mafcon[3]) == Decimal(0))) ):
					True
					wfile3.write(writeout)
					count3 = count3 + 1
				#--------
				if ( all([x == True for x in hlst]) and ((Decimal(mafcon[0]) <= Decimal(0.05)) or (Decimal(mafcon[0]) == Decimal(0))) and ((Decimal(mafcon[1]) <= Decimal(0.05)) or (Decimal(mafcon[1]) == Decimal(0))) and ((Decimal(mafcon[2]) <= Decimal(0.05)) or (Decimal(mafcon[2]) == Decimal(0))) and ((Decimal(mafcon[3]) <= Decimal(0.05)) or (Decimal(mafcon[3]) == Decimal(0))) ):
					True
					wfile2.write(writeout)
					count2 = count2 + 1
				wfile1.write(writeout)
				count1 = count1 + 1
				#print mafcon, writeout
			wfile3.close()
			wfile2.close()
			wfile1.close()

		print "Total_Number of all lines","\t",count1
		print "Total_Number of het filter lines","\t",count2
		print "Total_Number of het and MAF5% fileter lines","\t",count3

		print "Done seach for Exp-data of Human_BodyMap2.0...."

		return None

if __name__ == "__main__":
	parser = ArgumentParser(
	description="Script to pars the annovar .vcf and .txt files and create three output files, one raw file without filters(*_Variants.txt) and two filtered files with 5% Minor Allele Frequency and GATK PASS filter (*_Variants_Zyg_MAF5%.txt and *_Variants_Zyg_PASS_MAF5%.txt). Provide the sample IDs and their zygosity filter (HET or HOM) for each sample respectively.")

	parser.add_argument("-fid", "--fileid", dest="FileID",default=False, required=True, help="INPUT1= SampleID (HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT)" )
	parser.add_argument("-sid", "--sampleid", dest="SampleID",default=False, required=True, help="INPUT2= More than one sampleIDs should be ordered and separated with comma (B00H7EW,B00H7EX,B00H7EY)" )
	parser.add_argument("-f", "--filter", dest="Filtertype",default=False, required=True, help="INPUT3= Filtering the type of Zygosity (HET or HOM) should be comma separated (HET,HET,HOM)" )
	parser.add_argument("-d", "--dir", dest="Directory", default="False", required=True, help="INPUT4= Working directory and input file location full Path" )

	if len(sys.argv) == 2:
		parser.print_help()
		print "Example: ./Pars_ANNOVAR-outfile_allsamp_include-zygo.py -fid HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT -sid B00H7EW,B00H7EX,B00H7EY -f HET,HET,HOM -d /WORKDIR/FULLPATH/"
		sys.exit()

	args = parser.parse_args()

	FID = args.FileID
	SID = args.SampleID
	FILT = args.Filtertype
	DIR = args.Directory

	clF1 = SearchDB().ANNOVAR_pars(FID,SID,FILT,DIR)
	

