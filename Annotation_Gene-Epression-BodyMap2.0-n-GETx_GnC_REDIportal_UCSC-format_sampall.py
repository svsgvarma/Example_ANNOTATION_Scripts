#!/usr/bin/python


"""
-Script to add gene expression FPKM scores from Gene-Epression-BodyMap2.0 public resource for 16 different tissues, and GTEx annotation for each gene can be performed

Note: please change the coloumn name according to the exonic and  genename ( (value[24] == 'exonic') and line.strip().split()[25]) 
#[24]==exonic;[25]==genename; #cut -f25,26,30,31   

INPUT and OUTPUT file for the script:

python Annotation_Gene-Epression-BodyMap2.0-n-GETx_GnC_REDIportal_UCSC-format_sampall.py -i HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_Variants_Zyg_PASS_MAF5%_CRLD_RE-RG-LDB.txt -o HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_Variants_Zyg_PASS_MAF5%_CRLD_RE-RG-LDB_GExp_GnC_REDIportal_UCSC.txt -d ./Example_WORKDIR/ -s 3

#Gene_Names	adipose	adrenal	blood	brain	breast	colon	heart	kidney	liver	lung	lymph	ovary	prostate	skeletal_muscle	testes	thyroid

"""


import sys
import subprocess
import re
from argparse import ArgumentParser

#############------------

class fileHandler:
	def __init__(self):
		self.data = []
		#print "Calling fileHandler constructor"
	def open_file(self,readfl):
		self.rfile = open(readfl,'r').readlines()
		return self.rfile
	def write_file(self,writefl):
		self.wfile = open(writefl,'w')
		return self.wfile

class SearchDB(fileHandler):
	def __init__(self):
		self.data = []
		from collections import defaultdict
		self.ident_ranges_HMBM = defaultdict(list)

	def VariantCount_onGenes(self,readfl1,writedir,samp):
		"""
		Calling method to count variants on each gene...
		"""
		#----------------------
		samp2 = int(samp) * 2
		gene_name = 17+int(samp2)+2
		#----------------------
		with open(writedir+readfl1,'r') as f1, open(writedir+"temp-outfile0",'w') as output:
			first_line = f1.readline().strip()
			HeadDB1 = "#Gene based Variants count"
			output.write(first_line+"\t"+HeadDB1+"\n")
			#----------------------
			def Countall(alist):
				dict_count={}
				for i in alist:
					if not dict_count.has_key(i): dict_count[i]=1  #also: if not i in d
					else: dict_count[i]+=1
				return dict_count
			#----------------------
			List_Genes =[]
			for line in f1:
				value = line.strip().split()[gene_name]
				List_Genes.append(value)
			countfun = Countall(List_Genes)
			for line in open(writedir+readfl1,'r').readlines()[1:]:
				lineS = line.strip()
				value = line.strip().split()[gene_name]
				#print value+str("\t")+str(countfun[value])		
				output.write(lineS+"\t"+str(countfun[value])+"\n")
		print "Done counting variants on each gene...."

	def funczip(self,lt1,lt2):
		return '\t'.join(i + ', ' + j for i,j in zip(lt1,lt2))

	def Search_ExpDB1(self, GName):
		self.list1 = []
		DIR="/media/varma/My_Passport_2/Scripts-for-LocalDBsearch/Gene-Expression_BodyMap2.0_GETx/FPKMs_Human-BodyMap2.0_mRNA-seq-data_gene.matrix.csv"
		try:
			True
			cmdFls1 = "LANG=C grep -m 1 -w '"+str(GName)+"' "+str(DIR)+""
			grepout1 =  subprocess.check_output(cmdFls1, shell=True)
			self.list1 = grepout1.strip().split(",")
		except:
			True
			#print GName+" not present in DB"
			self.list1 = [str('"'+GName+'"'), '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
		return self.list1

	def Human_BodyMap20(self,writedir,samp):
		"""
		Calling HBM localsearch method
		"""
		#----------------------
		samp2 = int(samp) * 2
		gene_name = 17+int(samp2)+2
		#----------------------
		with open(writedir+"temp-outfile0",'r') as f1, open(writedir+"temp-outfile1",'w') as output:
			first_line = f1.readline().strip()
			HeadDB1 = "#Gene_Names_(Gene-Epression-BodyMap2.0_16-tissues)	Adipose	Adrenal	Blood	Brain	Breast	Colon	Heart	Kidney	Liver	Lung	Lymph	Ovary	Prostate	Skeletal_muscle	Testes	Thyroid"
			output.write(first_line+"\t"+HeadDB1+"\n")
			for line in f1:
				lineS = line.strip()
				value = line.strip().split()[gene_name]
				val1 = value.split(",")
				if len(val1) == 1:
					ExpG = "\t".join(self.Search_ExpDB1(val1[0]))
					output.write(lineS+"\t"+ExpG+"\n")
				elif len(val1) != 1:
					ExpG1 = self.Search_ExpDB1(val1[0])
					ExpG2 = self.Search_ExpDB1(val1[1])
					MergExpG = self.funczip(ExpG1, ExpG2)
					output.write(lineS+"\t"+MergExpG+"\n")
		print "Done seach for Exp-data of Human_BodyMap2.0...."

	def Search_ExpDB2(self, GName):
		self.list1 = []
		self.mergedlist = []
		DIR="/media/varma/My_Passport_2/Scripts-for-LocalDBsearch/Gene-Expression_BodyMap2.0_GETx/hg19_Expression-GTEx.txt"
		try:
			True
			cmdFls1 = "LANG=C grep -m 1 -w '"+str(GName)+"' "+str(DIR)+""
			grepout1 =  subprocess.check_output(cmdFls1, shell=True)
			self.list1 = re.split('\t|,', grepout1.strip())[:-1]
		except:
			True
			#print GName+" not present in DB"
			self.list1 = ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.','.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.','.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
		return self.list1

	def GeneExp_GETx(self,writedir,samp):
		"""
		Calling GETx localsearch method
		"""
		#----------------------
		samp2 = int(samp) * 2
		gene_name = 17+int(samp2)+2
		#----------------------

		with open(writedir+"temp-outfile1") as f1, open(writedir+"temp-outfile2",'w') as output:
			first_line = f1.readline().strip()
			HeadDB2 = "#chrom_(Gene-Epression-GTEx_53-tissues)	chromStart	chromEnd	name	score	strand	geneId	geneType	expCount	expScores_(-'Tissue','Number of RNASeq and Genotyped samples','Number of RNASeq Samples','Number of eGenes')_'Muscle - Skeletal','361','430','7508'	'Whole Blood','338','393','7332'	'Skin - Sun Exposed (Lower leg)','302','357','9069'	'Adipose - Subcutaneous','298','350','9026'	'Artery - Tibial','285','332','8535'	'Thyroid','278','323','10610'	'Lung','278','320','7774'	'Nerve - Tibial','256','304','10489'	'Esophagus - Mucosa','241','286','7784'	'Cells - Transformed fibroblasts','272','284','9413'	'Skin - Not Sun Exposed (Suprapubic)','196','250','5793'	'Esophagus - Muscularis','218','247','7322'	'Adipose - Visceral (Omentum)','185','227','4667'	'Artery - Aorta','197','224','6599'	'Heart - Left Ventricle','190','218','4814'	'Breast - Mammary Tissue','183','214','4411'	'Colon - Transverse','169','196','4853'	'Heart - Atrial Appendage','159','194','4203'	'Stomach','170','193','3818'	'Testis','157','172','9435'	'Pancreas','149','171','4568'	'Esophagus - Gastroesophageal Junction','127','153','2978'	'Colon - Sigmoid','124','149','3066'	'Adrenal Gland','126','145','3529'	'Artery - Coronary','118','133','2511'	'Brain - Cerebellum','103','125','4528'	'Liver','97','119','1644'	'Cells - EBV-transformed lymphocytes','114','118','3116'	'Brain - Caudate (basal ganglia)','100','117','2612'	'Brain - Cortex','96','114','2768'	'Brain - Nucleus accumbens (basal ganglia)','93','113','2202'	'Brain - Frontal Cortex (BA9)','92','108','2152'	'Prostate','87','106','1500'	'Brain - Cerebellar Hemisphere','89','105','3403'	'Spleen','89','104','2881'	'Pituitary','87','103','2260'	'Brain - Putamen (basal ganglia)','82','97','1653'	'Ovary','85','97','1615'	'Brain - Hypothalamus','81','96','1253'	'Vagina','79','96','859'	'Brain - Hippocampus','81','94','1164'	'Small Intestine - Terminal Ileum','77','88','1415'	'Brain - Anterior cingulate cortex (BA24)','72','84','1289'	'Uterus','70','83','941'	'Brain - Amygdala','62','72','null'	'Brain - Spinal cord (cervical c-1)','59','71','null'	'Brain - Substantia nigra','56','63','null'	'Minor Salivary Gland','51','57','null'	'Kidney - Cortex','26','32','null'	'Bladder','11','11','null'	'Cervix - Ectocervix','6','6','null'	'Fallopian Tube','6','6','null'	'Cervix - Endocervix','5','5','null'"
			
			output.write(first_line+"\t"+HeadDB2+"\n")
			for line in f1:
				lineS = line.strip()
				value = line.strip().split()[gene_name]
				val1 = value.split(",")
				if len(val1) == 1:
					ExpG = "\t".join(self.Search_ExpDB2(val1[0]))
					#print ExpG
					output.write(lineS+"\t"+ExpG+"\n")
				elif len(val1) != 1:
					ExpG1 = self.Search_ExpDB2(val1[0])
					ExpG2 = self.Search_ExpDB2(val1[1])
					MergExpG = self.funczip(ExpG1, ExpG2)
					#print MergExpG
					output.write(lineS+"\t"+MergExpG+"\n")
		print "Done seach for Exp-data of GETx...."

		return None
	def REDIportal_add(self,writedir,samp):
		"""
		Calling REDIportal localsearch method
		"""
		def srchdb3(RR):
			try:
				True
				#print "chr"+RR[0]+":"+RR[8]+"-"+RR[8]
				cmdFls1 = "tabix -f /media/varma/SAMSUNG/REDItools/REDIportal/table1_full.tsv.gz chr"+RR[0]+":"+RR[8]+"-"+RR[8]
				grepout =  subprocess.check_output(cmdFls1, shell=True)
			except:
				False
				grepout = "."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."
			return grepout

		with open(writedir+"temp-outfile2") as f1, open(writedir+"temp-outfile3",'w') as output:
			first_line = f1.readline().strip()
			HeadDB2 = "#Region_#CHROMOSOME NAME	Position_#GENOMIC COORDINATE	Ref_#REFERENCE NUCLEOTIDE	Ed_#EDITED NUCLEOTIDE	Strand_#ORIENTATION	db_#KNOWN IN ATLAS/RADAR/DARNED	type_#EDITING CLASSIFICATION IN ALU/REP/NONREP	dbsnp_#GENOMIC SNP	repeat_#REPEAT NAME	Func.wgEncodeGencodeBasicV19_#GENE LOCATION ACCORDING TO GENCODE V19 BY ANNOVAR WITH VALUES: exonic,splicing,ncRNA,UTR5,UTR3,intronic,upstream,downstream,intergenic	Gene.wgEncodeGencodeBasicV19_#GENE SYMBOL ACCORDING TO GENCODE V19	ExonicFunc.wgEncodeGencodeBasicV19_#EXONIC FUNCTION ACCORDING TO GENCODE V19 BY ANNOVAR WITH VALUES: nonsynonymous SNV, synonymous SNV, stopgain, stoploss	AAChange.wgEncodeGencodeBasicV19_#AMINO ACID CHANGES ACCORDING TO GENCODE V19 BY ANNOVAR	Func.refGene_#GENE LOCATION ACCORDING TO REFSEQ BY ANNOVAR WITH VALUES: exonic,splicing,ncRNA,UTR5,UTR3,intronic,upstream,downstream,intergenic	Gene.refGene_#GENE SYMBOL ACCORDING TO REFSEQ	ExonicFunc.refGene_#EXONIC FUNCTION ACCORDING TO REFSEQ BY ANNOVAR WITH VALUES: nonsynonymous SNV, synonymous SNV, stopgain, stoploss	AAChange.refGene_#AMINO ACID CHANGES ACCORDING TO REFSEQ BY ANNOVAR	Func.knownGene_#GENE LOCATION ACCORDING TO UCSC BY ANNOVAR WITH VALUES: exonic,splicing,ncRNA,UTR5,UTR3,intronic,upstream,downstream,intergenic	Gene.knownGene_#GENE SYMBOL ACCORDING TO UCSC	ExonicFunc.knownGene_#EXONIC FUNCTION ACCORDING TO UCSC BY ANNOVAR WITH VALUES: nonsynonymous SNV, synonymous SNV, stopgain, stoploss	AAChange.knownGene_#AMINO ACID CHANGES ACCORDING TO REFSEQ BY ANNOVAR	phastConsElements46way_#PHASTCON CONSERVATION ACROSS 46 ORGANISMS FROM UCSC (LOD/SCORE)	nSamples_#NUMBER OF EDITED SAMPLES	nTissues_#NUMBER OF EDITED TISSUES	nBodySites_#NUMBER OF EDITED BODY SITES"
			
			output.write(first_line+"\t"+HeadDB2+"\n")
			for g in f1:
				g1 = g.strip()
				gg = g1.split("\t")
				SCOREout = srchdb3(gg).strip()
				CDS = SCOREout.split()
				if len(CDS)==0:
					output.write(str(g1+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\n"))
				else:
					output.write(str(g1+"\t"+SCOREout+"\n"))
		print "Done adding REDIportal...."
		return None

	def UCSC_add(self,writefl1,writedir,samp):
		"""
		Convert to UCSC BED6 format file...
		###------
		#chrom	chromStart	chromEnd	name	score	strand	thickStart	thickEnd	itemRgb
		#(value[24] == 'exonic') # cut cut -f25,30,31
		#annovar_loc(Func.refGene)	1000g2015aug_all	ExAC_ALL
		#Note: Change this if the value is int---#Decimal(value[36]) <= "0.01")
		"""
		#----------------------
		samp2 = int(samp) * 2
		loc = 17+int(samp2)+1
		KG = 17+int(samp2)+6
		dbSNP = 17+int(samp2)+7
		#----------------------
		with open(writedir+"temp-outfile3",'r') as f1, open(writedir+writefl1,'w') as output:
			first_line = f1.readline().strip()
			HeadDB1 = "#UCSC-format-header_(track name='WGS_ADCA' description='ADCA_009-028-080' visibility=2 itemRgb=On)_chrom	chromStart	chromEnd	name	score	strand	thickStart	thickEnd	itemRgb"
			output.write(HeadDB1+"\t"+first_line+"\n")

			#flsp = first_line.split("\t")
			#print flsp[24],flsp[25],flsp[29],flsp[30]
			
			for line in f1:
				val1 = line.strip()
				value = val1.split("\t")
				if ((value[loc] == 'exonic')  and (value[KG] == ".") and (value[dbSNP] == ".")): 
					k1 = str("chr"+value[0]+"\t"+str(int(value[1])-1)+"\t"+value[2]+"\t"+"Pos1"+"\t"+value[14]+"\t"+"+"+"\t"+str(int(value[1])-1)+"\t"+value[2]+"\t"+"255,0,0")
					output.write(k1+"\t"+val1+"\n")
				elif (value[loc] == 'exonic'):
					k1 = str("chr"+value[0]+"\t"+str(int(value[1])-1)+"\t"+value[2]+"\t"+"Pos1"+"\t"+value[14]+"\t"+"+"+"\t"+str(int(value[1])-1)+"\t"+value[2]+"\t"+"255,0,55")
					output.write(k1+"\t"+val1+"\n")
				###----non-exonic -----
				elif ((value[loc] != 'exonic' ) and (value[KG] == ".") and (value[dbSNP] == ".")):
					k1 = str("chr"+value[0]+"\t"+str(int(value[1])-1)+"\t"+value[2]+"\t"+"Pos1"+"\t"+value[14]+"\t"+"+"+"\t"+str(int(value[1])-1)+"\t"+value[2]+"\t"+"0,0,255")
					output.write(k1+"\t"+val1+"\n")
				###----Rest -----
				else:
					k1 = str("chr"+value[0]+"\t"+str(int(value[1])-1)+"\t"+value[2]+"\t"+"Pos1"+"\t"+value[14]+"\t"+"+"+"\t"+str(int(value[1])-1)+"\t"+value[2]+"\t"+"0,55,255")
					output.write(k1+"\t"+val1+"\n")
		subprocess.check_output("rm "+writedir+"temp-outfile0", shell=True)
		subprocess.check_output("rm "+writedir+"temp-outfile1", shell=True)
		subprocess.check_output("rm "+writedir+"temp-outfile2", shell=True)
		subprocess.check_output("rm "+writedir+"temp-outfile3", shell=True)
		print "Done adding UCSC format...."
		return None

# file1: Input positions SNP
# write1: Output file

if __name__ == "__main__":
	parser = ArgumentParser(
	description="-Script to add gene expression FPKM scores from Gene-Epression-BodyMap2.0 public resource for 16 different tissues, and GTEx annotation for each gene can be performed.")

	parser.add_argument("-i", "--infile", dest="INPUTFile",default=False, required=True, help="Full input file name" )
	parser.add_argument("-o", "--outfile", dest="OUTPUTFile",default=False, required=True, help="Full output file name" )
	parser.add_argument("-d", "--dir", dest="Directory", default="False", required=True, help="Working directory; input and output file location full Path" )
	parser.add_argument("-s", "--samples", dest="SampleNumber",default=False, required=True, help="Number of samples included, ex: 3" )

	if len(sys.argv) == 2:
		parser.print_help()
		print "Example: ./Annotation_Gene-Epression-BodyMap2.0-n-GETx_GnC_REDIportal_UCSC-format_sampall.py -i HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_Variants_Zyg_PASS_MAF5%_CRLD_RE-RG-LDB.txt -o HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_Variants_Zyg_PASS_MAF5%_CRLD_RE-RG-LDB_GExp_GnC_REDIportal_UCSC.txt -d /WORKDIR/FULLPATH/ -s 3"
		sys.exit()

	args = parser.parse_args()

	INF = args.INPUTFile
	OUTF = args.OUTPUTFile
	DIR = args.Directory
	SAMP = args.SampleNumber

clF0 = SearchDB().VariantCount_onGenes(INF,DIR,SAMP)
clF1 = SearchDB().Human_BodyMap20(DIR,SAMP)
clF2 = SearchDB().GeneExp_GETx(DIR,SAMP)
clF3 = SearchDB().REDIportal_add(DIR,SAMP)
clF3 = SearchDB().UCSC_add(OUTF,DIR,SAMP)


