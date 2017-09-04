#!/usr/bin/python


"""
WholeGenomeSequencing_Repeat-Mask-Annotation

-Whole Genome Sequcing (WGS) of individual sample VCF file or multi-sample VCF file can be recruted for annotation pipleline, this allows to check for varaints falling with in the low-complexity regions

Note: please change the coloumn name according to the exonic (infll1[24] =="exonic") # cut -f25

INPUT_FILES:
python Annotation_checkfor_RE_RG_LDB-WES-WGS_sampall_CMH.py -i HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_Variants_Zyg_PASS_MAF5%_CRLD.txt -o HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_Variants_Zyg_PASS_MAF5%_CRLD_RE-RG-LDB.txt -d ./Example_WORKDIR/

"""


import sys
import subprocess
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
		self.ident_ranges_UCSC = defaultdict(list)
		self.ident_ranges_RG = defaultdict(list)
		self.ident_ranges_RG_LCR_HD = defaultdict(list)
		self.ident_ranges_RG_UM75 = defaultdict(list)
		
	def Lsearch_UCSC_RG_LCR_HD_UM75(self,readfl1,writedir):

		"""
		Calling UCSC localsearch method:
		-Script to filter out false negatives and true negatives of variants with in the Repeat-Masking(RM) regions, RM data was obtained from UCSC browser
		Calling RG localsearch method:
		-Checks and retain reliable variants which are falling within concordant regions of public resources dowloaded for ~200 WGS data
		"""

		DIR="/media/varma/My_Passport_2/Scripts-for-LocalDBsearch/RG-ReliableGenome/wtchg-rg-master/data/public/"
		with open("/media/varma/My_Passport_2/Scripts-for-LocalDBsearch/RG-ReliableGenome/Repeat-masking/rmsk-pos.txt",'r') as f2:
			for row in f2:
				ident, start, end = row.strip().split()
				start, end = int(start), int(end)
				self.ident_ranges_UCSC[ident[3:]].append((start, end))
		
		with open(DIR+"20160422_RG-win1000-score1_-3-RELIABLE-above0.5.bed",'r') as f2:
			for row in f2:
				ident, start, end = row.strip().split("\t")[0:3]
				start, end = int(start), int(end)
				self.ident_ranges_RG[ident].append((start, end))
		with open(DIR+"20160422_RG-win1000-score1_-3-RELIABLE-above0.5.bed.gz-min-LCR-min-HD.bed",'r') as f2:
			for row in f2:
				ident, start, end = row.strip().split("\t")[0:3]
				start, end = int(start), int(end)
				self.ident_ranges_RG_LCR_HD[ident].append((start, end))
		with open(DIR+"20160825_RG-win1000-score1_-3-RELIABLE-above0.5.-min-UM75.bed",'r') as f2:
			for row in f2:
				ident, start, end = row.strip().split("\t")[0:3]
				start, end = int(start), int(end)
				self.ident_ranges_RG_UM75[ident].append((start, end))
		
		with open(writedir+readfl1,'r') as f1, open(writedir+"temp-outfile1_UCSC_RG",'w') as output:
			first_line = f1.readline().strip()
			output.write(first_line+str("\t"+"RE_Repeating-Elements-by-RepeatMasker"+"\t"+"RG_ReliableGenome"+"\t"+"RG_ReliableGenome_min-LCR-min-HD"+"\t"+"RG_ReliableGenome_min-UM75")+"\n")
			for line in f1:
				Lstcat= []
				line1 = line.strip()
				ident, value = line.strip().split()[0:2]
				value = int(value)
				if not any(start <= value <= end for start, end in self.ident_ranges_UCSC[ident]):
					True
					Lstcat.append("NotIn_RE")
				else:
					True
					Lstcat.append("In_RE")
				if any(start <= value <= end for start, end in self.ident_ranges_RG[ident]):
					True
					Lstcat.append("RG")
				else:
					True
					Lstcat.append(".")
				if any(start <= value <= end for start, end in self.ident_ranges_RG_LCR_HD[ident]):
					True
					Lstcat.append("RG_min-LCR-HD")
				else:
					True
					Lstcat.append(".")
				if any(start <= value <= end for start, end in self.ident_ranges_RG_UM75[ident]):
					True
					Lstcat.append("RG_min-UM75")
				else:
					True
					Lstcat.append(".")
				#print str(str(value)+"\t"+Lstcat[0]+"\t"+Lstcat[1]+"\t"+Lstcat[2]+"\t"+Lstcat[3])
				output.write(str(line1+"\t"+Lstcat[0]+"\t"+Lstcat[1]+"\t"+Lstcat[2]+"\t"+Lstcat[3]+"\n"))
		output.close()
		print "Done seach for UCSC_RG_LCR_HD_UM75...."
		return None

	def Lsearch_WGS_WES(self,writefl1,writedir):
		"""
		Calling localsearch of  WGS and WES method
		-Final check against local 22 WGS and ~ 93 WES individual datasets for futher filtering step
		"""
		infl1 = self.open_file(writedir+"temp-outfile1_UCSC_RG") 
		infl_nhr = infl1[1:]
		file2 = self.write_file(writedir+writefl1)
		file2.write(str(infl1[0].strip())+str("\t"+"LocalDataBase(LDB)-search_WGS_WES")+"\n")
		totlns= len(infl_nhr)
		count1 = 0;
		while (count1 < totlns):
			infll_st = infl_nhr[count1].strip()
			infll1 = infl_nhr[count1].strip().split()
			#grep -w '1'$'\t''15118'
			def Lsearch_just_WGS(infll_st):
				###-------------
				#print "not present in WES"
				infll = infll_st.split()
				try:
					True
					cmdFls2 = "tabix -f /media/varma/My_Passport_1/Backup_Varma/RenduOsler1/39_WGS_VCF/39_12345_1-22_variants_GATK_step1_norm.vcf.gz "+infll[0]+":"+infll[1]+"-"+infll[1]
					grepout2 =  subprocess.check_output(cmdFls2, shell=True)
					if len(grepout2.split())==0:
						raise ValueError('empty string')
					else:
						file2.write(str(infll_st)+str("\t"+"In_LDB-WGS-H1")+"\n")
				except ValueError:
					True
					###-------------
					#print "not present in 39"
					try:
						True
						cmdFls3 = "tabix -f /media/varma/My_Passport_1/Backup_Varma/RenduOsler2/V6_WGS_VCF/V6_12345_1-22_variants_GATK_step1_norm.vcf.gz "+infll[0]+":"+infll[1]+"-"+infll[1]
						grepout3 =  subprocess.check_output(cmdFls3, shell=True)
						if len(grepout3.split())==0:
							raise ValueError('empty string')
						else:
							file2.write(str(infll_st)+str("\t"+"In_LDB-WGS-H2")+"\n")
					except ValueError:
						True
						###-------------
						#print "not present in 39 & V6_DBs"
						try:
							True
							cmdFls4 = "tabix -f /media/varma/SAMSUNG/HCM_WGS_VCF/HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_variants_SNPnINDELcalls.vcf.gz "+infll[0]+":"+infll[1]+"-"+infll[1]
							grepout4 =  subprocess.check_output(cmdFls4, shell=True)
							if len(grepout4.split())==0:
								raise ValueError('empty string')
							else:
								file2.write(str(infll_st)+str("\t"+"In_LDB-WGS-DCM")+"\n")
						except ValueError:
							True
							###-------------
							#print "not present in 39, V6_DBs and ADCA-F1"
							try:
								True
								cmdFls4 = "tabix -f /media/varma/Maxtor/ADCA_WGS_VCF/ADCA_009-028-080_chr1-MT_variants_SNPnINDELcalls.vcf.gz "+infll[0]+":"+infll[1]+"-"+infll[1]
								grepout4 =  subprocess.check_output(cmdFls4, shell=True)
								if len(grepout4.split())==0:
									raise ValueError('empty string')
								else:
									file2.write(str(infll_st)+str("\t"+"In_LDB-WGS-ADCA1")+"\n")
							except ValueError:
								True
								try:
									True
									cmdFls5 = "for i in /media/varma/Maxtor/MALRARWGSTEV/variants/B00HP0*_hs37d5_BOTH.merged.annot.vcf.gz; do tabix $i '"+str("chr"+infll[0]+":"+infll[1]+"-"+infll[1])+"' && break; done"
									grepout5 =  subprocess.check_output(cmdFls5, shell=True)
									if len(grepout5.split())==0:
										raise ValueError('empty string')
									else:
										file2.write(str(infll_st)+str("\t"+"In_LDB-WGS-CA2")+"\n")
								except ValueError:
									True
									###-------------
									#print "not present in 39, V6_DBs, ADCA1-F-3IND and ADCA2-F-9IND"
									file2.write(str(infll_st)+str("\t"+"NotIn_LDB-WGS")+"\n")

				return None
			def Lsearch_just_WES(infll_st):
				###-------------
				#print "not present in WES"
				infll = infll_st.split()
				try:
					cmdFls1 = "for i in /media/varma/My_Passport_2/SCA_WES_VCF/*_SNPnINDELcalls.vcf; do LANG=C grep -m 1 -w '"+str(infll[0]+'\t'+infll[1])+"' $i && break; done"
					grepout1 =  subprocess.check_output(cmdFls1, shell=True)
					file2.write(str(infll_st)+str("\t"+"In_LDB-WES_SCA")+"\n")
				except:
					True
					try:
						cmdFls2 = "for i in /media/varma/My_Passport_2/DDLD_WES_VCF/*_SNPnINDELcalls.vcf; do LANG=C grep -m 1 -w '"+str(infll[0]+'\t'+infll[1])+"' $i && break; done"
						grepout2 =  subprocess.check_output(cmdFls2, shell=True)
						file2.write(str(infll_st)+str("\t"+"In_LDB-WES_DDLD")+"\n")
					except:
						True
						try:
							cmdFls3 = "for i in /media/varma/My_Passport_2/DVT*_WES_VCF/*_SNPnINDELcalls.vcf; do LANG=C grep -m 1 -w '"+str(infll[0]+'\t'+infll[1])+"' $i && break; done"
							grepout3 =  subprocess.check_output(cmdFls3, shell=True)
							file2.write(str(infll_st)+str("\t"+"In_LDB-WES_DVT3")+"\n")
						except:
							True
							try:
								cmdFls4 = "for i in /media/varma/My_Passport_2/CMT_WES_VCF/*_SNPnINDELcalls.vcf; do LANG=C grep -m 1 -w '"+str(infll[0]+'\t'+infll[1])+"' $i && break; done"
								grepout4 =  subprocess.check_output(cmdFls4, shell=True)
								file2.write(str(infll_st)+str("\t"+"In_LDB-WES_CMT")+"\n")
							except:
								True
								WGSDB = Lsearch_just_WGS(infll_st)
				return None

			if infll1[24] =="exonic":
				True
				WESDB = Lsearch_just_WES(infll_st)
			else:
				True
				WGSDB = Lsearch_just_WGS(infll_st)
			count1 += 1
		file2.close()
		subprocess.check_output("rm "+writedir+"temp-outfile1_UCSC_RG", shell=True)
		return None

# file1: Input positions SNP
# write1: Output file 

if __name__ == "__main__":
	parser = ArgumentParser(
	description="Script to annotate with 1) WholeGenomeSequencing_Reliable Genome(RG) Annotation 2) UCSC Repeat Elements(RE) and 3) Local Database of Whole genome and Whole Exome search(LDB-WES-WGS).")

	parser.add_argument("-i", "--infile", dest="INPUTFile",default=False, required=True, help="Full input file name" )
	parser.add_argument("-o", "--outfile", dest="OUTPUTFile",default=False, required=True, help="Full output file name" )
	parser.add_argument("-d", "--dir", dest="Directory", default="False", required=True, help="Working directory; input and output file location full Path" )

	if len(sys.argv) == 2:
		parser.print_help()
		print "Example: ./Annotation_checkfor_RE_RG_LDB-WES-WGS_sampall_CMH.py -i HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_Variants_Zyg_PASS_MAF5%_CRLD.txt -o HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_Variants_Zyg_PASS_MAF5%_CRLD_RE-RG-LDB.txt -d /WORKDIR/FULLPATH/"
		sys.exit()

	args = parser.parse_args()

	INF = args.INPUTFile
	OUTF = args.OUTPUTFile
	DIR = args.Directory

	clF1 = SearchDB().Lsearch_UCSC_RG_LCR_HD_UM75(INF,DIR)
	SDB = SearchDB().Lsearch_WGS_WES(OUTF,DIR)




