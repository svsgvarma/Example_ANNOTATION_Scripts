#!/usr/bin/python


"""
#Script to check for Protein codon changes:

Annotating with correct coding protein (with mutated genetic codon) is needed to predict the function of the gene. When two or three consecutive variants falling with in the coding region (including synonymous and non-synonymous variants), this might change the coding AA. While looking at the combined effect, then the coding AA is not been properly annotated with traditional ANNOVAR tool, hence we tackle this problem using in-house scripts.
(Output file extension: *_Protein-CODON_changes.txt)

#####------Inputs-------
# Annotation_Protein-codon-changes.py INPUT1=FULLPATH+Inputfile(VCF) INPUT1=Outputfile1 INPUT3=Outputfile2 INPUT4=Working-Directory-fullpath/

# python ./Annotation_Protein-codon-changes.py -i HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_variants_GATK_step2_annovar.hg19_multianno.vcf -o HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_variants_GATK_Protein-CODON_changes.txt -d ./Example_WORKDIR/

"""

import sys
import re
import os
import tempfile
import commands
import subprocess
#import subprocess32
from subprocess import *
from subprocess import call
from argparse import ArgumentParser

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
		
	def Search_CODON(self,readfl1,workdir):
		"""
		Calling paried variants of protein coding
		"""
		GATK="/home/varma/softwares/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar"
		SOFT="/home/varma/softwares/"
		HG19="/home/varma/softwares/human_g1k_v37.fasta"
		SNPEFF="/home/varma/softwares/snpEff/snpEff.jar"
		WORK= workdir

		#######################

		cmd1 = "java -Xmx12g -jar "+GATK+" -T VariantFiltration -R "+HG19+" -V "+WORK+readfl1+" -o "+WORK+"snp_clsters2_ws3.vcf --clusterSize 2 --clusterWindowSize 3"
		
		cmd2 = "java -Xmx12g -jar "+GATK+" -T SelectVariants -R "+HG19+" -V "+WORK+"snp_clsters2_ws3.vcf -o "+WORK+"snp_clsters2_ws3_clstronly.vcf -select 'FILTER == SnpCluster'"
		
		cmd3 = "java -Xmx12g -jar "+SNPEFF+" -v GRCh37.75 -formatEff -lof -classic "+WORK+"snp_clsters2_ws3_clstronly.vcf > "+WORK+"snp_clsters2_ws3_clstronly_annt.vcf"
		
		cmd4 = "java -Xmx12g -jar "+GATK+" -T SelectVariants -R "+HG19+" -V "+WORK+"snp_clsters2_ws3_clstronly_annt.vcf -o "+WORK+"snp_clsters2_ws3_clstronly_annt_snv.vcf --selectTypeToExclude INDEL"
		
		cmd5 = "java -Xmx12g -jar "+GATK+" -T VariantsToTable -R "+HG19+" -V "+WORK+"snp_clsters2_ws3_clstronly_annt_snv.vcf -F CHROM -F POS -F REF -F ALT -F avsnp147 -F Func.refGene -F ExonicFunc.refGene -F Gene.refGene -F EFF -F 1000g2015aug_all -F ExAC_ALL -F gnomAD_exome_ALL -F gnomAD_genome_ALL -GF GT --allowMissingData --showFiltered -o "+WORK+"snp_clsters2_ws3_clstronly_annt_snv_clstronly.table"
		
		subprocess.check_output(cmd1, shell=True)
		subprocess.check_output(cmd2, shell=True)
		subprocess.check_output(cmd3, shell=True)
		subprocess.check_output(cmd4, shell=True)
		subprocess.check_output(cmd5, shell=True)
		####--------------------------------------
		def Change_zygo(ref, alt, zyg):
			"""
			program to convert zygosity code ref/alt to 0/1.
			Input Variables:
			
			ref = "A";alt = "G"
			zyg =  ['A/G', 'G/A', 'G/G', 'A/A', './.', 'G/.', './G', './A', 'A/.']
			chzyg = ['0/1', '1/0', '1/1', '0/0', './.', '1/.', './1', './0', '0/.']
			InputUsage Variables:
			chgz = Change_zygo(ref, alt, zyg)

			"""
			import re
			chg_zyg = {};i=1
			for cs in zyg:
				csp = re.split('[|/]',cs)
				if ((ref == csp[0]) and (ref == csp[1]) and ((csp[0] != ".") and (csp[1] != "."))):
					chg_zyg[i]= "0/0"
				elif ((ref != csp[0]) and (ref == csp[1]) and ((csp[0] != ".") and (csp[1] != "."))):
					chg_zyg[i]= "1/0"
				elif ((ref == csp[0]) and (ref != csp[1]) and ((csp[0] != ".") and (csp[1] != "."))):
					chg_zyg[i]= "0/1"
				elif ((ref != csp[0]) and (ref != csp[1]) and ((csp[0] != ".") and (csp[1] != "."))):
					chg_zyg[i]= "1/1"
				elif ((csp[0] == ".") and (csp[1] == ".")):
					chg_zyg[i]= cs
				elif ((csp[1] == ".")):
					if (ref == csp[0]):
						chg_zyg[i]= "0/."
					elif (ref != csp[0]):
						chg_zyg[i]= "1/."
				elif ((csp[0] == ".")):
					if (ref == csp[1]):
						chg_zyg[i]= "./0"
					elif (ref != csp[1]):
						chg_zyg[i]= "./1"
				i+=1
			return list(chg_zyg.values())
		####--------------------------------------

		with open(WORK+"snp_clsters2_ws3_clstronly_annt_snv_clstronly.table",'r') as f1, open(WORK+"temp_file1",'w') as output:
			first_line = f1.readline().strip()
			zyg_head = '\t'.join(first_line.split()[13:])
			output.write(first_line+"\t"+zyg_head+"\t"+str("protein_coding_EFF	AAChange	Refcodon	Altcodon")+"\n")
			for line in f1:
				line1 = line.strip()
				line2 = line1.split("\t")
				line3 = line2[8].split("|")
				if ((line2[5] =="exonic") and (line2[6] == "nonsynonymous_SNV" or line2[6] == "synonymous_SNV") and (line3[1] == "MISSENSE" or line3[1] == "SILENT")):
					True
					linesp = line3[2].split("/")
					ref = line2[2]; alt = line2[3]	
					zyg = line2[13:]
					chgz = Change_zygo(ref, alt, zyg)
					chgz_out = '\t'.join(chgz)
					wrt = str(line3[1]+"\t"+line3[3]+"\t"+linesp[0]+"\t"+linesp[1])
					output.write(line1+"\t"+chgz_out+"\t"+wrt+"\n")
	
	###--------
	import string
	##-------------------------------
	gencode = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T','AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R','CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P','CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R','GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A','GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G','TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L','TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
	
	# a function to translate a single codon
	def translate_codon(self,codon):
		return self.gencode.get(codon.upper(), '#')
		# a function to split a sequence into codons
	def split_into_codons(self,dna,frame):
		codons = []
		for i in range(frame - 1, len(dna)-2, 3):
			codon = dna[i:i+3]
			codons.append(codon)
		return codons
	# a function to translate a dna sequence in a single frame
	def translate_dna_single(self,dna, frame=1):
		codons = self.split_into_codons(dna, frame)
		amino_acids = ''
		for codon in codons:
			amino_acids = amino_acids + self.translate_codon(codon)
		return amino_acids
	
	def TWO_VAR(self,workdir):
		"""
		###---two variants(2VAR) codan changes---
		##----------------------------------
		"""
		lines =  open(workdir+"temp_file1","r").read().splitlines()
		writ2 = self.write_file(workdir+"temp_file2")
		#---
		midline_head =  lines[0].strip().split("\t")
		try:
			midline_head.remove(midline_head[8])
			midline_headcrp = '\t'.join(midline_head)
		except ValueError:
			print "empty string"
		writ2.write(str(midline_headcrp+"\t"+"Altcodon_merge-2VAR"+"\t"+"AA-change-2VAR")+"\n")
		#---
		i=0; TRcode =""; protcode = ""
		for i in range(len(lines)):
			True
			#---
			midline_crp =  lines[i].strip().split("\t")
			try:
				midline_crp.remove(midline_crp[8])
				midline_crpit = '\t'.join(midline_crp)
			except ValueError:
				print "empty string"
			#---
			try:
				beforeline = lines[i-1].strip()
				line0 = beforeline.split("\t")[-3]
				beforeline1 =  re.findall("\d+", line0)
				
				midline = lines[i].strip()
				line1 = midline.split("\t")[-3]
				midline1 =  re.findall("\d+", line1)
				
				nextline = lines[i+1].strip()
				line2 = nextline.split("\t")[-3]
				nextline1 =  re.findall("\d+", line2)

				REFbf=[];
				if ((beforeline1[0] == midline1[0]) or (midline1[0] == nextline1[0])):
					True 
					spREFbf=[];lscod1=[];lscod2=[]
					REF= lines[i].strip().split("\t")[-2]
					if (midline1[0] == beforeline1[0]):
						REFbf = lines[i-1].strip().split("\t")[-2]

					line11 = lines[i].strip().split("\t")[-1]
					if (midline1[0] == nextline1[0]):
						line22 = lines[i+1].strip().split("\t")[-1]
					
					if REFbf!=[]:
						for cod in REFbf:
							spREFbf.append(cod)
					for cod in line11:
						lscod1.append(cod)
					try:
						True
						for cod in line22:
							lscod2.append(cod)
					except:
						False
					if (((lscod1[0].islower()==True and lscod2[0].islower()==True) or (lscod1[1].islower()==True and lscod2[1].islower()==True) or (lscod1[2].islower()==True and lscod2[2].islower()==True)) and ((lscod1[0] == lscod2[0]) or (lscod1[1] == lscod2[1]) or (lscod1[2] == lscod2[2]))):
						True
						threeltr_code =[]
						if ((lscod1[0].isupper()==True and lscod2[0].islower()==True) or (lscod1[0] == lscod2[0])):
							True
							threeltr_code.append(lscod1[0])
						elif((lscod1[0].islower()==True and lscod2[0].isupper()==True) or (lscod1[0] == lscod2[0])):
							True
							threeltr_code.append(lscod2[0])
	
						if((lscod1[1].isupper()==True and lscod2[1].islower()==True) or (lscod1[1] == lscod2[1])):
							True
							threeltr_code.append(lscod1[1])
						elif((lscod1[1].islower()==True and lscod2[1].isupper()==True) or (lscod1[1] == lscod2[1])):
							True
							threeltr_code.append(lscod2[1])

						if((lscod1[2].isupper()==True and lscod2[2].islower()==True) or (lscod1[2] == lscod2[2])):
							True
							threeltr_code.append(lscod1[2])
						elif((lscod1[2].islower()==True and lscod2[2].isupper()==True) or (lscod1[2] == lscod2[2])):
							True
							threeltr_code.append(lscod2[2])
						#-----------------------
						if(len(threeltr_code)==3):
							True
							TRcode = ''.join(threeltr_code)
						else:
							True
							TRcode = 'Multiallelic-t1'
						#-----------------------
					else:
						True
						TRcode = '...'
					
					if(TRcode != line22):
						True
						protcode = self.translate_dna_single(TRcode)
					else:				
						True
						TRcode = "."; #protcode = "...";

					if (REF == REFbf):
						True
						protcode = "#####";TRcode = "Multiallelic-t2"	
					#print midline1, line11,line22,TRcode,protcode
					writ2.write(str(midline_crpit+"\t"+TRcode+"\t"+protcode)+"\n")
			except IndexError:
				False
				fstlst1 = midline_crpit.split("\t")[0]
				if (fstlst1 != "CHROM"):
					True
					writ2.write(str(midline_crpit+"\t"+"."+"\t"+".")+"\n")
		print "Done 2VAR of protein coding variants..."
	
	def THREE_VAR(self,workdir):
		"""
		###---three variants(3VAR) codan changes---
		##----------------------------------
		"""
		lines =  open(workdir+"temp_file2","r").read().splitlines()
		writ3 = self.write_file(workdir+"temp_file3")
		#---
		midline_head =  lines[0].strip()
		writ3.write(str(midline_head+"\t"+"Altcodon_merge-3VAR"+"\t"+"AA-change-3VAR")+"\n")
		#---
		i=0; TRcode =""; protcode = ""
		for i in range(len(lines)):
			True
			try:
				#---
				midline_crp =  lines[i].strip()
				#---
				beforeline = lines[i-1].strip()
				line0 = beforeline.split("\t")[-5]
				beforeline1 =  re.findall("\d+", line0)
				
				midline = lines[i].strip()
				line1 = midline.split("\t")[-5]
				midline1 =  re.findall("\d+", line1)
				
				nextline = lines[i+1].strip()
				line2 = nextline.split("\t")[-5]
				nextline1 =  re.findall("\d+", line2)

				line11=[];line22=[];writelst= [];
				if ((beforeline1[0] == midline1[0]) or (midline1[0] == nextline1[0])):
					True
					lscod1=[];lscod2=[]
					line11 = lines[i].strip().split("\t")[-2]
					if (midline1[0] == nextline1[0]):
						line22 = lines[i+1].strip().split("\t")[-2]
					
					for cod in line11:
						lscod1.append(cod)
					for cod in line22:
						lscod2.append(cod)
					#-----------------------
					if(len(lscod1)==3 and len(lscod2)==3):
						True
						threeltr_code =[];

						if ((lscod1[0].isupper()==True) and (lscod2[0].islower()==True) or (lscod1[0].isupper()==True and lscod2[0].isupper()==True) ):
							True
							threeltr_code.append(lscod1[0])
						elif ((lscod1[0].islower()==True) and (lscod2[0].isupper()==True)):
							True
							threeltr_code.append(lscod2[0])

						if ((lscod1[1].isupper()==True) and (lscod2[1].islower()==True) or (lscod1[1].isupper()==True and lscod2[1].isupper()==True)):
							True
							threeltr_code.append(lscod1[1])
						elif ((lscod1[1].islower()==True) and (lscod2[1].isupper()==True)):
							True
							threeltr_code.append(lscod2[1])
			
						if ((lscod1[2].isupper()==True) and (lscod2[2].islower()==True) or (lscod1[2].isupper()==True and lscod2[2].isupper()==True)):
							True
							threeltr_code.append(lscod1[2])
						elif ((lscod1[2].islower()==True) and (lscod2[2].isupper()==True)):
							True
							threeltr_code.append(lscod2[2])
						#-----------------------
						if(len(threeltr_code)==3):
							True
							TRcode = ''.join(threeltr_code)
						else:
							True
							TRcode = 'Multiallelic-t1'
						#-----------------------
					else:
						True
						TRcode = '.'
					
					if(TRcode != line22):
						True
						protcode = self.translate_dna_single(TRcode)
						if len(protcode) ==0:
							protcode = ".";TRcode = "."
					else:				
						True
						TRcode = "."; #protcode = "."

					if protcode == "" or TRcode == "":
						protcode = ".";TRcode = "."
					#-----------------------
					#print midline1,line11,line22,TRcode,protcode
					writ3.write(str(midline_crp+"\t"+TRcode+"\t"+protcode)+"\n")
					
					##----------------------
			except IndexError:
				False
		print "Done 3VAR of protein coding variants..."
		return None
		
	def PARS_OUT_VAR(self,readfl3,workdir):
		"""
		###---Pars variants(2VAR _ 3VAR) based on change of protein codons ---
		##----------------------------------
		"""
		lines =  open(workdir+"temp_file3","r").read().splitlines()
		writ4 = self.write_file(workdir+readfl3)

		#---
		midline_head =  lines[0].strip()
		midline_head1 = str('\t'.join(midline_head.split()[:8])+"\t"+'\t'.join(midline_head.split()[-8:])+"\t"+'\t'.join(midline_head.split()[8:-8]))
		writ4.write(str(midline_head1)+"\n")
		#---
		i=0; TRcode ="";
		for i in range(len(lines)):
			True
			try:
				#---
				midline_crp =  lines[i].strip()
				midline_crp1 = str('\t'.join(midline_crp.split()[:8])+"\t"+'\t'.join(midline_crp.split()[-8:])+"\t"+'\t'.join(midline_crp.split()[8:-8]))
				#---
				beforeline = lines[i-1].strip()
				line0 = beforeline.split("\t")[-7]
				beforeline1 =  re.findall("\d+", line0)
				
				midline = lines[i].strip()
				line1 = midline.split("\t")[-7]
				midline1 =  re.findall("\d+", line1)
				
				nextline = lines[i+1].strip()
				line2 = nextline.split("\t")[-7]
				nextline1 =  re.findall("\d+", line2)

				if ((beforeline1[0] == midline1[0]) or (midline1[0] == nextline1[0])):
					True 				
					lscod1=[];lscod2=[];line11=[];line22=[]; protcode = "0"; 
					AAchng = lines[i].strip().split("\t")
					line11 = lines[i].strip().split("\t")[-4]
					AAchng_NXT2 = nextline.split("\t")[-3]
					AAchng_NXT3 = nextline.split("\t")[-1] 
					AAchng13 = AAchng[-3]
					AAchng15 = AAchng[-1]
					if (line11 != "Multiallelic-t1" and line11 != "Multiallelic-t2" and (AAchng13 != "#####") ):
						True
						###################
						alpha2 = re.findall("[A-Z]", line1)
						alpha3 = re.findall("[A-Z]", line2)
						#print line1,line2,alpha2,alpha3,AAchng13,AAchng15,line11,midline1[0], nextline1[0],AAchng_NXT2, AAchng_NXT3
						
						if (((len(alpha2) == 2) and (len(alpha3) == 2))):
							True
							if (((alpha2[1] != AAchng13) and (alpha3[1] != AAchng13)) and ((alpha2[1] != AAchng15) and (alpha3[1] != AAchng15)) ):
								True
								writ4.write(str(midline_crp1)+"\n")
						elif ((len(alpha2) == 2) and (len(alpha3) == 1)):
							True
							if ((alpha2[1] != AAchng13) and (alpha2[1] != AAchng15) and (alpha2[1] != AAchng_NXT2) and (alpha2[1] != AAchng_NXT3) ):
								True
								writ4.write(str(midline_crp1)+"\n")

						elif ((len(alpha2) == 1) and (len(alpha3) == 2)):
							True
							if ( (alpha3[1] != AAchng13) and (alpha3[1] != AAchng15) and (alpha3[1] != AAchng_NXT2) and (alpha3[1] != AAchng_NXT3)):
								True
								writ4.write(str(midline_crp1)+"\n")

 				##----------------------
				####################
			except IndexError:
				False
				#---
				midline_crp =  lines[i].strip()
				midline_crp1 = str('\t'.join(midline_crp.split()[:8])+"\t"+'\t'.join(midline_crp.split()[-8:])+"\t"+'\t'.join(midline_crp.split()[8:-8]))
				#---
				try:
					midline = lines[i].strip()
					line1 = midline.split("\t")[-7]
					midline1 =  re.findall("\d+", line1)
					
					nextline = lines[i+1].strip()
					line2 = nextline.split("\t")[-7]
					nextline1 =  re.findall("\d+", line2)
					
					alpha2 = re.findall("[A-Z]", line1)
					AAchng = lines[i].strip().split("\t")
					AAchng13 = AAchng[-3]
					if (midline_crp.split("\t")[0] != "CHROM"):
						True
						if (midline1[0] == nextline1[0] and (alpha2[1] != AAchng13)):
							True
							writ4.write(str(midline_crp1)+"\n")
				except:
					True
					writ4.write(str(midline_crp1)+"\n")
		
		subprocess.check_output("rm "+workdir+"temp_file1", shell=True)
		subprocess.check_output("rm "+workdir+"temp_file2", shell=True)
		subprocess.check_output("rm "+workdir+"temp_file3", shell=True)
		subprocess.check_output("rm "+workdir+"snp_clsters2_ws3*", shell=True)
		print "Done Pars 2VAR and 3VAR of protein coding variants..."
		return None
		
		###------
		print "Done all pairs of coding variants..."
		return None
		
# file1: Input files
# write1: Output file

if __name__ == "__main__":
	parser = ArgumentParser(
	description="-Script to check for pair of protein codon changes. ")

	parser.add_argument("-i", "--infile", dest="INPUTFile",default=False, required=True, help="Full name of .vcf annovar input file" )
	parser.add_argument("-o", "--outfile", dest="OUTPUTFile",default=False, required=True, help="Full output file name" )
	parser.add_argument("-d", "--dir", dest="Directory", default="False", required=True, help="Working directory; input and output file location full Path" )

	if len(sys.argv) == 2:
		parser.print_help()
		print "Example: ./Annotation_Protein-codon-changes.py -i HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_variants_GATK_step2_annovar.hg19_multianno.txt -o HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_variants_GATK_Protein-CODON_changes.txt -d /WORKDIR/FULLPATH/"
		sys.exit()

	args = parser.parse_args()

	INF = args.INPUTFile
	OUTF = args.OUTPUTFile
	DIR = args.Directory


clF1 = SearchDB().Search_CODON(INF,DIR)
clF2 = SearchDB().TWO_VAR(DIR)
clF3 = SearchDB().THREE_VAR(DIR)
clF4 = SearchDB().PARS_OUT_VAR(OUTF,DIR)


