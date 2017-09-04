#!/usr/bin/python


"""
#Script for annotation...

#####------Inputs-------
# Annotation_withlocal-CADD-ReMM-LINSIGHT_sampall.py INPUT=Inputfile OUTPUT=Outputfile DIR=Directory-fullpath

#./Annotation_withlocal-CADD-ReMM-LINSIGHT_sampall.py -i HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_Variants_Zyg_PASS_MAF5%.txt -o HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_Variants_Zyg_PASS_MAF5%_CRLD.txt -d ./Example_WORKDIR/
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

	def Search_CADD(self,readfl1,writedir):
		"""
		Calling Search CADD localsearch method
		"""
		def srchdb1(self,RR):
			try:
				True
				cmdFls1 = "tabix -f "+writedir+readfl1+".tsv.gz "+RR[0]+":"+RR[8]+"-"+RR[8]
				grepout =  subprocess.check_output(cmdFls1, shell=True)
			except subprocess.CalledProcessError:
				False
				grepout = "."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\n"
			return grepout

		cmd1 = "gzip -c "+writedir+readfl1+" > "+writedir+readfl1+".gz"
		cmd2 = "/media/varma/SAMSUNG/CADD/v1.3/CADD_v1.3/bin/score.sh "+writedir+readfl1+".gz "+writedir+readfl1+".tsv.gz"
		cmd3 = "tabix -s 1 -b 2 -e 2 "+writedir+readfl1+".tsv.gz"
		cmd4 = "rm "+writedir+readfl1+".gz"
		subprocess.check_output(cmd1, shell=True)
		subprocess.check_output(cmd2, shell=True)
		subprocess.check_output(cmd3, shell=True)
		subprocess.check_output(cmd4, shell=True)

		with open(writedir+readfl1,'r') as f1, open(writedir+"temp-outfile1",'w') as output:
			first_line = f1.readline().strip()
			HeadDB1 = "CADD1.3_RawScore"+"\t"+"CADD1.3_PHRED"
			output.write(first_line+"\t"+HeadDB1+"\n")
			for g in f1:
				g1 = g.strip()
				gg = g1.split("\t")
				CADDout = srchdb1(self,gg)
				CDS = CADDout.split()
				if len(CDS)==0:
					output.write(str(g1+"\t"+"."+"\t"+"."+"\n"))
				else:
					True
					output.write(str(g1+"\t"+CDS[4]+"\t"+CDS[5]+"\n"))
		subprocess.check_output("rm "+writedir+readfl1+".tsv.gz*", shell=True)
		print "Done seach for CADD...."
		return None

	def Search_ReMM(self,writedir):
		"""
		Calling Search REMM localsearch method
		"""
		def srchdb2(RR):
			try:
				True
				#tabix -f remmData.tsv.gz 6:88757669-88757669
				cmdFls1 = "tabix -f /media/varma/SAMSUNG/ReMM_score/ReMM.v0.3.1.tsv.gz "+RR[0]+":"+RR[8]+"-"+RR[8]
				grepout =  subprocess.check_output(cmdFls1, shell=True)
			except:
				False
				grepout = "."+"\t"+"."+"\t"+"."+"\n"
			return grepout

		with open(writedir+"temp-outfile1",'r') as f1, open(writedir+"temp-outfile2",'w') as output:
			first_line = f1.readline().strip()
			HeadDB1 = "ReMM_Scores"
			output.write(first_line+"\t"+HeadDB1+"\n")
			for g in f1:
				g1 = g.strip()
				gg = g1.split("\t")
				REMMout = srchdb2(gg)
				CDS = REMMout.split()
				if len(CDS)==0:
					output.write(str(g1+"\t"+"."+"\n"))
				else:
					output.write(str(g1+"\t"+CDS[2]+"\n"))

		
		print "Done seach for REMM...."
		return None

	def Search_LINSIGHT(self,writedir):
		"""
		Calling Search LINSIGHT localsearch method
		"""
		def srchdb3(RR):
			try:
				True
				cmdFls1 = "tabix -f /media/varma/SAMSUNG/LINSIGHT/LINSIGHT.bed.gz chr"+RR[0]+":"+RR[8]+"-"+RR[8]
				grepout =  subprocess.check_output(cmdFls1, shell=True)
			except:
				False
				grepout = "."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\n"
			return grepout

		with open(writedir+"temp-outfile2",'r') as f1, open(writedir+"temp-outfile3",'w') as output:
			first_line = f1.readline().strip()
			HeadDB1 = " LINSIGHT_Scores"
			output.write(first_line+"\t"+HeadDB1+"\n")
			for g in f1:
				g1 = g.strip()
				gg = g1.split("\t")
				LINSIGHTout = srchdb3(gg)
				CDS = LINSIGHTout.split()
				if len(CDS)==0:
					output.write(str(g1+"\t"+"."+"\n"))
				else:
					output.write(str(g1+"\t"+CDS[4]+"\n"))
		print "Done seach for LINSIGHT...."
		return None

	def Search_DANN(self,writefl1,writedir):
		"""
		Calling Search DANN localsearch method
		"""
		def srchdb3(RR):
			try:
				True
				cmdFls1 = "tabix -f /media/varma/SAMSUNG/DANN/DANN_whole_genome_SNVs.tsv.bgz "+RR[0]+":"+RR[8]+"-"+RR[8]
				grepout =  subprocess.check_output(cmdFls1, shell=True)
			except:
				False
				grepout = "."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\n"
			return grepout

		with open(writedir+"temp-outfile3",'r') as f1, open(writedir+writefl1,'w') as output:
			first_line = f1.readline().strip()
			HeadDB1 = "DANN_Scores(Deep_learning_approach)"
			output.write(first_line+"\t"+HeadDB1+"\n")
			for g in f1:
				g1 = g.strip()
				gg = g1.split("\t")
				SCOREout = srchdb3(gg)
				CDS = SCOREout.split()
				if len(CDS)==0:
					output.write(str(g1+"\t"+"."+"\n"))
				else:
					output.write(str(g1+"\t"+CDS[4]+"\n"))
		subprocess.check_output("rm "+writedir+"temp-outfile1", shell=True)
		subprocess.check_output("rm "+writedir+"temp-outfile2", shell=True)
		subprocess.check_output("rm "+writedir+"temp-outfile3", shell=True)
		print "Done seach for DANN...."
		return None

# file1: Input positions SNP
# write1: Output file


if __name__ == "__main__":
	parser = ArgumentParser(
	description="Script to annotate with different scoring systems such as: 1) CADD-scores 2) ReMM-scores 3) LINSIGHT-scores and 4) DANN-scores.")

	parser.add_argument("-i", "--infile", dest="INPUTFile",default=False, required=True, help="Full input file name" )
	parser.add_argument("-o", "--outfile", dest="OUTPUTFile",default=False, required=True, help="Full output file name" )
	parser.add_argument("-d", "--dir", dest="Directory", default="False", required=True, help="Working directory; input and output file location full Path" )

	if len(sys.argv) == 2:
		parser.print_help()
		print "Example: ./Annotation_withlocal-CADD-ReMM-LINSIGHT_sampall.py -i HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_Variants_Zyg_PASS_MAF5%.txt -o HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_Variants_Zyg_PASS_MAF5%_CRLD.txt -d /WORKDIR/FULLPATH/"
		sys.exit()

	args = parser.parse_args()

	INF = args.INPUTFile
	OUTF = args.OUTPUTFile
	DIR = args.Directory

	clF1 = SearchDB().Search_CADD(INF,DIR)
	clF2 = SearchDB().Search_ReMM(DIR)
	clF3 = SearchDB().Search_LINSIGHT(DIR)
	clF3 = SearchDB().Search_DANN(OUTF,DIR)

