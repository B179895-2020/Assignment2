
#!/bin/usr/python3

import sys, subprocess, shutil, os
import string, re
import collections
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from collections import OrderedDict

print("\n###WHEN INPUTING THE TAXONOMIC GROUP AND THE PROTEIN FAMILY, PLEASE CHECK THEIR NAMES ARE OFFICIAL!!!\n")

protein= input("\nPlease enter the PROTEIN you are searching for\n")

def count1(name):
	num = int(subprocess.check_output("esearch -db protein -query '{0}' | xtract -pattern ENTREZ_DIRECT -element Count".format(name),shell=True))
	if num == 0:
		return 0
	return 1

while count1(protein) == 0:
	print("ERROR!\n***CANNOT FIND "+protein+" IN THE DATABASE, PLEASE CHECK IT AGAIN!***\n")
	protein=input("\nPlease enter the right name of the protein family:\n")

taxonomic= input("\nPlease enter the TAXONOMIC GROUP you are searching for\n")

while count1(taxonomic) == 0:
	print("ERROR!\n***CANNOT FIND THE "+taxonomic+" IN THE DATABASE, PLEASE CHECK IT AGAIN!***\n")
	taxonomic = input("\nPlease enter the right name of the taxonomic group:\n")

def count2(pro,tax):
	num=int(subprocess.check_output("esearch -db protein -query '{0}[prot]' -organism '{1}' | xtract -pattern ENTREZ_DIRECT -element Count".format(pro,tax),shell=True))
	if num == 0:
		return 0
	return 1

while count2(protein, taxonomic) == 0:
	print("ERROR!\n***CANNOT FIND \" {0} \" WITH \" {1} \" IN THE DATABASE\n".format(protein,taxonomic))
	protein = input("\nPlease enter a valid protein name again\n")
	while count1(protein) == 0:
		print("ERROR!\n***CANNOT FIND "+protein+" IN THE DATABASE, PLEASE CHECK IT AGAIN!***\n")
		protein=input("\nPlease enter the right name of the protein family:\n")
	taxonomic = input("\nPlease enter a valid taxonomic name again\n")
	while count1(taxonomic) == 0:
		print("ERROR!\n***CANNOT FIND THE "+taxonomic+" IN THE DATABASE, PLEASE CHECK IT AGAIN!***\n")
		taxonomic = input("\nPlease enter the right name of the taxonomic group:\n")

print("\n###Thanks, you have chosen " + protein + " in " + taxonomic + " and a new directory will be made for you!\n")

nowdir = os.getcwd()
os.mkdir('{0}/{1}_{2}'.format(nowdir,protein,taxonomic))
os.chdir('{0}/{1}_{2}'.format(nowdir,protein,taxonomic))

subprocess.call('esearch -db protein -query "{0}[prot]" -organism "{1}" | efetch -db protein -format docsum | xtract -pattern DocumentSummary -element Organism > {0}_{1}.txt'.format(protein,taxonomic), shell=True)

sequences=len(open('{0}_{1}.txt'.format(protein,taxonomic)).readlines())
species=len(set(open('{0}_{1}.txt'.format(protein,taxonomic)).readlines()))

if sequences >= 1000:
	print("\n###Sorry! Your starting sequence is over 1000. The program are going to exit. Please try other groups later\n")
	exit()

print(str("\n###Congratulations!\n###After searching, there is "+str(sequences)+" sequences in "+str(species)+" species.\n"))

reply1 = input("\nIf you choose to continue with the current dataset, please enter YES! Otherwise, you will exit.\n")

if not (re.search('[YES]',reply1.upper())):
	print("\n###Expect for your next try!\n")
	exit()

subprocess.call('esearch -db protein -query "{0}[prot]" -organism "{1}" | efetch -db protein -format fasta > {0}_{1}.fasta'.format(protein,taxonomic), shell=True)

subprocess.call('clustalo -i {0}_{1}.fasta -o {0}_{1}_align.fasta -v'.format(protein,taxonomic), shell =True)

reply2 = input("\nIf you choose to see the output of CLUSTALO, please enter YES! Otherwise, continue the next step.\n")

if (re.search('[YES]',reply2.upper())):
	subprocess.call('prettyplot -sequence {0}_{1}_align.fasta -graph svg'.format(protein,taxonomic), shell=True)
	subprocess.call('firefox prettyplot.svg&', shell=True)

print("\n###prettyplot.svg is saved as the output picture of clustalo!\n")

reply3=input("\nIf you want to scan the changes in detail using INFOALIGN, please enter YES! Otherwise, continue the next steps.\n")

if (re.search('[YES]',reply3.upper())):
	subprocess.call("infoalign -sequence {0}_{1}_align.fasta -only -name -outfile output1.infoalign".format(protein,taxonomic), shell=True)
	subprocess.call("infoalign -sequence {0}_{1}_align.fasta -only -change -outfile output2.infoalign".format(protein,taxonomic), shell=True)
	file1=open("output1.infoalign").read().rstrip("\n") 
	keys =list(file1.split("\n")) 
	file2=open("output2.infoalign").read().rstrip("\n") 
	list2=file2.split("\n") 
	values = list(map(float, list2))
	dic = dict(zip(keys, values))
	dic = sorted(dic.items(), key=lambda values:values[1])
	print("\n###This is %change of each sequence from low to hign.\n")
	print(dic)
	print("\n")

subprocess.call('cons -plurality 0.8 -sequence {0}_{1}_align.fasta -outseq {0}_{1}_consensus.fasta'.format(protein, taxonomic), shell=True) 

subprocess.call('makeblastdb -in {0}_{1}.fasta -dbtype prot -out {0}_{1}_database'.format(protein,taxonomic), shell=True)

subprocess.call('blastp -db {0}_{1}_database -query {0}_{1}_consensus.fasta -outfmt 7 -out {0}_{1}_blast.txt'.format(protein,taxonomic), shell = True)

lines = (i for i in open('{0}_{1}_blast.txt'.format(protein,taxonomic), 'r') if '#' not in i )
file = open('{0}_{1}_blast_lines.txt'.format(protein,taxonomic), 'w', encoding="utf-8")
file.writelines(lines)
file.close()

file = open("{0}_{1}_blast_lines.txt".format(protein,taxonomic)).read().rstrip('\n')
dic={}
lines = file.split('\n')
for i in lines:
	part=i.split('\t')
	acc=part[1]
	hsp=part[-1]
	dic[acc]=hsp

dic = sorted(dic.items(),key=lambda item:item[1],reverse=True)[:250]
dic250 = {}
for i in dic:
	dic250[i[0]] = i[1]

acc = open("acc.txt", "w")
acc.write("\n".join(dic250.keys())) 
acc.close()

subprocess.call('/localdisk/data/BPSM/Assignment2/pullseq -i {0}_{1}_align.fasta -n acc.txt > 250_similar_seq.fasta -v'.format(protein,taxonomic), shell=True)

reply4 = input("\nIf you choose to see the output of BLAST, please enter YES! Otherwise, continue the next step.\n")

if (re.search('[YES]',reply4.upper())):
	subprocess.call('plotcon -sequence 250_similar_seq.fasta -winsize 4 -graph svg', shell=True)
	subprocess.call('eog plotcon.svg&', shell = True)

print("\n###plotcon.svg is saved as the output picture!\n")







print('\n###NOW YOU HAVE FINISHED ALL THE STEPS!\n###FILES ARE SAVED IN THE DIRECTORY MADE FOR YOU, YOU CAN CHECK THEM WHENVER YOU WANT.\n\n###THANKS FOR USEING!###\n')
