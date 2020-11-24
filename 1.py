
#!/bin/usr/python3

import subprocess
import os
import re
import shutil

print("\n###THIS PROGRAM WILL HELP YOU TO FIND A DATABASE AND DO SOME COMPARATION.####\n")
print("\n\n###WHEN INPUTING THE TAXONOMIC GROUP AND THE PROTEIN FAMILY, PLEASE CHECK THEIR NAMES ARE OFFICIAL!!!\n")
print("\n###YOU CAN ONLY PUT LETTERS AND FIGURES IN THE INPUT LINES, OR YOU WILL GET WRONG. THANKS!!!###\n")

#STEP1:get input of the protein family and the taxonomic group, and do some checks.[TASK1]


protein= input("\nPlease enter the PROTEIN you are searching for\n")

#use ENTREZ_DIRECT to calculate the count. If there is no name in the database, let the user try again later.
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

#check if there is protein inputted with organism inputted in the database.
def count2(pro,tax):
	num=int(subprocess.check_output("esearch -db protein -query '{0}[prot]' -organism '{1}' | xtract -pattern ENTREZ_DIRECT -element Count".format(pro,tax),shell=True))
	if num == 0:
		return 0
	return 1

#use while loops to repeat, until get right inputs.
#check the name seperately, and both names again and again.
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

#make a new directory to store all the data.
print("\n###Thanks, you have chosen " + protein + " in " + taxonomic + " and a new directory will be made for you!\n")

nowdir = os.getcwd()
os.mkdir('{0}/{1}_{2}'.format(nowdir,protein,taxonomic))
os.chdir('{0}/{1}_{2}'.format(nowdir,protein,taxonomic))

#to get the counts of species and sequences, downloaded the data.
subprocess.call('esearch -db protein -query "{0}[prot]" -organism "{1}" | efetch -db protein -format docsum | xtract -pattern DocumentSummary -element Organism > {0}_{1}.txt'.format(protein,taxonomic), shell=True)

#there maybe some repeated sequnces in species, so use set.
sequences=len(open('{0}_{1}.txt'.format(protein,taxonomic)).readlines())
species=len(set(open('{0}_{1}.txt'.format(protein,taxonomic)).readlines()))

#allowable starting sequence set shouldn't have more than 1000.
if sequences >= 1000:
	print("\n###Sorry! Your starting sequence is over 1000. The program are going to exit. Please try other groups later\n")
	exit()

#ask the user if he/she want to continue with such dataset.
print(str("\n###Congratulations!\n###After searching, there is "+str(sequences)+" sequences in "+str(species)+" species.\n"))

reply1 = input("\nIf you choose to continue with the current dataset, please enter YES! Otherwise, you will exit.\n")

if (re.search('[YES]',reply1.upper())):
	print("\n###Please wait for the clustalo process, which may take several seconds...\n")
else:
	print("\n###Expect for your next try!\n")
	shutil.rmtree('{0}/{1}_{2}'.format(nowdir,protein,taxonomic))
	exit()

#STEP2:determine, and plot, the level of conservation between the protein sequences.(TASK2)


subprocess.call('esearch -db protein -query "{0}[prot]" -organism "{1}" | efetch -db protein -format fasta > {0}_{1}.fasta'.format(protein,taxonomic), shell=True)

#align all the protein sequences using clustalo
subprocess.call('clustalo -i {0}_{1}.fasta -o {0}_{1}_align.fasta -v'.format(protein,taxonomic), shell =True)

#ask the user if he/she want to visualise the result.
reply2 = input("\nIf you choose to see the output of CLUSTALO, please enter YES! Otherwise, continue the next step.\n")
subprocess.call('prettyplot -sequence {0}_{1}_align.fasta -graph svg'.format(protein,taxonomic), shell=True)

if (re.search('[YES]',reply2.upper())):
	subprocess.call('firefox prettyplot.svg&', shell=True)

print("\n###prettyplot.svg is saved as the output picture of clustalo!\n")

#another emboss analysis to show some biological information about the output.(TASK4)


#ask the user if he/she want to see more details about clustalo.
reply3=input("\nIf you want to scan the percent of changes in detail using INFOALIGN, please enter YES! Otherwise, continue the next steps.\n")

if (re.search('[YES]',reply3.upper())):
	#download name part and %change part, and make a dictionary after sort.
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

#use cons to calculate a consenus sequence from the multiple sequence alignment for futher blast.
subprocess.call('cons -plurality 0.8 -sequence {0}_{1}_align.fasta -outseq {0}_{1}_consensus.fasta'.format(protein, taxonomic), shell=True) 

#make database and process protein blast.
subprocess.call('makeblastdb -in {0}_{1}.fasta -dbtype prot -out {0}_{1}_database'.format(protein,taxonomic), shell=True)
subprocess.call('blastp -db {0}_{1}_database -query {0}_{1}_consensus.fasta -outfmt 7 -out {0}_{1}_blast.txt'.format(protein,taxonomic), shell = True)

#remove these comment lines.
lines = (i for i in open('{0}_{1}_blast.txt'.format(protein,taxonomic), 'r') if '#' not in i )
file = open('{0}_{1}_blast_lines.txt'.format(protein,taxonomic), 'w', encoding="utf-8")
file.writelines(lines)
file.close()

#make directory for the subject acc and bit score.
file = open("{0}_{1}_blast_lines.txt".format(protein,taxonomic)).read().rstrip('\n')
dic={}
lines = file.split('\n')
for i in lines:
	part=i.split('\t')
	acc=part[1]
	hsp=part[-1]
	dic[acc]=hsp

#select 250 the most similar sequences by their bit score.
dic = sorted(dic.items(),key=lambda item:item[1],reverse=True)[:250]
dic250 = {}
for i in dic:
	dic250[i[0]] = i[1]

#collect the 250 subject accs.
acc = open("subject_acc.txt", "w")
acc.write("\n".join(dic250.keys())) 
acc.close()

#given a list(250) of head ids to extract sequences from the fasta file in to a new file.
subprocess.call('/localdisk/data/BPSM/Assignment2/pullseq -i {0}_{1}_align.fasta -n subject_acc.txt > 250_similar_seq.fasta -v'.format(protein,taxonomic), shell=True)

#ask the user if he/she want to visualise the result.
reply4 = input("\nIf you choose to see the output of BLAST, please enter YES! Otherwise, continue the next step.\n")
subprocess.call('plotcon -sequence 250_similar_seq.fasta -winsize 4 -graph svg', shell=True)

if (re.search('[YES]',reply4.upper())):
	subprocess.call('eog plotcon.svg&', shell = True)

print("\n###plotcon.svg is saved as the output picture of blastp!\n")

#STEP3:scan protein sequence(s) of interest with motifs from the PROSITE database.


#use seqretsplit to seperate each sequence and put them into different files.
#before that, make a new directory to store these files.
print("\n###When seqretsplit, you need to press enter key to make it continue...\n")
os.mkdir("motif")
os.chdir("motif")
subprocess.call("cp ../250_similar_seq.fasta .", shell=True) 
subprocess.call("cp ../subject_acc.txt .", shell=True)
subprocess.call("seqretsplit -sequence 250_similar_seq.fasta -sformat fasta -osformat fasta",shell=True)
file1=open("subject_acc.txt")
sub_name = file1.read().splitlines()

#use patmatmotifs to scan a protein sequence with motifs from the PROSITE database.
print("\n###Trying to seek motifs in the subset from the PROSITE database, please wait...\n")
for i in sub_name:
	i=i.lower()
	subprocess.call("patmatmotifs -sequence {0}.fasta -outfile {0}.patmatmotifs -full".format(i), shell=True)
	file2 = open("{0}.patmatmotifs".format(i)).read().splitlines()
	for j in file2:
		if re.search('#',j):
			next
		elif re.search('Motif', j):
			print("\n###In the "+i+", find "+j)

print("\n###Everytime a motif found in a sequence, their names both showed in the screen!\n")

print('\n###NOW YOU HAVE FINISHED THIS PROGRAM! ALL THE FILES ARE SAVED IN THE DIRECTORY MADE FOR YOU.\n\n###THANKS FOR USEING!###\n')

