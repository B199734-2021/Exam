python3
#!/usr/local/bin/python3
import os, subprocess, sys, re, shutil
from Bio.Blast import NCBIWWW
from Bio import SeqIO
import numpy as np
subprocess.call("pip3 install pandas", shell = True)
import pandas as pd
### create a work place for user
print("Welcome to use this pipeline. This pipeline is for Blast analysis. It can offer blastn, blastp, tblstn, tblastx. It is a interaction design. First We need you input  a work place.")
print("It is a interaction design. We need you input something. You should make a work place and specify protein family and taxonomic at begening. And according to some output to determine whether continue or exit in the progress. Let'begin!") 
work_place=str(input("Please make a dir for work:"))
os.mkdir(work_place) 
os.chdir(work_place)
###################################################################################################
# step 1: ask user whether make a own blast database
###################################################################################################
own_database_or_not=input("Would you like make a own database?/nyes or no?/n")
if own_database_or_not startswith("y"):
    def databasemake(protein_or_DNA="",Organism,protein_family="",gene_name=""):
         print("You have provide the database type: "+ protein_or_DNA+ ", Organism: "+Organism+", protein_family: "+protein_family+", gene_name: "+gene_name+". Please wait...")
         subprocess.call("esearch -db {0} -query {1} +"[ORGN] AND "+ {2}+{3}+"NOT PARTIAL" | efetch -format fasta >> download.fasta".format(protein_or_DNA,Organism,protein_family,gene_name), shell=True)
         if protein_or_DNA.lower startswith ("p"):
             dbtype=prot
         if protein_or_DNA.lower startswith ("d"):
             dbtype=nucl
         subprocess.call("makeblastdb -in download.fasta -dbtype {0} -out Target_database".format(dbtype),shell=True)
    details={}
    details["protein_or_DNA"] = input ("Please specify the type of database you want. protein or DNA?")
    details["Organism"]=input("Please specify the Organism.")
    details["prorein_family"]=input ("please specify the name of protein family. If you choose DNA you can press return.")
    details["gene_name"]=input ("please specify the gene name. If you choose protein you can press return.")
    databasemake(*list(details.values()))
    print("Database making, please wait...")

#####################################################################
# step 2: get input from user to blast
#####################################################################
query_type=input("please choose the type of your submission for blast./n The raw sequence press 1./nThe accession number press 2./nThe one query fasta file press 3./n The more than one query fasta file press 4.")   
query=input("please input the sequence or fasta file or accession number you want to blast.")
seq_type=input("please tell us the data type you submit, protein or DNA?")
#### get all submission data type as sequence type
if int(query_type) == 1:
    record=query
if int(query_type) == 2:
    #### get the sequence 
    subprocess.call("esearch {0} | efetch -format fasta >> query.fasta".format(query), shell=True)
    record=SeqIO.read(“query.fasta”,format= “fasta”)
if int(query_type) == 3:
    record=SeqIO.read(“query.fasta”,format= “fasta”)
if int(query_type) == 4:
    record=next(SeqIO.parse(“query”,“fasta”))
#### check the choice by user by calculate the ATCG contents to ensure the query type
length=len(str(record))
a_count = record.upper().count('A')
t_count = record.upper().count('T')
c_count = record.upper().count('C')
g_count = record.upper().count('G')
ATCGcontenes=((a_count+t_count+c_count+g_count)/length)*100
if ATCGcontenes >= 80:
    seq_type=DNA
else:
    seq_type=protein
#####################################################
# step 3 blast
#####################################################
databasetype=input("pleass input the data base type you want. DNA or protein?")
if databasetype.lower startswith ("p"):
    db=protein
if databasetype.lower startswith ("d"):
    db=DNA
### if user specify own database
if own_database_or_not.lower startswith("y"):
    if dbtype=prot:
        db=protein
    if dbtype=nucl:
        db=DNA
if db=="protein" && seq_type=="protein":
    blastchoice="blastp"
if db=="protein" && seq_type=="DNA":
    blastchoice="blastx"
if db=="DNA" && seq_type=="protein":
    blastchoice="tblastn"
if db=="DNA" && seq_type=="DNA":  
    blastchoice="blastn"
result_handle=NCBIWWW.qblast(“blastchoice”, “nt”, record.format(“seq”))
if own_database_or_not.lower startswith("y"):
    subprocess.call("{0} -db {1} -query {2} -outfmt 7 -out blastoutput.txt".format(blastchoice,Target_database,record),shell=True)
else:
    subprocess.call("{0} -db nr -query {1} -outfmt 7 -out blastoutput.txt".format(blastchoice,record),shell=True)

#### check the choice by user by calculate the ATCG contents to ensure the query type
    length=len(str(query))
    a_count = query.upper().count('A')
    t_count = query.upper().count('T')
    c_count = query.upper().count('C')
    g_count = query.upper().count('G')
    ATCGcontenes=((a_count+t_count+c_count+g_count)/length)*100
    if ATCGcontenes >= 80:
        seq_type=DNA
    else:
        seq_type=protein

### use esearch and efetch get the sequences, accession number, species
#####################################################
# step 4 analysis for blast result
#####################################################
df = pd.read_csv('blastoutput.txt', sep="\t", na_values=['-'])
print("Here is the top10 blast bitscore information")
df.sort_valuse("bitscroe",ascending=False).head(10)
df[df["mismatches"]>=20]
print("Here is the sequences more than 20 mismatches.")
df[df[alignment_length]>100]
print("Here is the sequences have more than 100 alignment length.")
print("All the work is finished! \nBye~")
