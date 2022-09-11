#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 27 23:30:32 2021

@author: terynmorgan
"""

import numpy as np 
import pandas as pd 
import os
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import pairwise2
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
from scipy.stats import normaltest
from scipy.stats import kruskal
from scipy import mean

#PART 1: develop metadata for the COVID seq genomes
#extracts file names from folder 'Genome dataset'
genome_data=[] #holds the file names from directory
path='Genome dataset' #path of directory to genomes
fileList = os.listdir(path) #gets file list of directory
for line in fileList: #appends file names to genome_data list
    genome_data.append(line)
genome_data.sort() #sorts in alphabetical ascending order
del genome_data[0] #corrupted file had to be removed

#extracted sequences form file name
genomes=[] #holds all the genomic sequences 
for idx in enumerate(genome_data):#get a counter and the value from the iterable
    #open path for each txt genomic file {path} from genome_data using counter int from enumerate
    seq=open('Genome dataset/{path}'.format(path=genome_data[idx[0]]),'r')
    #appends each file content (genome) into list, descards unnecessary characters
    #appends string of text format genome
    genomes.append(seq.read().replace('\n', '').replace('\ufeff', '')) 
#genome[22]= 'MN997409.txt' needs following code to fix formatting:
new22=''
for i in genomes[22]:
    if i!=' ' and i.isnumeric()==False:
        new22+=i
genomes[22]=new22.upper()

#for each genome, for loop gets GC content, length, and molecular weight 
    #adds to cooresponding list
GC_content=[]
base_pairs=[]
mol_weight=[]
for idx, i in enumerate(genomes): #loops over list that holds all genomic sequences
    GC_content.append(GC(i))
    base_pairs.append(len(i))
    X = ProteinAnalysis(genomes[idx])
    mol_weight.append(X.molecular_weight())
GC_mean= mean(GC_content) #37.74504258968151
    
#creates name labels to be inserted into dataframe
names=[]
names2=[]
for i in range(1,43):
    names.append('Sequence {}'.format(i))
    names2.append(str(i))

#creates dataframe to hold the metadata
    #will contain Name (Sequence 1,...Sequence 43), base_pairs (len), GC_Content
df=pd.DataFrame({'Name': names, 'Base Pairs Length': base_pairs,'GC Content': GC_content, 'Molecular Weight': mol_weight})
df.to_csv('COVID metadata.csv', index=False)

#Graphing GC Content
plt.rcdefaults()
fig, ax = plt.subplots()
width = 0.5
ind = np.arange(len(names2))
plt.barh([i*1.5 for i in ind],GC_content,width,align='center')
ax.invert_yaxis()
ax.set_yticks(ind*1.5)
ax.set_yticklabels(names2)
ax.set_title('Sequence Analysis: GC_Content')
ax.set_ylabel('Sequences')
ax.set_xlabel('GC_content')
plt.show()

#PART 2: Pairwise Sequence Alignment: following function finds aligments score
#n=genome sequence
#returns list of simularity matches btwn inputted gene to every other gene in genome_data
def pairwise(n):
    perc_match=[]
    for seq in genomes:
        #finds the alignment score between inputted genome and given genome in genomes list 
        alignments=pairwise2.align.globalxx(n,seq,one_alignment_only=True, score_only=True)
        #converts this alignment score to a percent: shows Simularity % between 2 sequences
        perc_match.append(alignments/len(n) *100)
    return perc_match

    #***THIS FOR LOOP TOOK OVER 3 HOURS TO RUN SO BEWARE***
    #***GOING THORUGH 42X42 SEQUENCE FOR PAIRWISE ANALYSIS***
# ==================================================================
#pairwise sequence alignment can only compare two genes to each other at a time 
#for loop traverses genomes and runs pairwise() for each 
#puts each df into a list
L=[]
for i in genomes:
    L.append(pairwise(i))
#combines list of lists into df  
pairwise_df=pd.DataFrame(L)
#renames column names
pairwise_df.set_axis(names, axis=1, inplace=True)
pairwise_df.set_axis(names,inplace=True)
pairwise_df.to_csv('Pairwise Sequence Alignment2.csv')

#plots as heatmap
sns.heatmap(pairwise_df, cmap='RdBu')
plt.title('Pairwise Sequence Alignment Heatmap')
plt.savefig('Heatmap.png')
plt.show()

#plots as grouped bar graph
f = plt.figure(figsize=(12, 6), dpi=80)
pairwise_df.plot(kind='bar',ax=f.gca())
plt.title('Pairwise Sequence Alignment')
plt.xlabel('Sequences')
plt.ylabel('Pairwise Alignment Score')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()

#PART 3: Amino Acid Distribution
#takes a genome sequence, makes a Seq object, and returns the translated strand 
def gen_protein_seq(n):
    seq=Seq(n)
    return str(seq.translate)

#takes the translated strand, returns dict of amino acid freq.
def amino_acid_composition(protein_seq):
    #dic of amino acids
    amino_acids_dict = {'Y': 0, 'A': 0, 'Q': 0, 'D': 0, 'C': 0, 'G': 0, 'V': 0, 'T': 0, 'E': 0, 'N': 0, 
                       'K': 0, 'R': 0, 'S': 0, 'I': 0, 'H': 0, 'M': 0, 'F': 0, 'L': 0, 'W': 0, 'P': 0}
    for amino_acid in amino_acids_dict:
        #counts how many times the translated seq has the specified amino acid
        #divides by the len of the seq and converts to a normalized frequncy (*100)  
        amino_acids_dict[amino_acid] = protein_seq.count(amino_acid)/len(protein_seq)*100
    return amino_acids_dict
#traverses genomes and runs gen_protein_seq() to get translated seq
#then runs amino_acid_composition() for each and adds to a list
aa_dicts=[]
for i in genomes:
    aa_dicts.append(amino_acid_composition(gen_protein_seq(i)))

#makes new list of dicts that gets rid of amino acids w/ freq of 0
aa_dic=[]
for i in aa_dicts:
    aa_dic.append({x:y for x,y in i.items() if y!=0})
df3=pd.DataFrame(aa_dic)

#Amino Acid Distribution as grouped bar chart
f = plt.figure()
plt.title('Amino Acid Distribution', color='black')
#uses df3 to create plot
df3.T.plot(kind='bar', ax=f.gca())
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel('Amino Acids') 
plt.ylabel('Normalized Frequency(%)') 
plt.show()

#Part 4
#creates list into an array of all the T values from aa_dic
T_vals=[]
for i in aa_dic:
    T_vals.append(i.get('T'))
arr_T=np.array(T_vals)

#Box Plot to show if T Amino Acid contains an outlier
sns.boxplot(data=arr_T, orient="h", palette="Set2")
plt.title('T Amino Acid Boxplot')
plt.xlabel("T Normalized Frequency (%)")
plt.show()

#After removing outlier in AA of T
del T_vals[T_vals.index(max(T_vals))]
arr_T2=np.array(T_vals)
#since equal_var=False, doesn't assume pops have the same variance
    #gives identical results when sample sizes are variances are equal
#two-sided test for the H0 that 2 independent samples have identical average
    #doesn't indicate statistical significance
    #can't reject null hypothesis that outlier sign. impact on AA distribution
stat, p = ttest_ind(arr_T, arr_T2,equal_var=False)
print('stat=%.3f, p=%.3f' % (stat, p)) #stat=0.282, p=0.779
#tests if array came from normal distribution
    #they have differing, still not stat. sign. that outlier affecting distribution 
    #values for each may indicate T doesn't have a normal distribution
    #but comparing the two values, there isn't a large difference
stat2, p2=normaltest(arr_T) #pvalue=0.1584445119563328
stat3, p3=normaltest(arr_T2) #pvalue=0.22536615996733753
#can also be seen in kruskal test: 
    #determines if stat sign btwn two groups distribution
stat4, p4 = kruskal(arr_T, arr_T2)
print('stat=%.3f, p=%.3f' % (stat2, p2)) #stat=0.037, p=0.847

#makes a data frame so all the p-values can be viewed easiliy  
data={'T-Test': p,'Kruskal-Wallis': p4,
      'Normal Distribution W/ Outlier': p2,'Normal Distribution W/O Outlier': p3}
stats_df=pd.DataFrame(data, index=['p-value'])
stats_df.to_csv("Statistical Analysis.csv")
