#!/usr/bin/env python
# coding: utf-8

# # Rosalind Problems

# ### These problems are found at rosalind.info, a platform for learning bioinformatics through problem solving. The problems below deal with solving bioinformatics problems using python, while additional problems dealing with sequencing can be found at https://github.com/agalianese/Sequencing.

# In[50]:


import Bio
from Bio import SeqIO
from Bio import Seq
import pandas as pd
import numpy as np


# In[15]:


#Computing GC Content

#input file and it's type (fasta, genbank, etc.)
file_in = open("C:/Users/mc/Downloads/rosalind_gc.txt", 'r')
file_type = 'fasta'


#create a function
def GC_percent(file_in, file_type):
    """This function takes in a file and its type (fasta, db, genbank) and returns
    which sequence id has the highest GC percent and the percent"""
    
    #beginning percent
    best_percent = 0
    
    #filter through the file by 'record' or sequence
    for record in SeqIO.parse(file_in, file_type):
        
        #count of how much of the sequence is a G or C, versus A or T
        GC_count = 0
        
        #for loop of each base per sequence
        for seq in record.seq:
            
            #if the base is a G or C, it increases the count
            if seq == 'G' or seq == 'C':
                GC_count += 1
                
            #logs the percent for each sequence
            percent = (GC_count / len(record.seq))
            
        #if the current sequence has a greater GC percent, then it becomes the new 'best' and its percent and id are logged 
        if percent > best_percent:
            best_percent = percent
            best_id = record.id
            
    #function finishes by printing the best id name and its associated percent GC
    print(best_id)
    print(str(100 * best_percent) + '%')
    
    
#use the function by inputing the file location, and its type    
GC_percent(file_in, file_type)


# In[45]:


#Counting Point Mutations

#input the original and mutated string
original = 'CGGGACCTTCCGCTTTCGATTATTTGACCTATCGACCAAGGGTGACTAGTTCTAATTTTTTCCGGGGAGCGATAATATGCGTCAAGCCTACAAAGTTACAGCGCAAAACCGTCACAGTACCCTAGGTGCTGAACCAAAATTATGACCAATTATGTTACGTGCTGATAGTGGCGACATAGCACTGGGTTAGATGAGACGGACGTTCGGCACAACGGACAGGTATAAGATGCCCTCTATCCCTCGCTTACAATGGCCGTTGTTGATGCGAGTTACTGTACGTAGCGACACGGAGCAGTTCCGCGAGACCTACAGATTCTGTACCCTTTACCGCTGACGTTCATCGTATAAATCGGCAATTTGTACTTTACACCTCTGTAATTTATCATGTGCACCATATCCCATTGACGCGTCATGTTAGATAAGGCCTACCTGGTGGAGCACATTCTCGCCAATCTTGATATTTGCCTCCTTGCAGGAACCTAGCTATATAAATCCTGTGCAAAAAATAGTTCTGCGCGTAGCGGTAACGAACTGTGAACTTGAACATAACAAGGTTCGATGCGCCTAAGGGCTCTGCTCACAGCTAGTTGGATACCTAAAGTAGCGATGATAAAATTAGTAGGTGATTCCTTGTCGGAATAGACGGTTTTACCCACGGCTCAAACAATCTAGGAGTCAATGCTGCTAGTCGGAAGTCCGTTGCTTCGTGAAATGATGCGCACACGATCCCTAGGCCAGAATTACATTCCGTCTATGGTTTTACTTGCTTGTGGACTATGTCTATAACTCATACCTTCGGTTCATATTGTTGTAACTTCTCCAGGCACTATTTTGGCCGACACCACAAGCACTGGTAGAGCCAGACAGTTATGTGTGGAGTCCCTGAGAAGCATTAATGAGTTGCCGTAATTGCTGATCCCCCTCTTCCACTACATCGGCCTGAACGGCGGAGATAC'
copy = 'TAGAACATTCCAGCGTAGTCAAGATCCAGCTCCCGCGCTTTGACTCCAATACTAAGTCAATCCTGTCTTCGATACTCGTGGTATAGGCCACATAGCCTCTGGTCGTAGCTCTTACCAATAAATGGGAGCTGACCCACAGATACCGCCCGATAGGATCCGGACTGATCACCGTGACATACCGTTAGCTTTAATGTTACGGACGGACCGCGGGTTTCTAACTTATTGACTCCCCTCCGTCCCTCCCACGCTATCGCTGTAGGTGTTTACCGTGAGTTGTTATTTCGGTACGGGGAGGTTCCTGGAGACGGAAGGAAAACGAACCATTTCCACTAGCCCCTCATGGACTTAACGCTCCATATGTTTTTCGCTCCCCTGAATTACATTAAGGCTGCAATCTGTTTTTGATGCAGGGCGTAATTTCAGATCTAGCCGGTCGCCCATATTCTCGCCAATTTGGTGTGTAGTCTTACACCATGCTAATGTCTGTACTAATCGTCTTCGCAAAATCTACCTGCGTGGATGTGTTTTGGACTGGAAACTTGATAATCTCCATCCACGCTCCACCTAATGATTCTGCGAACCGAGATTGCTCGAAATAAATGAGTGCGAGAAAACTTCATTGGTTACTAGCCGTCCCCACTAAAGTTAGTAAACCCTTATCTACCAACATAATCCTTCCCGGAGCTGTAGAGTTTTCTGCCCTTGCATGCAAGGATGCGCCTACCCCTCTTCGCTGGGAAGATAGAGTCTGAAAGGCTATAACTTTCTTGCGGAATTGGAATCCCACTTATACGAACCGTTTCACGGGTTGCGGCTCATCCCGGGTCTAGAGCGCCAGACCCAATATCCCATGGACCAGGGATCGCGCTACGTGGTGGGTCGTTTAAAAGCGGTCCTGAATTAACACTCACGGCGAACCGCCTACTGAATTTCATTGCGAACACGGGGCTAGAATT'

def Point_Mutation_Counter(original, copy):
    """Function that takes in an original sequence and its mutated copy, and reports how many differences there are.
    Also reports which basepairs are being mutated.
    """
    #begin counters
    mismatch_counter = 0
    a_count = 0
    t_count = 0
    g_count = 0
    c_count = 0

    #iterate over the original and the copy sequence
    for i in range(len(original)):

        #if the original base does not match the copy, increase the counter
        if original[i] != copy[i]:
            mismatch_counter += 1
            
            #computes which bases are being mutated 
            if original[i] == 'A':
                a_count += 1
            elif original[i] == 'T':
                t_count += 1
            elif original[i] == 'G':
                g_count += 1
            elif original[i] == 'C':
                c_count += 1
            else:
                print("Other base pair reported")
                break

    #return pandas dataframe composed of the amount of mismatches, and which basepairs were mutated
    return_df = pd.DataFrame(
    {'Mismatches': mismatch_counter,
     "Number of A's Mutated" : a_count, 
     "Number of T's Mutated" : t_count,
     "Number of G's Mutated" : g_count,
     "Number of C's Mutated" : c_count,}, index = [1])
    
    #print the final count
    return(return_df)

Point_Mutation_Counter(original, copy)


# In[17]:


#dicionary for keys for the rna and values for the aa counterpart to be used below

RNA_to_AA= {"UUU" : "F","CUU" : "L","AUU" : "I","GUU" : "V","UUC" : "F","CUC" : "L","AUC" : "I","GUC" : "V","UUA" : "L","CUA" : "L","AUA" : "I",     
"GUA" : "V","UUG" : "L","CUG" : "L","AUG" : "M","GUG" : "V","UCU" : "S","CCU" : "P","ACU" : "T","GCU" : "A","UCC" : "S","CCC" : "P",
"ACC" : "T","GCC" : "A","UCA" : "S","CCA" : "P","ACA" : "T","GCA" : "A","UCG" : "S","CCG" : "P","ACG" : "T","GCG" : "A","UAU" : "Y",
"CAU" : "H","AAU" : "N","GAU" : "D","UAC" : "Y","CAC" : "H","AAC" : "N","GAC" : "D","UAA" : "Stop","CAA" : "Q","AAA" : "K","GAA" : "E","UAG" : "Stop",
"CAG" : "Q","AAG" : "K","GAG" : "E","UGU" : "C","CGU" : "R","AGU" : "S","GGU" : "G","UGC" : "C","CGC" : "R","AGC" : "S","GGC" : "G","UGA" : "Stop",
"CGA" : "R","AGA" : "R","GGA" : "G","UGG" : "W","CGG" : "R","AGG" : "R","GGG" : "G"}


# In[21]:


#Translating RNA into Protein

#read in the RNA strand to be translated
RNA_strand = open("C:/Users/mc/Downloads/rosalind_prot.txt", 'r')
RNA_strand = RNA_strand.read()


#RNA to Protein Function
def RNA_to_Protein(RNA_strand):
    '''This function will take in a RNA strand and convert it into its protein counterpart using Amino Acid codons'''
    
    #begin with an empty protein string
    protein = ''
    
    #codons are a 3 nucleotide sequence that can be translated into a single letter
    #this will iterate over the RNA strand in 3 letter chunks
    for i in range(0,len(RNA_strand),3):
        
        #each codon is a 3 letter long section of the RNA strnd
        codon = RNA_strand[i:i+3]
        
        #using the RNA to AA dictionary, the amino acid is the translated codon
        amino_acid = RNA_to_AA[codon]
        
        #if the amino acid is not a stop codon, it will add the single letter code to the protein strand
        if amino_acid != 'Stop':
            protein = protein + amino_acid
            
        #if the codon is a stop codon, the loop will break
        else: break
            
    #function returns the translated protein
    return(protein)

#returns the translated protein as a string
RNA_to_Protein(RNA_strand)


# In[25]:


#Finding a Motif in DNA

#motif and sequence to be analyzed
motif = 'CGTGAGTGCTCCCGTGAGTTGTAGCGTGAGTATGCGTGAGTCGTGAGTGGTGGCGTGAGTCACCACGTGAGTCAAACGTGAGTGCGTGAGTCGTGAGTCGTGAGTTCGTGAGTCCACCGTGAGTGCCATCTAAGGCAGTGGCGTGAGTCCACGTGAGTTTCGTGAGTACTTCCGTGAGTTCGTGAGTTCGTGAGTCGTGAGTCCGTGAGTCGTGAGTGTCACGTGAGTGCCGTGAGTCGTGAGTCGTGAGTCGTGAGTGCCGTGAGTACAGCGTGAGTCTGCGTGAGTTGCGTGAGTCGTGAGTAGCCGTGAGTCGTGAGTCCCTAAAGCGTGAGTCGTGAGTCCCGTGAGTGCGCAAACCGTGAGTACCGCGTGAGTACGTGAGTCCGTGAGTCGTGAGTCGTGAGTCGTGAGTCCGTGAGTATCGTGAGTGGACGTGAGTCGACGTGAGTAGTCGTGAGTGAAGCGTGAGTCGTGAGTCGTGAGTAGTCGTGAGTGCGTGAGTTTCTCGTGAGTCGTGAGTGCGTGAGTGGATTGGCGTGAGTCGTGAGTACGTGAGTGTCGGCAGTCCGTGAGTCGTGAGTCGTGAGTAGCTCGTGAGTAGGGTGCCGTGAGTAAGCACAACGCGTGAGTGCGTGAGTGTTCGTGAGTGGTGATCGTGAGTTCGTGAGTCGTGAGTCGTGAGTCACTAGTCCGTGAGTCGTGAGTGTGAGAGGGCGTGAGTCGTGAGTCGGTTCTCGTGAGTTACGTGAGTTCGTGAGTCGTGAGTCGTGAGTCACGTGAGTCCCGTGAGTTCACGTGAGTCGCGTGAGTACTGCGTGAGTCGCGTGAGTCGTGAGTCCCGTGAGTACGTGAGTCGTGAGTCCCATGCAGGCAAT'
seq = 'CGTGAGTCG'

def Motif_Finder(motif, seq):
    """Identifies a chosen sequence in DNA, RNA, or protein motif, and reports the location it identifies the sequence"""
    
    #creates an empty string to place the locations
    print_str = ''
    
    #iterates over the initial motif
    for i in range(0, len(motif)):
        
        #if the sequence is in the motif, it will log the location as the beginning basepair location
        if motif[i:i+len(seq)] == seq:
            
            #changes the starting location from 0 to 1 for easier counting
            index = i + 1
            
            #adds to the location string 
            print_str = print_str + " " + str(index)
            
    #returns the string full of the location indexes
    return(print_str)    


Motif_Finder(motif, seq)


# In[49]:


#Introduction to the Bioinformatics Armory

my_seq = "TATTGGGATCTCGGCATTATTGTGGACCACCCCATGGTCCACCTTTAAGCCGCCTTAGGTGATCGTCAAAATGTGAATGACCTGCTGCGATGGGTAGACGCCAGTCCAGCAATCAGGGAGACGAAAGTTGCTCGACTGACCCGTGTGGGTTACCCGTTGCAACAACCTCGCCAAAGGACTTCCGCGGCCAGACTAAATGCCTACCGTTTAGTGCTCCTCTTGGTGCACGCTGCTGAGAAAATTGCTCATGTCCTGGTGGCTGGGCCCCCGACAGGAGGTTAGGGCTTAATATAATGCCGCGTAGGGGTTGTATGTGGCCCATTATAATTGAACTCCTCGCAAGCTGAAAGTAATGCGTACCTCGCGTTGATGTTGTGTCGCAGCACCCGTTGTCTTCTCCGTACTGAGCAGTGCGGTTAGGCCGCGAAGCATTTTCATCAGTGATCGCCTAAATATATTTTGTATCTCTATTTACTGCCTCGGTGGGTCTTGTCCTGCTATGGATAAACGAGACGGTCTCCCATGAAAAGCGTTTTATGCCAGACTTCTACGGCCCCGCCTTATCAGCTTGCTGGGACGTCTCGCAGCACCTTTTCTGGGAGCCGTTCAGTTCAGTTTTCCTCGCGCTTAATGCGCCGATGAGAGGCTGAGTAACACATCGTATGTATCGTACGCTATTGGTGTTCTGATCGCGTGAGCCCGGGTGTTTGCCTAGTCATACGAAAGACGGAAGGGAGGGGCGCTTGTGGGTCGGGCAAGGATATCCATCAAGAGAGCACTGGTCCCCAAGGCTGCAACTGTATTACATAATGGATTTGTGAGATAAGAATCTACTACCAATACTTACAAGGTGCGGGCTGCTAACATTGCTGTAGGCGACGAACTGCGCGGTGGACGACCTCAGCATGCGCTCCAGGTGTGCCGATTTT"
    
def Base_Counter(my_seq):
    """Counts the number of each base pair in a sequence"""
    
    #counts the number of times each base pair is present in a sequence
    acount = str(my_seq.count("A"))
    ccount = str(my_seq.count("C"))
    gcount = str(my_seq.count("G"))
    tcount = str(my_seq.count("T"))

    #returns a pandas datafrane of the counts with an index starting at 1
    return_df = pd.DataFrame(
    {"A's" : acount, 
     "T's" : tcount,
     "G's" : gcount,
     "C's" : ccount,}, index = [1])
    
    #return the final counts
    return(return_df)


Base_Counter(my_seq)


# In[ ]:




