#!/usr/bin/env python3
# coding: utf-8

"""
2021/01/18
author:guisen chen
email:thecgs001@foxmail.com
"""

import re
import sys
from collections import defaultdict

mrna_length_dict = defaultdict(list)
cds_length_dict = defaultdict(list)
exon_length_dict = defaultdict(list)

with open(sys.argv[1],'r') as f:
    for line in f:
        if line.startswith('#') or line == "\n": #remove '#' and '\n'
            continue
        line = line.rstrip().split("\t")
        if line[2] == "mRNA":
            mrna_id = re.search("ID=(.*?);",line[8]).group(1).strip()
            per_gene_per_mrna_length = int(line[4]) - int(line[3]) + 1
            mrna_length_dict[mrna_id].append(per_gene_per_mrna_length)
        elif line[2] == "CDS":
            mrna_id = re.search("Parent=(.*?);",line[8]+";").group(1)
            per_gene_per_cds_length = int(line[4]) - int(line[3]) + 1
            cds_length_dict[mrna_id].append(per_gene_per_cds_length)
        elif line[2] == "exon":
            mrna_id = re.search("Parent=(.*?);",line[8]+";").group(1)
            per_gene_per_exon_length = int(line[4]) - int(line[3]) + 1
            exon_length_dict[mrna_id].append(per_gene_per_exon_length)

#print(len(cds_length_dict))
#print(len(exon_length_dict))
#print(len(mrna_length_dict))

def Average_ORF_length():
    """
    intron_length = mrna_length - exon_length
    ORF_length = intron_length + CDS_length
    """
    temp = []
    for mrna_id,exon_length_list in exon_length_dict.items():
        intron_length = int(mrna_length_dict[mrna_id][0]) - int(sum(exon_length_list))
        temp.append(sum(cds_length_dict[mrna_id]) + intron_length)
    Average_ORF_length = sum(temp)/len(mrna_length_dict)
    return Average_ORF_length

def mrna_number():
    mrna_number = len(mrna_length_dict)
    return mrna_number

def Average_mrna_length():
    """
    sum(mrna_length)/mrna_number
    """
    temp = []
    for mrna_id,mrna_length in mrna_length_dict.items():
        temp.append(mrna_length[0])
    Average_mrna_length = sum(temp)/len(mrna_length_dict)
    return Average_mrna_length

def per_gene_Average_CDS_length():
    """
    sum(cds_length)/mrna_number
    """
    temp = []
    for mrna_id,cds_length in cds_length_dict.items():
        temp.append(sum(cds_length))
    per_gene_Average_CDS_length = sum(temp)/len(mrna_length_dict)
    return per_gene_Average_CDS_length

def Average_CDS_length():
    """
    sum(cds_length)/cds_number
    """
    temp = []
    number = []
    for mrna_id,cds_length in cds_length_dict.items():
        temp.append(sum(cds_length))
        number.append(len(cds_length))
    Average_CDS_length = sum(temp)/sum(number)
    return Average_CDS_length

def per_gene_Average_cds_number():
    """
    sum(per_gene_cds_number)/mrna_number
    """
    temp = []
    for gene_id,cds_length_list in cds_length_dict.items():
        temp.append(len(cds_length_list))
    per_gene_Average_cds_number = sum(temp)/len(mrna_length_dict)
    return per_gene_Average_cds_number


def per_gene_Average_exon_length():
    """
    sum(exon_length_per_gene)/mrna_number
    """
    temp = []
    for mrna_id,exon_length in exon_length_dict.items():
        temp.append(sum(exon_length))
    per_gene_Average_exon_length = sum(temp)/len(mrna_length_dict)
    return per_gene_Average_exon_length

def Average_exon_length():
    """
    sum(exon_length)/exon_number
    """
    temp = []
    number = []
    for mrna_id,exon_length in exon_length_dict.items():
        temp.append(sum(exon_length))
        number.append(len(exon_length))
    Average_exon_length = sum(temp)/sum(number)
    return Average_exon_length

def Average_exons_length_per_gene():
    """
    sum(per_gene_exon_all_length/per_gene_exon_number)/mrna_number
    """
    temp = []
    for mrna_id,exon_length_list in exon_length_dict.items():
        temp.append(sum(exon_length_list)/len(exon_length_list))
    Average_exons_length_per_gene = sum(temp)/len(mrna_length_dict)
    return Average_exons_length_per_gene

def per_gene_Average_exons_number():
    """
    sum(per_gene_exon_number)/mrna_number
    """
    temp = []
    for gene_id,exon_length_list in exon_length_dict.items():
        temp.append(len(exon_length_list))
    per_gene_Average_exons_number = sum(temp)/len(mrna_length_dict)
    return per_gene_Average_exons_number

def Average_intron_length():
    """
    intron_length = mrna_length - exon_length
    Average_intron_length = intron_length/mrna_number
    """
    temp = []
    for mrna_id,exon_length_list in exon_length_dict.items():
        intron_length = int(mrna_length_dict[mrna_id][0]) - int(sum(exon_length_list))
        temp.append(intron_length)
    Average_intron_length = sum(temp)/len(mrna_length_dict)
    return Average_intron_length

def cds_Average_intron_length():
    """
    intron_length = mrna_length - cds_length
    Average_intron_length = intron_length/mrna_number
    """
    temp = []
    for mrna_id,cds_length_list in cds_length_dict.items():
        intron_length = int(mrna_length_dict[mrna_id][0]) - int(sum(cds_length_list))
        temp.append(intron_length)
    cds_Average_intron_length = sum(temp)/len(mrna_length_dict)
    return cds_Average_intron_length


#print(f"mrna_number: {mrna_number()}") 
#print("Average_ORF_length: %.2f" %Average_ORF_length())
#print("Average_gene_length: %.2f" %Average_gene_length())
#print("Average_mrna_length: %.2f" %Average_mrna_length())
#print("per_gene_Average_CDS_length: %.2f" %per_gene_Average_CDS_length()) 
#print("Average_CDS_length: %.2f" %Average_CDS_length())
#print("per_gene_Average_cds_number: %.2f" %per_gene_Average_cds_number())
#print("per_gene_Average_exon_length: %.2f" %per_gene_Average_exon_length()) 
#print("Average_exon_length: %.2f" %Average_exon_length()) 
#print("Average_exons_length_per_gene: %.2f" %Average_exons_length_per_gene())
#print("per_gene_Average_exons_number: %.2f" %per_gene_Average_exons_number())   
#print("Average_intron_length: %.2f" %Average_intron_length())
#print("cds_Average_intron_length: %.2f" %cds_Average_intron_length())

print("Number of genes \t Average transcript length(bp) \t Average CDS length(bp) \t Average exon length(bp) \t Average intron length(bp) \t Average exons per gene") #head
print("%d \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f" %(mrna_number(),Average_mrna_length(),per_gene_Average_CDS_length(),Average_exon_length(),Average_intron_length(),per_gene_Average_exons_number()))
