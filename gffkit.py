#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import re
import pyfastx
import argparse
import numpy as np
from collections import defaultdict

def seq_compliment(seq):
    intab = "ATGCatgc"
    outtab = "TACGtacg"
    trantab = str.maketrans(intab, outtab)
    return seq.translate(trantab)

def seq_translate(DNA_seq):
    DNA_seq = DNA_seq.upper()
    Genetic_code = {'TTT': 'F', 'CTT': 'L', 'ATT': 'I', 'GTT': 'V',                      'TTC': 'F', 'CTC': 'L', 'ATC': 'I', 'GTC': 'V',                      'TTA': 'L', 'CTA': 'L', 'ATA': 'I', 'GTA': 'V',                      'TTG': 'L', 'CTG': 'L', 'ATG': 'M', 'GTG': 'V',                      'TCT': 'S', 'CCT': 'P', 'ACT': 'T', 'GCT': 'A',                      'TCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',                      'TCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A',                      'TCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A',                      'TAT': 'Y', 'CAT': 'H', 'AAT': 'N', 'GAT': 'D',                      'TAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',                      'TAA': '*', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E',                      'TAG': '*', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',                      'TGT': 'C', 'CGT': 'R', 'AGT': 'S', 'GGT': 'G',                      'TGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G',                      'TGA': '*', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',                      'TGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'}
    protein = ''
    for site in range(0, len(DNA_seq), 3):
        if "N" in DNA_seq[site:site+3]:
            protein += "X"
        if "N" not in DNA_seq[site:site+3]:
            if len(DNA_seq[site:site+3]) == 3:
                protein += Genetic_code[DNA_seq[site:site+3]]
            if len(DNA_seq[site:site+3]) != 3:
                pass
    return protein

def seq(args):
    # Open a geneome file, from genome file create a index file, and create a grep dict.
    genome_input = pyfastx.Fastx(args.fastx)
    genome_fa = defaultdict(str)
    for name,seq,comment in genome_input:
        genome_fa[name] = seq

    # Open a gff file, and parse it.
    gff_input = open(args.gff,'r')
    Info = defaultdict(list)
    gene_like_list = []#like
    gene2RNA = defaultdict(list)#like
    gene2mRNA = defaultdict(list)
    exon = defaultdict(list)
    intron = defaultdict(list)
    CDS = defaultdict(list)
    for line in gff_input:
        if line.startswith('#') or line == "\n": #remove '#' and '\n'
            continue
        line_list = line.split("\t")
        # API-like
        if args.feature == 'gene-like':
            if line_list[2] == args.str:
            #such as "region","pseudogene",'cDNA_match'
                ID = re.search("ID=(.*?)[;,\n]",line_list[8]).group(1).strip()
                Info[ID] = line_list
                gene_like_list.append(ID)
        if args.feature == "RNA-like":
            if line_list[2] == args.str:
            #such as "transcript","lnc_RNA","tRNA","snoRNA","snRNA","rRNA","guide_RNA","V_gene_segment","C_gene_segment"
                RNA_ID = re.search("ID=(.*?)[;,\n]",line_list[8]).group(1).strip()       
                Info[RNA_ID] = line_list
                ID = re.search("Parent=(.*?)[;,\n]",line_list[8]).group(1).strip()
                if ID in gene2RNA:
                    gene2RNA[ID].append(RNA_ID)
                else:
                    gene2RNA[ID] = [RNA_ID]
            if line_list[2] == 'exon':
                exon_ID = re.search("ID=(.*?)[;,\n]",line_list[8]).group(1).strip()
                RNA_ID = re.search("Parent=(.*?)[;,\n]",line_list[8]).group(1).strip()
                Info[exon_ID] = line_list
                if RNA_ID in exon:
                    exon[RNA_ID].append(exon_ID)
                else:
                    exon[RNA_ID] = [exon_ID]    
        else:
            if line_list[2] == 'gene':
                gene_ID = re.search("ID=(.*?)[;,\n]",line_list[8]).group(1).strip()
                Info[gene_ID] = line_list       
            if line_list[2] == 'mRNA':
                mRNA_ID = re.search("ID=(.*?)[;,\n]",line_list[8]).group(1).strip()       
                Info[mRNA_ID] = line_list
                gene_ID = re.search("Parent=(.*?)[;,\n]",line_list[8]).group(1).strip()
                if gene_ID in gene2mRNA:
                    gene2mRNA[gene_ID].append(mRNA_ID)
                else:
                    gene2mRNA[gene_ID] = [mRNA_ID]            
            if line_list[2] == 'exon':
                exon_ID = re.search("ID=(.*?)[;,\n]",line_list[8]).group(1).strip()
                mRNA_ID = re.search("Parent=(.*?)[;,\n]",line_list[8]).group(1).strip()
                Info[exon_ID] = line_list
                if mRNA_ID in exon:
                    exon[mRNA_ID].append(exon_ID)
                else:
                    exon[mRNA_ID] = [exon_ID]            
            if line_list[2] == 'intron':
                intron_ID = re.search("ID=(.*?)[;,\n]",line_list[8]).group(1).strip()
                mRNA_ID = re.search("Parent=(.*?)[;,\n]",line_list[8]).group(1).strip()
                Info[intron_ID] = line_list
                if mRNA_ID in intron:
                    intron[mRNA_ID].append(intron_ID)
                else:
                    intron[mRNA_ID] = [intron_ID]            
            if line_list[2] == 'CDS':
                mRNA_ID = re.search("Parent=(.*?)[;,\n]",line_list[8]).group(1).strip()
                if mRNA_ID in CDS:
                    CDS[mRNA_ID].append(line_list)
                else:
                    CDS[mRNA_ID] = [line_list]    
    gff_input.close()

    # If feature is "gene-like", the sequences were extracted according to gff file 1, 4 ,5 and 7 column.
    if args.feature == 'gene-like':
        fasta = open(args.out,'w+')
        for gene_ID in gene_like_list:
            #print(f">{gene_ID}.gene")
            fasta.write(f">{gene_ID}.{args.str}\n")
            start = int(Info[gene_ID][3:5][0])
            end = int(Info[gene_ID][3:5][1])
            chr_ID = Info[gene_ID][0]
            if Info[gene_ID][6] == '-':
                seq = genome_fa[chr_ID][start-1:end]                    
                seq = seq_compliment(seq)[::-1]
                #print(seq)
                fasta.write(seq+"\n")
            if Info[gene_ID][6] == '+':
                seq = genome_fa[chr_ID][start-1:end]
                #print(seq)
                fasta.write(seq+"\n")
        fasta.close()

    # If feature is "RNA-like", the sequences were extracted according to gff file 1, 4 ,5 and 7 column.
    if args.feature == 'RNA-like':
        fasta = open(args.out,'w+')
        for gene_ID,RNA_ID in gene2RNA.items():
            for r in RNA_ID:
                #print(f">{m}.exon")
                fasta.write(f">{r}.{args.str}.all.exon\n")
                pos = np.sort(np.array([Info[e][3:5] for e in exon[r]]).astype(int),axis=0)
                chr_ID = Info[r][0]
                if Info[r][6] == '-':
                    pos = pos[np.argsort(-pos[:,0])]
                    for start,end in pos:
                        seq = genome_fa[chr_ID][start-1:end]
                        seq = seq_compliment(seq)[::-1]
                        #print(seq)
                        fasta.write(seq+"\n")
                if Info[r][6] == '+':
                    for start,end in pos:
                        seq = genome_fa[chr_ID][start-1:end]
                        #print(seq)
                        fasta.write(seq+"\n")
        fasta.close()

    # If feature is "gene", the sequences were extracted according to gff file 4,5 and 7 column.
    if args.feature == 'gene':
        fasta = open(args.out,'w+')
        for gene_ID,mRNA_ID in gene2mRNA.items():
            #print(f">{gene_ID}.gene")
            fasta.write(f">{gene_ID}.gene\n")
            start = int(Info[gene_ID][3:5][0])
            end = int(Info[gene_ID][3:5][1])
            chr_ID = Info[gene_ID][0]
            if Info[gene_ID][6] == '-':
                seq = genome_fa[chr_ID][start-1:end]                    
                seq = seq_compliment(seq)[::-1]
                #print(seq)
                fasta.write(seq+"\n")
            if Info[gene_ID][6] == '+':
                seq = genome_fa[chr_ID][start-1:end]
                #print(seq)
                fasta.write(seq+"\n")
        fasta.close()

    # If feature is "mRNA", the sequences were extracted according to gff file 4,5 and 7 column.
    if args.feature == 'mRNA':
        fasta = open(args.out,'w+')
        for gene_ID,mRNA_ID in gene2mRNA.items():
            for m in mRNA_ID:
                #print(f">{m}.mRNA")
                fasta.write(f">{m}.mRNA\n")
                start = int(Info[m][3:5][0])
                end = int(Info[m][3:5][1])
                chr_ID = Info[m][0]
                if Info[m][6] == '-':
                    seq = genome_fa[chr_ID][start-1:end]                    
                    seq = seq_compliment(seq)[::-1]
                    #print(seq)
                    fasta.write(seq+"\n")
                if Info[m][6] == '+':
                    seq = genome_fa[chr_ID][start-1:end]
                    #print(seq)
                    fasta.write(seq+"\n")
        fasta.close()

    # If feature is "exon", the sequences were extracted according to gff file 4,5 and 7 column. Name is mRNA id, and seq of every row is an exon.
    if args.feature == 'exon':
        fasta = open(args.out,'w+')
        for gene_ID,mRNA_ID in gene2mRNA.items():
            for m in mRNA_ID:
                #print(f">{m}.exon")
                fasta.write(f">{m}.exon\n")
                pos = np.sort(np.array([Info[e][3:5] for e in exon[m]]).astype(int),axis=0)
                chr_ID = Info[m][0]
                if Info[m][6] == '-':
                    pos = pos[np.argsort(-pos[:,0])]
                    for start,end in pos:
                        seq = genome_fa[chr_ID][start-1:end]
                        seq = seq_compliment(seq)[::-1]
                        #print(seq)
                        fasta.write(seq+"\n")
                if Info[m][6] == '+':
                    for start,end in pos:
                        seq = genome_fa[chr_ID][start-1:end]
                        #print(seq)
                        fasta.write(seq+"\n")
        fasta.close()

    # if feature is "intron", the sequences were extracted according to gff file 4,5 and 7 column. Name is mRNA id, and seq of every row is a intron.
    if args.feature == 'intron':
        fasta = open(args.out,'w+')
        for gene_ID,mRNA_ID in gene2mRNA.items():
            for m in mRNA_ID:
                #print(f">{m}.intron")
                fasta.write(f">{m}.intron\n")
                pos = np.sort(np.array([Info[i][3:5] for i in intron[m]]).astype(int),axis=0)
                chr_ID = Info[m][0]
                if Info[m][6] == '-':
                    pos = pos[np.argsort(-pos[:,0])]
                    for start,end in pos:
                        seq = genome_fa[chr_ID][start-1:end]
                        seq = seq_compliment(seq)[::-1]
                        #print(seq)
                        fasta.write(seq+"\n")
                if Info[m][6] == '+':
                    for start,end in pos:
                        seq = genome_fa[chr_ID][start-1:end]
                        #print(seq)
                        fasta.write(seq+"\n")
        fasta.close()

    # If feature is "upstream", the sequences were extracted upstream of mRNA.
    if args.feature == 'upstream':
        fasta = open(args.out,'w+')
        length = int(args.length)
        for gene_ID,mRNA_ID in gene2mRNA.items():
            for m in mRNA_ID:
                fasta.write(f">{m}.upstream\t{length}bp\n")
                chr_ID = Info[m][0]
                if Info[m][6] == '-':
                    pos = int(Info[gene_ID][3:5][1])
                    if len(genome_fa[chr_ID]) < pos+1+length:
                        seq = genome_fa[chr_ID][pos+1:]
                        seq = seq_compliment(seq)[::-1]
                        fasta.write(seq+"\n")
                    if len(genome_fa[chr_ID]) >= pos+1+length:
                        seq = genome_fa[chr_ID][pos+1:pos+1+length]
                        seq = seq_compliment(seq)[::-1]
                        fasta.write(seq+"\n")              
                if Info[m][6] == '+':
                    pos = int(Info[gene_ID][3:5][0])
                    if 0 <= pos-1-length:
                        seq = genome_fa[chr_ID][pos-1-length:pos-1]
                        fasta.write(seq+"\n")
                    if 0 > pos-1-length:
                        seq = genome_fa[chr_ID][:pos-1]
                        fasta.write(seq+"\n")                    
        fasta.close()

    # If feature is "CDS", the sequences were extracted according to gff file 4, 5, 7 and 8 column.
    # If feature is "pep", the CDS translates to protein
    if args.feature == 'pep' or args.feature =='CDS':
        fasta = open(args.out,'w+')
        for gene_ID,mRNA_ID in gene2mRNA.items():
            for m in mRNA_ID:
                #print(f">{m}.")
                fasta.write(f"\n>{m}\n")
                pos_size = np.array([c[3:5]+[c[7]] for c in CDS[m]])
                pos_size = pos_size.astype(int)
                chr_ID = Info[m][0]
                if Info[m][6] == '+':
                    pos_size = pos_size[np.argsort(pos_size[:,0])]
                    pos1 = pos_size[0,0:2]
                    pos = pos_size[1:,0:2]
                    size = int(pos_size[0,2])
                    if size == 0:
                        seq = genome_fa[chr_ID][int(pos1[0])-1:int(pos1[1])]
                        #print(seq)
                        fasta.write(seq)
                    if size == 1:
                        seq = genome_fa[chr_ID][int(pos1[0]):int(pos1[1])]
                        #print(seq)
                        fasta.write(seq)
                    if size == 2:
                        seq = genome_fa[chr_ID][int(pos1[0])+1:int(pos1[1])]
                        #print(seq)
                        fasta.write(seq)
                    for start,end in pos:
                        seq = genome_fa[chr_ID][start-1:end]
                        #print(seq)
                        fasta.write(seq)              
                if Info[m][6] == '-':
                    pos_size = pos_size[np.argsort(-pos_size[:,0])]
                    pos1 = pos_size[0,0:2]
                    pos = pos_size[1:,0:2]
                    size = int(pos_size[0,2])
                    if size == 0:
                        seq = genome_fa[chr_ID][int(pos1[0])-1:int(pos1[1])]
                        seq = seq_compliment(seq)[::-1]
                        #print(seq)
                        fasta.write(seq)
                    if size == 1:
                        seq = genome_fa[chr_ID][int(pos1[0]):int(pos1[1])]
                        seq = seq_compliment(seq)[::-1]
                        #print(seq)
                        fasta.write(seq)
                    if size == 2:
                        seq = genome_fa[chr_ID][int(pos1[0])+1:int(pos1[1])]
                        seq = seq_compliment(seq)[::-1]
                        #print(seq)
                        fasta.write(seq)
                    for start,end in pos:
                        seq = genome_fa[chr_ID][start-1:end]
                        seq = seq_compliment(seq)[::-1]
                        #print(seq)
                        fasta.write(seq)
        fasta.write("\n")
        fasta.close()
        if args.feature == 'pep':
            import os
            fa = pyfastx.Fastx(args.out)
            os.remove(args.out)
            fasta = open(args.out,'a')
            for name, seq, comment in fa:
                fasta.write(f">{name}\n{seq_translate(seq)}\n")
            fasta.close()

    print('finished')

def stat(args):
    mrna_length_dict = defaultdict(list)
    cds_length_dict = defaultdict(list)
    exon_length_dict = defaultdict(list)
    with open(args.gff,'r') as f:
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
    print(f"mrna_number: {mrna_number()}") 
    #print("Average_ORF_length: %.2f" %Average_ORF_length())
    #print("Average_gene_length: %.2f" %Average_gene_length())
    print("Average_mrna_length: %.2f" %Average_mrna_length())
    print("per_gene_Average_CDS_length: %.2f" %per_gene_Average_CDS_length()) 
    #print("Average_CDS_length: %.2f" %Average_CDS_length())
    #print("per_gene_Average_cds_number: %.2f" %per_gene_Average_cds_number())
    #print("per_gene_Average_exon_length: %.2f" %per_gene_Average_exon_length()) 
    print("Average_exon_length: %.2f" %Average_exon_length()) 
    #print("Average_exons_length_per_gene: %.2f" %Average_exons_length_per_gene())
    #print("per_gene_Average_exons_number: %.2f" %per_gene_Average_exons_number())   
    print("Average_intron_length: %.2f" %Average_intron_length())
    print("cds_Average_intron_length: %.2f" %cds_Average_intron_length())    

    #print("Number of genes \t Average transcript length(bp) \t Average CDS length(bp) \t Average exon length(bp) \t Average intron length(bp) \t Average exons per gene") #head
    #print("%d \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f" %(mrna_number(),Average_mrna_length(),per_gene_Average_CDS_length(),Average_exon_length(),Average_intron_length(),per_gene_Average_exons_number()))
def filter(args):
    def parse_gff(file_gff):
        gene_ID_dict = defaultdict(list) # all mRNA id corresponding to gene id 
        ID_info_dict = defaultdict(list) # info of gene, mRNA, CDS and exon
        CDS_length_dict = defaultdict(list) # all CDS length corresponding to mRNA id 
        the_longest_mRNA_ID_dict = defaultdict(str)

        for line in open(file_gff,'r'):

            if line.startswith('#') or line == "\n": #remove '#' and '\n'
                continue
            line_list = line.split("\t")

            if line_list[2] == "gene": 
                gene_ID = re.search("ID=(.*?)[;,\n]",line_list[8]).group(1).strip()
                ID_info_dict[gene_ID] = line_list[0:8] + [f'ID={gene_ID};']

            if line_list[2] == "mRNA" or line_list[2] == "transcript": 
                mRNA_ID = re.search("ID=(.*?)[;,\n]",line_list[8]).group(1)
                gene_ID = re.search("Parent=(.*?)[;,\n]",line_list[8]).group(1).strip()
                if gene_ID in gene_ID_dict:
                    gene_ID_dict[gene_ID].append(mRNA_ID)
                else:
                    gene_ID_dict[gene_ID] = [mRNA_ID]
                ID_info_dict[mRNA_ID].append(line_list[0:8] + [f'ID={mRNA_ID};Parent={gene_ID};'])

            if line_list[2] == "CDS": 
                CDS_ID = re.search("ID=(.*?)[;,\n]",line_list[8]).group(1)
                mRNA_ID = re.search("Parent=(.*?)[;,\n]",line_list[8]).group(1).strip()
                ID_info_dict[mRNA_ID].append(line_list[0:8] + [f'ID={CDS_ID};Parent={mRNA_ID};'])
                CDS_length = int(line_list[4]) - int(line_list[3]) + 1
                CDS_length_dict[mRNA_ID].append(CDS_length)

            if line_list[2] == "exon":
                exon_ID = re.search("ID=(.*?)[;,\n]",line_list[8]).group(1)
                mRNA_ID = re.search("Parent=(.*?)[;,\n]",line_list[8]).group(1).strip()
                ID_info_dict[mRNA_ID].append(line_list[0:8] + [f'ID={exon_ID};Parent={mRNA_ID};'])

        for gene_ID,mRNA_ID_list in gene_ID_dict.items():
            i = 0
            for mRNA_ID in mRNA_ID_list:
                CDS_length = sum(CDS_length_dict[mRNA_ID])
                if CDS_length > i:
                    i = CDS_length
                    the_longest_mRNA_ID_dict[gene_ID] = mRNA_ID

        gff_file = open(args.out,'a')
        gff = csv.writer(gff_file,delimiter='\t')
        for gene_ID, mRNA_ID in the_longest_mRNA_ID_dict.items():
            for line in ([ID_info_dict[gene_ID]]+ID_info_dict[mRNA_ID]):
                gff.writerow(line)
        gff_file.close()
        return print('finished')
    parse_gff(args.gff)

def add(args):
    gff_input = open(args.gff,'r')
    gff_output = open(args.out,'w')

    Info = defaultdict(list)
    gene2mRNA = defaultdict(list)
    CDS = defaultdict(list)
    exon = defaultdict(list)

    for line in gff_input:
        if line.startswith('#') or line == "\n": #remove '#' and '\n'
            continue

        line_list = line.split("\t")
        #print(line_list)
        if line_list[2] == 'gene':
            gene_ID = re.search("ID=(.*?)[;,\n]",line_list[8]).group(1).strip()
            Info[gene_ID] = line_list

        if line_list[2] == 'mRNA':
            mRNA_ID = re.search("ID=(.*?)[;,\n]",line_list[8]).group(1).strip()       
            Info[mRNA_ID] = line_list
            gene_ID = re.search("Parent=(.*?)[;,\n]",line_list[8]).group(1).strip()
            if gene_ID in gene2mRNA:
                gene2mRNA[gene_ID].append(mRNA_ID)
            else:
                gene2mRNA[gene_ID] = [mRNA_ID]

        if line_list[2] == 'CDS':
            mRNA_ID = re.search("Parent=(.*?)[;,\n]",line_list[8]).group(1).strip()
            if mRNA_ID in CDS:
                CDS[mRNA_ID].append(line_list)
            else:
                CDS[mRNA_ID] = [line_list]

        if line_list[2] == 'exon':
            exon_ID = re.search("ID=(.*?)[;,\n]",line_list[8]).group(1).strip()
            mRNA_ID = re.search("Parent=(.*?)[;,\n]",line_list[8]).group(1).strip()
            Info[exon_ID] = line_list
            if mRNA_ID in exon:
                exon[mRNA_ID].append(exon_ID)
            else:
                exon[mRNA_ID] = [exon_ID]

    gff_input.close()

    #print(Info)
    #print(gene2mRNA)
    #print(CDS)
    #print(exon)

    for gene_ID, mRNA_ID in gene2mRNA.items():
        #print(Info[gene_ID])
        gff_output.write("\t".join(Info[gene_ID]))
        for m in mRNA_ID:        
            temp = []
            #print(Info[m])
            gff_output.write("\t".join(Info[m]))
            for c in CDS[m]:
                gff_output.write("\t".join(c))            
            for e in exon[m]:            
                #print('\t'.join(Info[e]))
                gff_output.write("\t".join(Info[e]))
                temp.append(int(Info[e][3]))
                temp.append(int(Info[e][4]))
            temp = sorted(temp)
            if len(temp) <= 2:
                gff_output.write('\n')
                continue
            if len(temp) > 2 : # have Iron
                l = int(len(temp)/2)
                del(temp[0])
                temp.pop()
                start = (np.array(temp[::2]) + 1).tolist() 
                end = (np.array(temp[1::2]) -1).tolist()
                temp = [[i,j] for i,j in zip(start,end)]
                if Info[e][6] == '+':
                    for i,j in zip(temp,range(1,l)):
                        t = Info[e][0:2] + ["intron"] + [str(i[0])] + [str(i[1])] + Info[e][5:8] + [f'ID={m}.intron{j};Parent={m}\n']
                        gff_output.write('\t'.join(t))
                if Info[e][6] == '-':
                    for i,j in zip(temp,range(l-1,0,-1)):
                        t = Info[e][0:2] + ["intron"] + [str(i[0])] + [str(i[1])] + Info[e][5:8] + [f'ID={m}.intron{j};Parent={m}\n']
                        gff_output.write('\t'.join(t))                                
            gff_output.write('\n')
    gff_output.close()

parser = argparse.ArgumentParser(description='gffkit: a cross-platform and toolkit for gff file manipulation',usage='gffkit [command]',add_help=False,epilog='date:2022/04/05 author:guisen chen email:thecgs001@foxmail.com')
parser.add_argument('-h',action='help',help='Help')
parser.add_argument('-v',action='version',help='Version:1.00')
parent_parser = parser.add_subparsers()
seq_parser = parent_parser.add_parser('seq',help='Extract sequence',usage='gffkit seq -fastx genome.fa -gff file.gff -out out.fa -feature upstream -len 2000',add_help=False)
seq_parser.add_argument('-fastx',metavar='genome.fa',required=True,help='Open a plain or gzipped fasta file')
seq_parser.add_argument('-gff',metavar='file.gff',required=True,help = 'Open a gff file')
seq_parser.add_argument('-out',metavar='out.fa',required=True,help = 'Output a fasta file')
seq_parser.add_argument('-feature',metavar='pep',help = 'Only the "gene, mRNA, exon, intron, upstream, CDS, pep, gene-like, RNA-like" fields are supported',required=True)
seq_parser.add_argument('-len',metavar='num',type=int,default=2000,help='Use only when "-feature upstream", default 2000bp')
seq_parser.add_argument('-str',metavar='str',type=str,help='Use only when "-feature gene-like or RNA-link", such as "-str tRNA"')
seq_parser.add_argument('-h',action='help',help='Help')
seq_parser.set_defaults(func=seq)
stat_parser = parent_parser.add_parser('stat',help='stat gff file',usage='gffkit stat -gff file.gff',add_help=False)
stat_parser.add_argument('-gff',metavar='file.gff',required=True,help = 'Open a gff file')
stat_parser.add_argument('-h',action='help',help='Help')
stat_parser.set_defaults(func=stat)
add_parser = parent_parser.add_parser('add',help='Add an intron feature to the gff file',usage='gffkit add -gff file.gff -out out.gff',add_help=False)
add_parser.add_argument('-gff',metavar='file.gff',required=True,help = 'Open a gff file')
add_parser.add_argument('-out',metavar='out.gff',required=True,help = 'Output a gff file ')
add_parser.add_argument('-h',action='help',help='Help')
add_parser.set_defaults(func=add)
filter_parser = parent_parser.add_parser('filter',help='Get a the longest transcription',usage='gffkit filter -gff file.gff -out out.gff', add_help=False)
filter_parser.add_argument('-gff',metavar='file.gff',required=True,help = 'Open a gff file')
filter_parser.add_argument('-out',metavar='out.gff',required=True,help = 'Output a gff file')
filter_parser.add_argument('-h',action='help',help='Help')
filter_parser.set_defaults(func=add)
args = parser.parse_args()
args.func(args)


# In[26]:


import argparse
parser = argparse.ArgumentParser(description='gffkit: a cross-platform and toolkit for gff file manipulation',usage='gffkit [command]',add_help=False,epilog='date:2022/04/05 author:guisen chen email:thecgs001@foxmail.com')

parser.add_subparsers()

