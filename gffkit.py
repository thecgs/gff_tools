#!/usr/bin/env python
# coding: utf-8

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
        if args.feature == 'parent':
            if line_list[2] == args.str:
            #such as "region","pseudogene",'cDNA_match'
                ID = re.search("ID=(.*?)[;,\n]",line_list[8]).group(1).strip()
                Info[ID] = line_list
                gene_like_list.append(ID)
        if args.feature == "child":
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
            if line_list[2] == 'gene' or line_list[2] == 'pseudogene':
                gene_ID = re.search("ID=(.*?)[;,\n]",line_list[8]).group(1).strip()
                Info[gene_ID] = line_list       
            if line_list[2] == 'mRNA' or line_list[2] == 'transcript':
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
    if args.feature == 'parent':
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
    if args.feature == 'child':
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
        length = int(args.len)
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
    gene_length_dict = defaultdict(list)
    mRNA_length_dict = defaultdict(list)
    CDS_length_dict = defaultdict(list)
    exon_length_dict = defaultdict(list)
    parent_set = set()
    feature_list = defaultdict()
    with open(args.gff,'r') as f:
        for line in f:
            if line.startswith('#') or line == "\n": #remove '#' and '\n'
                continue
            line = line.rstrip().split("\t")
            ID = re.search("ID=(.*?)[;,\n]",line[8]).group(1).strip()
            feature = line[2]
            if feature not in feature_list:
                feature_list[feature] = [ID]
            if feature in feature_list:
                feature_list[feature].append(ID)
            if re.search("Parent=(.*?)[;,\n]",line[8]) == None:
                pass
            else:
                parent = re.search("Parent=(.*?)[;,\n]",line[8]).group(1).strip()
                parent_set.add(parent)
            if line[2] == 'gene' or line[2] == 'pseudogene':
                gene_id = re.search("ID=(.*?)[;,\n]",line[8]).group(1).strip()
                gene_length = int(line[4]) - int(line[3]) + 1
                gene_length_dict[gene_id].append(gene_length)             
            if line[2] == "mRNA" or line[2] == "transcript":
                mRNA_id = re.search("ID=(.*?)[;,\n]",line[8]).group(1).strip()
                #gene_id = re.search("Parent=(.*?)[;,\n]",line[8]).group(1).strip()
                mRNA_length = int(line[4]) - int(line[3]) + 1
                mRNA_length_dict[mRNA_id].append(mRNA_length)
            if line[2] == "CDS":
                mRNA_id = re.search("Parent=(.*?);",line[8]+";").group(1)
                CDS_length = int(line[4]) - int(line[3]) + 1
                CDS_length_dict[mRNA_id].append(CDS_length)
            if line[2] == "exon":
                mRNA_id = re.search("Parent=(.*?);",line[8]+";").group(1)
                exon_length = int(line[4]) - int(line[3]) + 1
                exon_length_dict[mRNA_id].append(exon_length)
            

    def gene_number():
        gene_number = len(gene_length_dict)
        return gene_number
    
    def mRNA_number():
        mRNA_number = len(mRNA_length_dict)
        return mRNA_number
    
    def gene_average_length():
        """
        sum(gene_length)/gene_number
        """
        temp = []
        for gene_id,gene_length in gene_length_dict.items():
            temp.append(gene_length[0])
        gene_average_length = sum(temp)/len(gene_length_dict)
        return gene_average_length

    def mRNA_average_length():
        """
        sum(mRNA_length)/mRNA_number
        """
        temp = []
        for mRNA_id,mRNA_length in mRNA_length_dict.items():
            temp.append(mRNA_length[0])
        mRNA_average_length = sum(temp)/mRNA_number()
        return mRNA_average_length
    
    def CDS_average_length():
        """
        sum(CDS_length)/CDS_number
        """
        temp = []
        number = []
        for mRNA_id,CDS_length in CDS_length_dict.items():
            temp.append(sum(CDS_length))
            number.append(len(CDS_length))
        CDS_average_length = sum(temp)/sum(number)
        return CDS_average_length
    
    def exon_average_length():
        """
        sum(exon_length)/exon_number
        """
        temp = []
        number = []
        for mRNA_id,exon_length in exon_length_dict.items():
            if mRNA_id in mRNA_length_dict:
                temp.append(sum(exon_length))
                number.append(len(exon_length))
            if mRNA_id not in mRNA_length_dict:
                pass
        exon_average_length = sum(temp)/sum(number)
        return exon_average_length
    
    def intron_average_length():
        """
        intron_length = mRNA_length - exon_length
        intron_average_length = intron_length/mRNA_number
        """
        temp = []
        number = []
        for mRNA_id,exon_length in exon_length_dict.items():
            if mRNA_id in mRNA_length_dict:
                intron_length = int(mRNA_length_dict[mRNA_id][0]) - int(sum(exon_length))
                temp.append(intron_length)
                number.append(len(exon_length)-1)
            if mRNA_id not in mRNA_length_dict:
                pass
        intron_average_length = sum(temp)/sum(number)
        return intron_average_length
    
    def per_mRNA_CDS_average_length():
        """
        sum(per_mRNA_all_CDS_length)/(per_mRNA_CDS_number*mrna_number)
        """
        temp = []
        for mrna_id,CDS_length_list in CDS_length_dict.items():
            temp.append(sum(CDS_length_list)/len(CDS_length_list))
        per_mRNA_CDS_average_length = sum(temp)/mRNA_number()
        return per_mRNA_CDS_average_length
    
    def per_mRNA_exon_average_length():
        """
        sum(per_mRNA_all_exon_length)/(per_mRNA_exon_number*mRNA_number)
        """
        temp = []
        for mRNA_id,exon_length in exon_length_dict.items():
            if mRNA_id in mRNA_length_dict:
                temp.append(sum(exon_length)/len(exon_length))
            if mRNA_id not in mRNA_length_dict:
                pass
        per_mRNA_exon_average_length = sum(temp)/mRNA_number()
        return per_mRNA_exon_average_length
    
    def per_mRNA_intron_average_length():
        """
        intron_length = mRNA_length - exon_length
        sum(per_mRNA_all_intron_length)/(per_mRNA_intron_number*mRNA_number)
        """
        temp = []
        for mRNA_id,exon_length in exon_length_dict.items():
            if mRNA_id in mRNA_length_dict:
                intron_length = int(mRNA_length_dict[mRNA_id][0]) - int(sum(exon_length))
                if len(exon_length)-1 == 0:
                    pass
                if len(exon_length)-1 != 0:
                    temp.append(intron_length/(len(exon_length)-1))
            if mRNA_id not in mRNA_length_dict:
                pass
        per_mRNA_intron_average_length = sum(temp)/mRNA_number()
        return per_mRNA_intron_average_length 
       
    def per_mRNA_all_CDS_average_length():
        """
        sum(per_mRNA_all_CDS_length)/mRNA_number
        """
        temp = []
        for mRNA_id,CDS_length in CDS_length_dict.items():
            temp.append(sum(CDS_length))
        per_mRNA_all_CDS_average_length = sum(temp)/mRNA_number()
        return per_mRNA_all_CDS_average_length
    
    def per_mRNA_all_exon_average_length():
        """
        sum(per_mRNA_all_exon_length)/mRNA_number
        """
        temp = []
        for mRNA_id,exon_length in exon_length_dict.items():
            if mRNA_id in mRNA_length_dict:
                temp.append(sum(exon_length))
            if mRNA_id not in mRNA_length_dict:
                pass
        per_mRNA_all_exon_average_length = sum(temp)/mRNA_number()
        return per_mRNA_all_exon_average_length

    def per_mRNA_all_intron_average_length():
        """
        sum(per_mRNA_all_intron_length)/mRNA_number
        """
        temp = []
        for mRNA_id,exon_length in exon_length_dict.items():
            if mRNA_id not in mRNA_length_dict:
                pass
            if mRNA_id in mRNA_length_dict:
                intron_length = int(mRNA_length_dict[mRNA_id][0]) - int(sum(exon_length))            
                temp.append(intron_length)
        per_mRNA_all_intron_average_length = sum(temp)/mRNA_number()
        return per_mRNA_all_intron_average_length
    
    def per_mRNA_CDS_average_number():
        """
        sum(per_mRNA_CDS_number)/mRNA_number
        """
        temp = []
        for mRNA_id,CDS_length in CDS_length_dict.items():
            temp.append(len(CDS_length))
        per_mRNA_CDS_average_number = sum(temp)/mRNA_number()
        return per_mRNA_CDS_average_number

    def per_mRNA_exon_average_number():
        """
        sum(per_mRNA_exon_number)/mRNA_number
        """
        temp = []
        for mRNA_id,exon_length in exon_length_dict.items():
            temp.append(len(exon_length))
        per_gene_average_exons_number = sum(temp)/len(mRNA_length_dict)
        return per_gene_average_exons_number
    
    def per_mRNA_intron_average_number():
        """
        sum(per_mRNA_intron_number)/mRNA_number
        """
        temp = []
        for mRNA_id,exon_length in exon_length_dict.items():            
            temp.append(len(exon_length)-1)
        per_gene_average_intron_number = sum(temp)/len(mRNA_length_dict)
        return per_gene_average_intron_number
    print("\n\
-------------------------------------------------------------------------------------------------------\n\
A formula to calculate:\n\
-------------------------------------------------------------------------------------------------------\n\
gene_average_length = sum(gene_length)/gene_number\n\
mRNA_average_length = sum(mRNA_length)/mRNA_number\n\
CDS_average_length= sum(CDS_length)/CDS_number\n\
exon_average_length= sum(exon_length)/exon_number\n\
intron_length = mRNA_length - exon_length\n\
intron_average_length= sum(intron_length)/intron_number\n\
per_mRNA_CDS_average_length = sum(per_mRNA_all_CDS_length)/(per_mRNA_CDS_number*mRNA_number)\n\
per_mRNA_exon_average_length = sum(per_mRNA_all_exon_length)/(per_mRNA_exon_number*mRNA_number)\n\
per_mRNA_intron_average_length = sum(per_mRNA_all_intron_length)/(per_mRNA_intron_number*mRNA_number)\n\
per_mRNA_all_CDS_average_length = sum(per_mRNA_all_CDS_length)/mRNA_number\n\
per_mRNA_all_exon_average_length = sum(per_mRNA_all_exon_length)/mRNA_number\n\
per_mRNA_all_intron_average_length = sum(per_mRNA_all_intron_length)/mRNA_number\n\
per_mRNA_CDS_average_number = sum(per_mRNA_CDS_number)/mRNA_number\n\
per_mRNA_exon_average_number = sum(per_mRNA_exon_number)/mRNA_number\n\
per_mRNA_intron_average_number = sum(per_mRNA_intron_number)/mRNA_number\n\
-------------------------------------------------------------------------------------------------------\n\
note: 1. XXX_length is value of [(col5 - col4) + 1] in row of XXX\n\
      2. gene(gene and pseudogene) \n\
      3. mRNA(mRNA and transcript) \n\
-------------------------------------------------------------------------------------------------------\n\n")
    print("Table 1")
    print("----------------------------------------------------")
    print("{:<17} {:<10} {:<0}".format("feature","number","Whether there are child",)) 
    print("----------------------------------------------------")
    for f,i in feature_list.items():  
        if i[0] in parent_set:
        #print(f"{f}:have child")        
            print("{:<17} {:<10} {:<0}".format(f,len(feature_list[f])-1,'have child'))
        if i[0] not in parent_set:
        #print(f"{f}:no clild")
            print("{:<17} {:<10} {:<0}".format(f,len(feature_list[f])-1,'no child'))
    print("----------------------------------------------------\n")       
    print("Table 2")
    print("----------------------------------------------------")
    print("{:<40} {:<0}".format("value","number"))
    print("----------------------------------------------------")
    print("{:<40} {:<0}".format("gene_number(gene and pseudogene)",gene_number()))
    print("{:<40} {:<0}".format("mRNA_number(mRNA and transcript)",mRNA_number()))
    print("{:<40} {:<0}".format("gene_average_length",round(gene_average_length(),2)))
    print("{:<40} {:<0}".format("mRNA_average_length",round(mRNA_average_length(),2)))
    print("{:<40} {:<0}".format("CDS_average_length",round(CDS_average_length(),2)))
    print("{:<40} {:<0}".format("exon_average_length",round(exon_average_length(),2)))
    print("{:<40} {:<0}".format("intron_average_length",round(intron_average_length(),2)))
    print("{:<40} {:<0}".format("per_mRNA_CDS_average_length",round(per_mRNA_CDS_average_length(),2)))
    print("{:<40} {:<0}".format("per_mRNA_exon_average_length",round(per_mRNA_exon_average_length(),2)))
    print("{:<40} {:<0}".format("per_mRNA_intron_average_length",round(per_mRNA_intron_average_length(),2)))
    print("{:<40} {:<0}".format("per_mRNA_all_CDS_average_length",round(per_mRNA_all_CDS_average_length(),2)))
    print("{:<40} {:<0}".format("per_mRNA_all_exon_average_length",round(per_mRNA_all_exon_average_length(),2)))
    print("{:<40} {:<0}".format("per_mRNA_all_intron_average_length",round(per_mRNA_all_intron_average_length(),2)))
    print("{:<40} {:<0}".format("per_mRNA_CDS_average_number",round(per_mRNA_CDS_average_number(),2)))
    print("{:<40} {:<0}".format("per_mRNA_exon_average_number",round(per_mRNA_exon_average_number(),2)))
    print("{:<40} {:<0}".format("per_mRNA_intron_average_number",round(per_mRNA_intron_average_number(),2)))
    print("----------------------------------------------------")
    
def filter(args):
    def parse_gff(file_gff):
        gene_ID_dict = defaultdict(list) # all mRNA id corresponding to gene id 
        ID_info_dict = defaultdict(list) # info of gene, mRNA, CDS and exon
        exon_length_dict = defaultdict(list) # all CDS length corresponding to mRNA id 
        the_longest_mRNA_ID_dict = defaultdict(str)
        for line in open(file_gff,'r'):
            if line.startswith('#') or line == "\n": #remove '#' and '\n'
                continue
            line_list = line.split("\t")
            if line_list[2] == "gene" or line_list[2] == "pseudogene": 
                gene_ID = re.search("ID=(.*?)[;,\n]",line_list[8]).group(1).strip()
                ID_info_dict[gene_ID] = list(line)
            if line_list[2] == "mRNA" or line_list[2] == "transcript": 
                mRNA_ID = re.search("ID=(.*?)[;,\n]",line_list[8]).group(1)
                gene_ID = re.search("Parent=(.*?)[;,\n]",line_list[8]).group(1).strip()
                if gene_ID in gene_ID_dict:
                    gene_ID_dict[gene_ID].append(mRNA_ID)
                else:
                    gene_ID_dict[gene_ID] = [mRNA_ID]
                ID_info_dict[mRNA_ID].append(list(line))
            if line_list[2] == "exon": 
                exon_ID = re.search("ID=(.*?)[;,\n]",line_list[8]).group(1)
                mRNA_ID = re.search("Parent=(.*?)[;,\n]",line_list[8]).group(1).strip()
                ID_info_dict[mRNA_ID].append(list(line))
                exon_length = int(line_list[4]) - int(line_list[3]) + 1
                exon_length_dict[mRNA_ID].append(exon_length)
            if line_list[2] == "CDS" or line_list[2] == "intron" or  line_list[2] == "polyA_sequence" or line_list[2] == "polyA_site" or line_list[2] == "three_prime_UTR" or line_list[2] == "five_prime_UTR":
                #exon_ID = re.search("ID=(.*?)[;,\n]",line_list[8]).group(1)
                mRNA_ID = re.search("Parent=(.*?)[;,\n]",line_list[8]).group(1).strip()
                ID_info_dict[mRNA_ID].append(list(line))
                
        for gene_ID,mRNA_ID_list in gene_ID_dict.items():
            i = 0
            for mRNA_ID in mRNA_ID_list:
                exon_length = sum(exon_length_dict[mRNA_ID])
                if exon_length > i:
                    i = exon_length
                    the_longest_mRNA_ID_dict[gene_ID] = mRNA_ID

        gff_file = open(args.out,'w')
        #gff = csv.writer(gff_file,delimiter='\t')
        for gene_ID, mRNA_ID in the_longest_mRNA_ID_dict.items():
            for line in ([ID_info_dict[gene_ID]]+ID_info_dict[mRNA_ID]+['\n']):
                for l in line:
                    gff_file.write(l)
        gff_file.close()
        return print('finished')
    parse_gff(args.gff)
def add(args):
    gff_input = open(args.gff,'r')
    gff_output = open(args.out,'w')
    Info = defaultdict(list)
    gene2mRNA = defaultdict(list)
    CDS = defaultdict(list)
    other = defaultdict(list)
    exon = defaultdict(list)

    for line in gff_input:
        if line.startswith('#') or line == "\n": #remove '#' and '\n'
            continue

        line_list = line.split("\t")
        #print(line_list)
        if line_list[2] == 'gene' or line_list[2] == 'pseudogene':
            gene_ID = re.search("ID=(.*?)[;,\n]",line_list[8]).group(1).strip()
            Info[gene_ID] = line_list

        if line_list[2] == 'mRNA' or line_list[2] == 'transcript':
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
                
        else:
            if re.search("Parent=(.*?)[;,\n]",line_list[8]) != None:
                mRNA_ID = re.search("Parent=(.*?)[;,\n]",line_list[8]).group(1).strip()
                if mRNA_ID in other:
                    other[mRNA_ID].append(line_list)
                else:
                    other[mRNA_ID] = [line_list]
            if re.search("Parent=(.*?)[;,\n]",line_list[8]) == None:
                pass

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
            for o in other[m]:
                gff_output.write("\t".join(o))
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
    print('finished')

parser = argparse.ArgumentParser(description='gffkit: a cross-platform and toolkit for gff file manipulation',usage='gffkit [command]',add_help=False,epilog='date:2022/04/05 author:guisen chen email:thecgs001@foxmail.com')
parser.add_argument('-h',action='help',help='Help')
parser.add_argument('-v',action='version',help='Version:1.00')
parent_parser = parser.add_subparsers()
seq_parser = parent_parser.add_parser('seq',help='Extract sequence',usage='gffkit seq -fastx genome.fa -gff file.gff -out out.fa -feature upstream -len 2000',add_help=False)
seq_parser.add_argument('-fastx',metavar='genome.fa',required=True,help='Open a plain or gzipped fasta file')
seq_parser.add_argument('-gff',metavar='file.gff',required=True,help = 'Open a gff file')
seq_parser.add_argument('-out',metavar='out.fa',required=True,help = 'Output a fasta file')
seq_parser.add_argument('-feature',metavar='pep',help = 'Only the "gene, mRNA, exon, intron, upstream, CDS, pep, parent, child" fields are supported',required=True)
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
filter_parser.set_defaults(func=filter)
args = parser.parse_args()
args.func(args)
