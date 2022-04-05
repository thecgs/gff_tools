#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import re
import pyfastx
import argparse
import numpy as np
from collections import defaultdict
parser = argparse.ArgumentParser(description='Extract the sequence from the GFF file',usage='python from_gff_Extract_seq -genome genome.fa -gff file.gff -feature gene/mRNA/exon/intron/promoter/CDS/pep -out out.fa',add_help=False,epilog='date:2022/04/05 author:guisen chen email:thecgs001@foxmail.com')
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument('-genome','--genome',metavar='[genome.fa]',help='genome.fa',required=True)
required.add_argument('-gff','--gff',metavar='[file.gff]',help='file.gff',required=True)
required.add_argument('-feature','--feature',metavar='[gene/mRNA/exon/intron/promoter/CDS/pep]',help='gene/mRNA/exon/intron/promoter/CDS/pep',required=True)
required.add_argument('-out','--output',metavar='[out.fa]',help='out.fa',required=True)
optional.add_argument('-h','--help',action='help',help='show this help message and exit')
optional.add_argument('-v','--version',action='version',version='v1.00')
args = parser.parse_args()

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

fasta = open(args.output,'w+')
gff_input = open(args.gff,'r')
genome_input = pyfastx.Fastx(args.genome)

#Genomes are too big and consume too much time and memory
#with open(args.genome,'r') as f:
#    genome_fa = defaultdict(str)
#    for line in f:
#        if line.startswith('#') or line == "\n": #remove '#' and '\n'
#            continue
#        if line.startswith('>'):
#            chr_ID = line.replace('>','').split()[0]    
#        else:
#            genome_fa[chr_ID] += line.replace('\n','').strip()
#    #print(genome_fa)
##genome_input.close()

genome_fa = defaultdict(str)
for name,seq,comment in genome_input:
    genome_fa[name] = seq

#parser gff file
Info = defaultdict(list)
gene2mRNA = defaultdict(list)
exon = defaultdict(list)
intron = defaultdict(list)
CDS = defaultdict(list)
for line in gff_input:
    if line.startswith('#') or line == "\n": #remove '#' and '\n'
        continue        
    line_list = line.split("\t")
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
#print(Info)
#print(gene2mRNA)
#print(CDS)
#print(exon)

if args.feature == 'gene':
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
    #fasta.close()

if args.feature == 'mRNA':
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
    #fasta.close()

if args.feature == 'exon':
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
    #fasta.close()

if args.feature == 'intron':
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
    #fasta.close()

if args.feature == 'promoter':
    for gene_ID,mRNA_ID in gene2mRNA.items():
        for m in mRNA_ID:
            fasta.write(f">{m}.promoter\n")
            chr_ID = Info[m][0]
            if Info[m][6] == '-':
                pos = int(Info[gene_ID][3:5][1])
                if len(genome_fa[chr_ID]) < pos+2001:
                    seq = genome_fa[chr_ID][pos+1:]
                    seq = seq_compliment(seq)[::-1]
                    fasta.write(seq+"\n")
                if len(genome_fa[chr_ID]) >= pos+2001:
                    seq = genome_fa[chr_ID][pos+1:pos+2001]
                    seq = seq_compliment(seq)[::-1]
                    fasta.write(seq+"\n")              
            if Info[m][6] == '+':
                pos = int(Info[gene_ID][3:5][0])
                if 0 <= pos-2001:
                    seq = genome_fa[chr_ID][pos-2001:pos-1]
                    fasta.write(seq+"\n")
                if 0 > pos-2001:
                    seq = genome_fa[chr_ID][:pos-1]
                    fasta.write(seq+"\n")                    
    #fasta.close()

#print(CDS)
if args.feature == 'pep' or 'CDS':
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
        fa = pyfastx.Fastx(args.output)
        os.remove(args.output)
        fasta = open(args.output,'a')
        for name, seq, comment in fa:
            fasta.write(f">{name}\n{seq_translate(seq)}\n")
        fasta.close()

print('finished')

