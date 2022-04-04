#!/usr/bin/env python
# coding: utf-8

# In[247]:


import re
import argparse
import numpy as np
from collections import defaultdict

parser = argparse.ArgumentParser(description='Add an intron feature to the GFF file',usage='python3 add_intron_to_gff.py -i [input] -o [output]',add_help=False,epilog='date:2022/04/04 author:guisen chen email:thecgs001@foxmail.com')
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument('-i','--input',metavar='[input_file]',help='input_file',required=True)
required.add_argument('-o','--output',metavar='[output_file]',help='output_file',required=True)
optional.add_argument('-h','--help',action='help',help='show this help message and exit')
optional.add_argument('-v','--version',action='version',version='v1.00')
args = parser.parse_args()

gff_input = open(args.input,'r')
gff_output = open(args.output,'w')

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


# In[ ]:




