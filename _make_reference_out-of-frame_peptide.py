import re
import random
codontable = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'
    }
def cut_text(text,lenth): 
    textArr = re.findall('.{'+str(lenth)+'}', text) 
    return textArr

fw = open("reference_out-of-frame_peptide.fa","a")  ### output search protein sequence file (0,+1,-1) frame protein
fr = open("CDS_DNA_human.fa","r")                   ### input CDS_DNA_human.fa (need search protein) 

dic_cds={}
for line in fr:
    if ">" in line:
        name=line.strip().split(">")[1]
    else:
        seq=line.strip()[100:] ### remove 5'utr keep 3'utr 100 bp
        dic_cds[name]=seq

from collections import defaultdict
dic=defaultdict(list)

for k,v in dic_cds.items():
    name=k
    seq=v
    plus=seq[1:]
    minus=seq[2:]
    seqence=cut_text(seq,3)
    seqence1=cut_text(plus,3)
    seqence2=cut_text(minus,3)
    new_list = [codontable[a] for a in seqence]
    aaseq="".join(new_list)
    new_list1 = [codontable[a] for a in seqence1]
    aaseq1="".join(new_list1)
    new_list2 = [codontable[a] for a in seqence2]
    aaseq2="".join(new_list2)
    fw.write(">"+name+":"+"normal"+"\n"+aaseq+"\n"+">"+name+":"+"plus"+"\n"+aaseq1+"\n"+">"+name+":"+"minus"+"\n"+aaseq2+"\n")