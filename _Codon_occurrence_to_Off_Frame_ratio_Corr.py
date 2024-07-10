#The correlation coefficient between 
#codon content and Off frame ratio
###python _Codon_occurrence_to_Off_Frame_ratio_Corr.py -r reference -i off-frame_data -o outputfile
import pandas as pd
import csv
import re
import argparse
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
def process_data(input_longest_file,fsr_file):
    seq_dict = {}
    with open(input_longest_file, 'r') as f:
        lines = f.readlines()
        for i in range(len(lines)):
            if lines[i][0] == '>':
                gene_name = lines[i][1:].strip().split("_")[0]#.split("-")[1]
                seq_dict[gene_name] = ''
            else:
                seq_dict[gene_name] += lines[i].strip()
    dic_cds={}
    dic_seqlen={}
    for gene_name, seq in seq_dict.items():
        seq=seq[100:-100].upper()
        dic_cds[gene_name]=seq
        dic_seqlen[gene_name]=len(seq)
    def cut_text(text,lenth): 
        textArr = re.findall('.{'+str(lenth)+'}', text) 
        return textArr

    codon_seq = {}
    for k, v in dic_cds.items():
        x = str(k)
        codon=cut_text(v,3)
        codon_seq[x]=codon
    dic_frs={}
    with open(fsr_file, 'r') as file:
        for i, line in enumerate(file):
            a=line.strip().split("\t")
            name=a[0]#.split("-")[1]#.split(":")[1]#.upper()
            frs=float(a[1])
            dic_frs[name]=frs
    header = ["NAME","lenth","data"]
    for co in codontable.keys():
        header.append(co)
    datas =[]
    for k1,v1 in codon_seq.items():
            se = pd.Series(v1)
            if k1 in dic_frs.keys():
                frq = dict(se.value_counts(normalize=True))
                frq.update({"NAME": k1}) 
                frq.update({"lenth": dic_seqlen[k1]}) 
                frq.update({"data": dic_frs[k1]}) 
                datas.append(frq)
    df = pd.DataFrame(datas)
    return df

def calculate_correlation(input_file, output_file):
    df = input_file
    df = df.drop(['TGA', "TAA", "TAG"], axis=1)
    df_corr = df.corr(numeric_only=True)  
    df_corr.to_csv(output_file, index=False)
    return df_corr

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some files.')
    parser.add_argument('-r', type=str, required=True, help='The path to the input longest file')
    parser.add_argument('-i', type=str, required=True, help='The path to the fsr file')
    parser.add_argument('-o', type=str, required=True, help='The path to the output file')
    args = parser.parse_args()

    input1=process_data(args.r, args.i)
    input2=calculate_correlation(input1, args.o)
    name=args.i[:-4]
    df = input2
    result = dict(zip(df.columns, df['data']))
    result.pop("data")
    result.pop("lenth")

    import matplotlib.pyplot as plt
    plt.figure(figsize=(14, 7))
    sorted_data = dict(sorted(result.items(), key=lambda item: item[1], reverse=True))
    sorted_keys = list(sorted_data.keys())
    sorted_values = list(sorted_data.values())
    color_list = ['green' if result[key]>0 else 'red' for key in sorted_keys]
    plt.bar(sorted_keys, sorted_values, color=color_list)
    plt.xticks(rotation=90)
    plt.ylim(-0.4, 0.4)
    plt.ylabel('CSC')
    plt.title('Bar Chart with '+name)

    plt.savefig(name+".pdf")
    plt.show()
