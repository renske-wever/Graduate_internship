from Bio.Seq import Seq
from Bio import SeqIO
import sys


sequence_file = sys.argv[1] #~/final_selection_bins/MESP_7-16_bin_2/prodigal/proteins_bin_2.fa
mutation_file = sys.argv[2] #dn-ds_prodigal_7_16_bin_2.txt
output_file = sys.argv[3] #dn-ds_prodigal_genes_7_16_bin_2.txt
out_fasta = sys.argv[4]
genes = {}
gene_sequences = {}
mutation_file = open(mutation_file, "r")
next(mutation_file)

for line in mutation_file:
    line = line.strip()
    line = line.split("\t")
    gene_id = line[5]#[4] with gff from rast
    dn_ds = line[3]
   
    if gene_id in genes.keys():
        if dn_ds == "non-synonymous":
            genes[gene_id][0]+=1
        if dn_ds == "synonymous":
            genes[gene_id][1]+=1
    else:
        if dn_ds == "non-synonymous":
            genes[gene_id] = [1,0]
        if dn_ds == "synonymous":
            genes[gene_id] = [0,1]
            
fasta_sequences = SeqIO.parse(open(sequence_file),'fasta')
for fasta in fasta_sequences:
    name, sequence = str(fasta.id), str(fasta.seq)
    for gene_id in genes:
        if name == gene_id:
            gene_sequences[name]=sequence

out_file = open(out_fasta,"w+")
for key in gene_sequences:
    out_file.write(">" + key + "\n" + gene_sequences[key] + "\n")
out_file.close()

output = open(output_file,"w+")
output.write("ID" + "\t" + "synonymous"+ "\t" +"non-synonymous" + "\n" )
for key in genes:
    output.write(key + "\t" + str(genes[key][0]) + "\t" + str(genes[key][1])  + "\n")
output.close()
