import sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Data import CodonTable
table = CodonTable.ambiguous_dna_by_id[1]

def main():
    sequence_file = sys.argv[1]
    mutation_file = sys.argv[2]
    gene_file = sys.argv[3]
    output_file = sys.argv[4]
    output_dict = {}
    
    sequence_dict = read_sequence_file(sequence_file) 
    mutation_dict = get_mutation(mutation_file,sequence_dict,output_dict)
    gene_dict = get_genes(gene_file)
    data_dict = mutation_identification(mutation_dict,gene_dict)
    output_dict = get_Codons(data_dict)
    write_output_file(output_dict,output_file)

    
# creates dictionary with ID and sequence
def read_sequence_file(sequence_file):
    sequence_dict = {}
    
    fasta_sequences = SeqIO.parse(open(sequence_file),'fasta')
    for fasta in fasta_sequences:
        sequence_dict[str(fasta.id)] = str(fasta.seq)
            
    fasta_sequences.close
    return(sequence_dict)

def get_genes(genefile):
    gene_counter = 0
    gene_file = open(genefile, 'r')
    gene_dict = {}
    
    next(gene_file)
    for line in gene_file:
        line = line.strip()
        line = line.split("\t")
        gene_ID = line[0]        
        START = int(line[3])
        STOP = int(line[4])
        STRAND = line[6]
        DESCRIPTION = line[8]
        if line[2] =="CDS":
            gene_counter += 1
            gene_dict[gene_counter] = [gene_ID,START,STOP,STRAND,DESCRIPTION]

    return gene_dict

# matches the mutation with the proper sequence
def get_mutation(mutation_input,sequence_dict,output_dict):
    data_dict = {}
   
    for sequence_key in sequence_dict:
        mutation_count = 0
        mutation_file = open(mutation_input, "r")
        next(mutation_file)
        for line in mutation_file:
            line = line.strip()
            line = line.split("\t")
            ORG = line[2]
            ID = line[0]
            POS = line[1]
            ALT = line[3]
            if ID == sequence_key:
                    mutation_ID = sequence_key +"_"+ str(mutation_count)
                    sequence = sequence_dict[sequence_key]
                    data_dict[mutation_ID] = ID,sequence,POS,ORG,ALT
                    mutation_count += 1
        mutation_file.close()    
    return data_dict

def mutation_identification(mutation_dict,gene_dict):
    hypothetical = 0
    positive = 0
    negative = 0
    annotated_gene = 0
    final_dict = {}
    for key in mutation_dict:
        for key_2 in gene_dict:
            if gene_dict[key_2][0] == mutation_dict[key][0]:
                description = gene_dict[key_2][4]
                
                sequence = mutation_dict[key][1]
                pos = int(mutation_dict[key][2])
                start = int(gene_dict[key_2][1])
                stop = int(gene_dict[key_2][2])
                strand = gene_dict[key_2][3]
                org = mutation_dict[key][3]
                alt = mutation_dict[key][4]
                gen_range = range(start,stop+1)
                
                if len(org) == 1 and len(alt) > 1 and "," not in alt:
                    pass
            
                if len(org) > 1:
                    pass

                if len(org) == 1 and "," in alt:
                    if pos in gen_range and "hypothetical protein" in description:
                        hypothetical += 2
                    if pos in gen_range and "hypothetical protein" not in description:
                        annotated_gene +=2
                        gene_name = gene_dict[key_2][4].split(";")
                        gene_id = gene_name[0][3:]
                        if strand == "+":
                            positive +=2
                            final_dict[key] = sequence,pos,org,alt,gene_id
                        if strand == "-":
                            negative +=2
                            sequence = str(Seq(sequence).reverse_complement())
                            org = str(Seq(org).reverse_complement())
                            alterations = alt.split(",")
                            alt_1 = str(Seq(alterations[0]).reverse_complement())
                            alt_2 = str(Seq(alterations[1]).reverse_complement())
                            alt = alt_1 + "," + alt_2
                            final_dict[key] = sequence,pos,org,alt,gene_id


                if len(org) == 1 and len(alt) == 1 and "," not in alt: 
                    if pos in gen_range and "hypothetical protein" in description:
                            hypothetical += 1
                    if pos in gen_range and "hypothetical protein" not in description:
                        annotated_gene +=1
                        gene_name = gene_dict[key_2][4].split(";")
                        gene_id = gene_name[0][3:]
                        if strand == "+":
                            positive +=1
                            final_dict[key] = sequence,pos,org,alt,gene_id
                        if strand == "-":
                            negative +=1
                            sequence = str(Seq(sequence).reverse_complement())
                            org = str(Seq(org).reverse_complement())
                            alt = str(Seq(alt).reverse_complement())
                            final_dict[key] = sequence,pos,org,alt,gene_id

    non_coding = (len(mutation_dict)-hypothetical-annotated_gene)
    print(" mutation in \n - annotated genes " + str(annotated_gene) + "\n - non_coding regions " +
          str(non_coding) + "\n - hypothetical proteins " + str(hypothetical) + "\n - the positive strand " +
          str(positive) + "\n - the negative strand " + str(negative))
    return final_dict
                
# gets original codon and codon with the new mutation and translates both to amino acid
def get_Codons(data_dict):
    AA_dict = {}
    for key in data_dict:
        
        sequence = data_dict[key][0]
        mutation_POS = data_dict[key][1]
        ORG = data_dict[key][2]
        ALT = data_dict[key][3]
        POS = int(mutation_POS)-1
        GENE_ID = data_dict[key][4]
        
        #first positon of codon
        if POS in range (0,len(sequence),3):
            REF_codon = sequence[POS:(POS+3)]
            REF_RNA = Seq(REF_codon).transcribe()
            REF_AA = REF_RNA.translate()
            
            if len(ORG) == 1 and len(ALT) > 1 and "," not in ALT:
                pass
            
            if len(ORG) > 1:
                pass
                
            if len(ORG) == 1 and "," in ALT:
                count = 0
                alterations = ALT.split(",")
                for i in alterations:
                    new_key = key + "_" + str(count)
                    ALT_codon = i + REF_codon[1:]
                    ALT_RNA = Seq(ALT_codon).transcribe()
                    ALT_AA = ALT_RNA.translate()
                    if str(REF_AA) == str(ALT_AA):
                        AA_dict[new_key] = str(REF_AA),str(ALT_AA),"synonymous",GENE_ID
                    if str(REF_AA) != str(ALT_AA):
                        AA_dict[new_key] = str(REF_AA),str(ALT_AA),"non-synonymous",GENE_ID
                    count += 1
                    
            if len(ORG) == 1 and len(ALT) == 1 and "," not in ALT:      
                ALT_codon = REF_codon[1:3]+ ALT
                ALT_RNA = Seq(ALT_codon).transcribe()
                ALT_AA = ALT_RNA.translate()
                if str(REF_AA) == str(ALT_AA):
                    AA_dict[key] = str(REF_AA),str(ALT_AA),"synonymous",GENE_ID
                if str(REF_AA) != str(ALT_AA):
                    AA_dict[key] = str(REF_AA),str(ALT_AA),"non-synonymous",GENE_ID

        # second position of codon
        if POS in range (1,len(sequence),3):
            REF_codon = sequence[(POS-1):(POS+2)]
            REF_RNA = Seq(REF_codon).transcribe()
            REF_AA = REF_RNA.translate()
            
            if len(ORG) == 1 and len(ALT) > 1 and "," not in ALT:
                pass
            
            if len(ORG) > 1:
                pass
                
            if len(ORG) == 1 and "," in ALT:
                count = 0
                alterations = ALT.split(",")
                for i in alterations:
                    new_key = key + "_" + str(count)
                    ALT_codon = REF_codon[0]+ i + REF_codon[2]
                    ALT_RNA = Seq(ALT_codon).transcribe()
                    ALT_AA = ALT_RNA.translate()
                    if str(REF_AA) == str(ALT_AA):
                        AA_dict[new_key] = str(REF_AA),str(ALT_AA),"synonymous",GENE_ID
                    if str(REF_AA) != str(ALT_AA):
                        AA_dict[new_key] = str(REF_AA),str(ALT_AA),"non-synonymous",GENE_ID
                    count += 1            
            if len(ORG) == 1 and len(ALT) == 1 and "," not in ALT:       
                ALT_codon = REF_codon[1:3]+ ALT
                ALT_RNA = Seq(ALT_codon).transcribe()
                ALT_AA = ALT_RNA.translate()
                if str(REF_AA) == str(ALT_AA):
                    AA_dict[key] = str(REF_AA),str(ALT_AA),"synonymous",GENE_ID
                if str(REF_AA) != str(ALT_AA):
                    AA_dict[key] = str(REF_AA),str(ALT_AA),"non-synonymous",GENE_ID

        # third position of codon
        if POS in range (2,len(sequence),3):
            REF_codon = sequence[(POS-2):(POS+1)]
            REF_RNA = Seq(REF_codon).transcribe()
            REF_AA = REF_RNA.translate()

            if len(ORG) == 1 and len(ALT) > 1 and "," not in ALT:
                pass
            
            if len(ORG) > 1:
                pass
                
            if len(ORG) == 1 and "," in ALT:
                count = 0
                alterations = ALT.split(",")
                for i in alterations:
                    new_key = key + "_" + str(count)
                    ALT_codon = REF_codon[1:3]+ i
                    ALT_RNA = Seq(ALT_codon).transcribe()
                    ALT_AA = ALT_RNA.translate()
                    if str(REF_AA) == str(ALT_AA):
                        AA_dict[new_key] = str(REF_AA),str(ALT_AA),"synonymous",GENE_ID
                    if str(REF_AA) != str(ALT_AA):
                        AA_dict[new_key] = str(REF_AA),str(ALT_AA),"non-synonymous",GENE_ID
                    count += 1
            if len(ORG) == 1 and len(ALT) == 1 and "," not in ALT:       
                ALT_codon = REF_codon[1:3]+ ALT
                ALT_RNA = Seq(ALT_codon).transcribe()
                ALT_AA = ALT_RNA.translate()
                if str(REF_AA) == str(ALT_AA):
                    AA_dict[key] = str(REF_AA),str(ALT_AA),"synonymous",GENE_ID
                if str(REF_AA) != str(ALT_AA):
                    AA_dict[key] = str(REF_AA),str(ALT_AA),"non-synonymous",GENE_ID

    return AA_dict

def write_output_file(output_dict,output_file):
    output = open(output_file,"w+")
    output.write("ID" + "\t" + "ORG AA" + "\t" + " ALT AA" + "\t" + "synonymous/non-synonymous" + "\t" + "gene ID" + "\n" )
    for key in output_dict:
        output.write(key + "\t" + output_dict[key][0] + "\t" + output_dict[key][1]  + "\t" + output_dict[key][2] + "\t" + output_dict[key][3] + "\n")
    output.close()

main()
