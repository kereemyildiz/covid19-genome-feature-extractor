from configs import configs
import json
import random


def translate_nucleotides_to_amino_acids(nucleotide_sequence, codon_mapper):
    
    aa_seq = ""
    nucleotide_sequence = "".join(nucleotide_sequence)
    # Iterate over the nucleotide sequence in steps of 3 (size of a codon)
    for i in range(0, len(nucleotide_sequence), 3):
        codon = nucleotide_sequence[i:i+3]
        # Translate the codon to an amino acid if possible, and add it to the sequence
        if codon in codon_mapper:
            aa = codon_mapper[codon]
            aa_seq += aa

    return aa_seq

def find_protein_region(index, protein_regions):
    
    for protein, (start, end) in protein_regions.items():
        if start <= index <= end:
            return protein
        
    return None

def get_sequence():
    
    with open('genome.txt', 'r') as file:
        # Skip the first line
        next(file)
        
        # ignore new lines
        genome = file.read().replace('\n', '')
    
        # genome = [random.choice(['A', 'T', 'G', 'C']) for _ in range(30000)]
        
        genome += 'A' # in order to make the length divisible by 3 (temporarily) 
        
        return genome
    
    
def get_elapsed_day():
    
    return None

def get_depth():
    
    return None    
    

def main():
    nucleotides = ['A', 'T', 'G', 'C']
    k = 30 # window size
    mid_point = k // 2
    
    config_file = configs() 
    protein_regions = config_file['protein regions']
    
    codon_mapping_file = open('codon_aa_mapping.json')
    codon_mapper = json.load(codon_mapping_file)
    
    seq = get_sequence()
    
    features = []
    aa_seq = translate_nucleotides_to_amino_acids(seq[mid_point:-mid_point], codon_mapper)
    count1 = 0
    count2 = 0
    
    for idx in range(len(seq)-k):
        
        window = seq[idx:idx+k]
        aa_idx = idx//3
        aa_nucleotides = list(seq[mid_point:-mid_point][aa_idx*3:aa_idx*3+3])
        mutation_position = idx % 3
        
        sample_data = []
        
        for nucleotide in nucleotides:
            
            aa_nucleotides[mutation_position] = nucleotide     
            new_aa = codon_mapper[''.join(aa_nucleotides)]

            sample_data.append([
                *list(window), # 1...k => window
                window[mid_point], # k+1 => center of k-mer
                nucleotide, # k+2 => nucleotide after mutation
                idx, # k+3 => index
                config_file['nucleotide sub. matrix'][window[mid_point]][nucleotide], # k+4 PAM score
                aa_seq[aa_idx], # k+5 => aa before mutation
                new_aa, # k+6 => aa after mutation 
                config_file['AA PAM matrix'][aa_seq[aa_idx]][new_aa], # k+7 => PAM score at amino acid level 
    #             getElapsedDay(), # k+8 => elapsed day
    #             getDepth(), # k+9 => depth ?
                int(aa_seq[aa_idx] == new_aa), # k+10 => synonymous or non-syn
                find_protein_region(idx, protein_regions), # k+11 => protein region
                config_file['AA features']['hydrophobicity'][new_aa],
                config_file['AA features']['polarity'][aa_seq[aa_idx]],
                config_file['AA features']['polarity'][new_aa],
                config_file['AA features']['iso-electric point'][aa_seq[aa_idx]],
                config_file['AA features']['iso-electric point'][new_aa],
                config_file['AA features']['volume'][aa_seq[aa_idx]],
                config_file['AA features']['volume'][new_aa],
                config_file['AA features']['molecular weight'][aa_seq[aa_idx]],
                config_file['AA features']['molecular weight'][new_aa],
                config_file['AA features']['pKa'][aa_seq[aa_idx]],
                config_file['AA features']['pKa'][new_aa],
                config_file['AA features']['pKb'][aa_seq[aa_idx]],
                config_file['AA features']['pKb'][new_aa],
                config_file['AA features']['pKx'][aa_seq[aa_idx]],
                config_file['AA features']['pKx'][new_aa],
                config_file['AA features']['pl'][aa_seq[aa_idx]],
                config_file['AA features']['pl'][new_aa]
                
            ])
        features += sample_data
        
    # print(features)
            
        
main()