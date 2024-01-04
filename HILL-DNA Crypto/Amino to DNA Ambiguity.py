#1ST part of the code is a theory we would use in order to eliminate ambiguity issue 
def encode_sequence(seq):
    encoding = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}
    encoded_values = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        binary_codon = ''.join(encoding[nucleotide] for nucleotide in codon)
        encoded_values.append(int(binary_codon, 2))
    return encoded_values

def decode_sequence(encoded_values):
    decoding = {'00': 'A', '01': 'C', '10': 'G', '11': 'T'}
    decoded_seq = ''
    for value in encoded_values:
        binary_codon = format(value, '06b')  # Ensure it's 6 bits long
        for i in range(0, len(binary_codon), 2):
            decoded_seq += decoding[binary_codon[i:i+2]]
    return decoded_seq

# Example 
sequence = "CGCGCGTTCGTT"
encoded_values = encode_sequence(sequence)
print(f"Encoded values: {encoded_values}")
decoded_sequence = decode_sequence(encoded_values)
print(f"Decoded sequence: {decoded_sequence}")


#2nd PART of the code shows how many possibilities we could get from one AMINO SEQUENCE.
def amino_to_dna(amino_acid):
    codon_table = {
        'A': ['GCT', 'GCC', 'GCA', 'GCG'],
        'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'N': ['AAT', 'AAC'],
        'D': ['GAT', 'GAC'],
        'C': ['TGT', 'TGC'],
        'Q': ['CAA', 'CAG'],
        'E': ['GAA', 'GAG'],
        'G': ['GGT', 'GGC', 'GGA', 'GGG'],
        'H': ['CAT', 'CAC'],
        'I': ['ATT', 'ATC', 'ATA'],
        'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
        'K': ['AAA', 'AAG'],
        'M': ['ATG'],
        'F': ['TTT', 'TTC'],
        'P': ['CCT', 'CCC', 'CCA', 'CCG'],
        'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
        'T': ['ACT', 'ACC', 'ACA', 'ACG'],
        'W': ['TGG'],
        'Y': ['TAT', 'TAC'],
        'V': ['GTT', 'GTC', 'GTA', 'GTG'],
        '*': ['TAA', 'TAG', 'TGA']  # Stop codons
    }
    return codon_table.get(amino_acid, [])

def generate_dna_sequences(amino_acid_sequence):
    if not amino_acid_sequence:
        return ['']
    first_amino_acid = amino_acid_sequence[0]
    rest_amino_acids = amino_acid_sequence[1:]
    possible_dna_sequences = amino_to_dna(first_amino_acid)
    remaining_dna_sequences = generate_dna_sequences(rest_amino_acids)
    return [dna + remainder for dna in possible_dna_sequences for remainder in remaining_dna_sequences]

def main():
    amino_acid_sequence = input("Enter the amino acid sequence: ").upper()
    dna_sequences = generate_dna_sequences(amino_acid_sequence)
    print(f"All possible DNA sequences that encode the amino acid sequence {amino_acid_sequence}:")
    for dna in dna_sequences:
        print(dna)

if __name__ == "__main__":
    main()
