import numpy as np
from math import gcd

# Convert text to its 8-bit binary representation
def text_to_8bit_binary(text):
    return ''.join(format(ord(c), '08b') for c in text)

# Convert binary representation to DNA sequence
def binary_to_adn(binary):
    mapping = {'00': 'A', '01': 'C', '10': 'G', '11': 'T'}
    while len(binary) % 2 != 0:
        binary += '0'
    return ''.join(mapping[binary[i:i+2]] for i in range(0, len(binary), 2))

# Convert DNA sequence to amino acid sequence using the codon table
def adn_to_amino_acid(adn_sequence):
    codon_table = { 
'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'AGC': 'S', 'AGT': 'S',
    # Phenylalanine (F)
    'TTC': 'F', 'TTT': 'F',
    # Leucine (L)
    'TTA': 'L', 'TTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    # Tyrosine (Y)
    'TAC': 'Y', 'TAT': 'Y',
    # Cysteine (C)
    'TGC': 'C', 'TGT': 'C',
    # Tryptophan (W)
    'TGG': 'W',
    # Proline (P)
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    # Histidine (H)
    'CAC': 'H', 'CAT': 'H',
    # Glutamine (Q)
    'CAA': 'Q', 'CAG': 'Q',
    # Arginine (R)
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'AGA': 'R', 'AGG': 'R',
    # Isoleucine (I)
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I',
    # Methionine (M)
    'ATG': 'M',
    # Threonine (T)
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    # Asparagine (N)
    'AAC': 'N', 'AAT': 'N',
    # Lysine (K)
    'AAA': 'K', 'AAG': 'K',
    # Serine (S, again)
    'AGC': 'S', 'AGT': 'S',
    # Valine (V)
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    # Alanine (A)
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    # Aspartic Acid (D)
    'GAC': 'D', 'GAT': 'D',
    # Glutamic Acid (E)
    'GAA': 'E', 'GAG': 'E',
    # Glycine (G)
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G'
    }  
    amino_acid_sequence = []
    for i in range(0, len(adn_sequence), 3):
        codon = adn_sequence[i:i+3]
        amino_acid = codon_table.get(codon, '')
        if amino_acid:
            amino_acid_sequence.append(amino_acid)
        else:
            print(f"No amino acid found for DNA sequence: {codon}")
    return ' '.join(amino_acid_sequence)


# Print a matrix in a readable format
def print_matrix(matrix):
    for row in matrix:
        print(' '.join(str(elem) for elem in row))
    print()
# Convert a matrix to a string representation
def matrix_to_string(matrix):
    """Convert a matrix to a string representation."""
    return '\n'.join([' '.join(map(str, row)) for row in matrix])

# Encrypt the amino acid sequence using the Hill cipher
def hill_encrypt(amino_acid_sequence, key_matrix):
    # Get matrix dimensions
    n, _ = key_matrix.shape
    encrypted_sequence = []
    for i in range(0, len(amino_acid_sequence), n):
        block = amino_acid_sequence[i:i + n]
        block_numbers = [ord(aa) - ord('A') for aa in block]
        block_matrix = np.array(block_numbers).reshape(n, -1)
        
        # Print the block matrix before encryption
        print("Block Matrix Before Encryption:")
        print(matrix_to_string(block_matrix))
        # Perform matrix multiplication and modulo operations
        encrypted_matrix = np.dot(key_matrix, block_matrix) % 26
        
        # Print the encrypted matrix
        print("Encrypted Matrix:")
        print(matrix_to_string(encrypted_matrix))
        # Convert the encrypted matrix back to a string of letters
        encrypted_block = ''.join([chr(ord('A') + num) for sublist in encrypted_matrix for num in sublist])
        encrypted_sequence.append(encrypted_block)
    return ' '.join(encrypted_sequence)


def matrix_to_string(matrix):
    """Convert a matrix to a string representation."""
    return '\n'.join([' '.join(map(str, row)) for row in matrix])

# Decrypt the encrypted sequence using the Hill cipher
def hill_decrypt(encrypted_sequence, key_matrix):
    determinant = int(round(np.linalg.det(key_matrix))) % 26
    
    # Check if determinant has a modular inverse
    if gcd(determinant, 26) != 1:
        print("Error: Determinant is not invertible modulo 26.")
        return ""
    # Compute modular inverse of the determinant
    inverse_determinant = pow(determinant, -1, 26)
    # Compute the adjugate matrix
    adjugate_matrix = (inverse_determinant * np.linalg.det(key_matrix) * np.linalg.inv(key_matrix)).astype(int) % 26
    
    n, _ = key_matrix.shape
    decrypted_sequence = []
    for block in encrypted_sequence:
        block_matrix = np.array([ord(c) - ord('A') for c in block]).reshape(n, -1)
        # Print the block matrix before decryption
        print("Block Matrix Before Decryption:")
        print(matrix_to_string(block_matrix))
        decrypted_matrix = np.dot(adjugate_matrix, block_matrix) % 26
        # Print the decrypted matrix
        print("Decrypted Matrix:")
        print(matrix_to_string(decrypted_matrix))
        decrypted_block = ''.join([chr(ord('A') + int(round(num))) for sublist in decrypted_matrix for num in sublist])
        decrypted_sequence.append(decrypted_block)
    return ' '.join(decrypted_sequence)


# Print a section header for better organization in the output
def print_section_header(title):
    print("\n" + "="*40)
    print(title)
    print("="*40)
# User input
plaintext = input("INSERT YOUR PLAINTEXT: ").lower()
matrix_size = int(input("INSERT MATRIX SIZE (2, 3, or 4): "))
key_matrix_str = input(f"INSERT YOUR {matrix_size}X{matrix_size} MATRIX (format: 'a b ...' without quotes): ")
key_matrix_list = list(map(int, key_matrix_str.split()))
key_matrix = np.array(key_matrix_list).reshape(matrix_size, matrix_size)
# Convert plaintext to different representations
binary_text = text_to_8bit_binary(plaintext)
dna_sequence = binary_to_adn(binary_text)
amino_acid_sequence = adn_to_amino_acid(dna_sequence)
# Display binary, DNA, and Amino Acid sequences
print_section_header("Input Sequences")
print(f"Binary Sequence: {binary_text}")
print(f"DNA Sequence: {dna_sequence}")
print(f"Amino Acid Sequence: {amino_acid_sequence}")
# Encrypt using Hill Cipher
encrypted_sequence = hill_encrypt(amino_acid_sequence.split(), key_matrix)
# Display encrypted sequence
print_section_header("Encryption Process")
print(f"Encrypted Sequence: {encrypted_sequence}")
# Decrypt using Hill Cipher
decrypted_sequence = hill_decrypt(encrypted_sequence.split(), key_matrix)
# Display decrypted sequence
print_section_header("Decryption Process")
print(f"Decrypted Sequence: {decrypted_sequence}")





