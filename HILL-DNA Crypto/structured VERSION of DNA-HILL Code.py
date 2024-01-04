import numpy as np

def text_to_8bit_binary(text):
    """Convert text to its 8-bit binary representation."""
    return ''.join(format(ord(c), '08b') for c in text)

def binary_to_adn(binary):
    """Convert binary string to DNA sequence."""
    mapping = {'00': 'A', '01': 'C', '10': 'G', '11': 'T'}
    while len(binary) % 2 != 0:
        binary += '0'
    return ''.join(mapping[binary[i:i+2]] for i in range(0, len(binary), 2))

def adn_to_amino_acid(adn_sequence):
    """Convert DNA sequence to amino acid sequence."""
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
    return ' '.join(codon_table.get(adn_sequence[i:i+3], '') for i in range(0, len(adn_sequence), 3))

def matrix_to_string(matrix):
    """Convert a matrix to a string representation."""
    return '\n'.join([' '.join(map(str, row)) for row in matrix])

def get_key_matrix(matrix_size):
    """Get the key matrix from user input."""
    while True:
        try:
            key_matrix_str = input(f"INSERT YOUR {matrix_size}X{matrix_size} MATRIX (format: 'a b ...' without quotes): ")
            key_matrix_list = list(map(int, key_matrix_str.split()))
            key_matrix = np.array(key_matrix_list).reshape(matrix_size, matrix_size)
            return key_matrix
        except ValueError:
            print("Invalid input. Please enter the matrix values again.")

def hill_encrypt(amino_acid_sequence, key_matrix):
    """Encrypt an amino acid sequence using the Hill Cipher."""
    n, _ = key_matrix.shape
    encrypted_sequence = []
    for i in range(0, len(amino_acid_sequence), n):
        block = amino_acid_sequence[i:i + n]
        block_numbers = [ord(aa) - ord('A') for aa in block]
        block_matrix = np.array(block_numbers).reshape(n, -1)
        encrypted_matrix = np.dot(key_matrix, block_matrix) % 26
        encrypted_block = ''.join([chr(ord('A') + num) for sublist in encrypted_matrix for num in sublist])
        encrypted_sequence.append(encrypted_block)
    return ' '.join(encrypted_sequence)

def hill_decrypt(encrypted_sequence, key_matrix):
    """Decrypt an encrypted sequence using the Hill Cipher."""
    determinant = int(round(np.linalg.det(key_matrix))) % 26
    if np.gcd(determinant, 26) != 1:
        raise ValueError("Error: Determinant is not invertible modulo 26.")
    inverse_determinant = pow(determinant, -1, 26)
    adjugate_matrix = (inverse_determinant * np.linalg.det(key_matrix) * np.linalg.inv(key_matrix)).astype(int) % 26
    
    n, _ = key_matrix.shape
    decrypted_sequence = []
    for block in encrypted_sequence:
        block_matrix = np.array([ord(c) - ord('A') for c in block]).reshape(n, -1)
        decrypted_matrix = np.dot(adjugate_matrix, block_matrix) % 26
        decrypted_block = ''.join([chr(ord('A') + int(round(num))) for sublist in decrypted_matrix for num in sublist])
        decrypted_sequence.append(decrypted_block)
    return ' '.join(decrypted_sequence)

def main():
    plaintext = input("INSERT YOUR PLAINTEXT: ").lower()
    matrix_size = int(input("INSERT MATRIX SIZE (2, 3, or 4): "))
    key_matrix = get_key_matrix(matrix_size)
    binary_text = text_to_8bit_binary(plaintext)
    dna_sequence = binary_to_adn(binary_text)
    amino_acid_sequence = adn_to_amino_acid(dna_sequence)
    print("Input Sequences")
    print(f"Binary Sequence: {binary_text}")
    print(f"DNA Sequence: {dna_sequence}")
    print(f"Amino Acid Sequence: {amino_acid_sequence}")
    encrypted_sequence = hill_encrypt(amino_acid_sequence.split(), key_matrix)
    print("\nEncryption Process")
    print(f"Encrypted Sequence: {encrypted_sequence}")
    decrypted_sequence = hill_decrypt(encrypted_sequence.split(), key_matrix)
    print("\nDecryption Process")
    print(f"Decrypted Sequence: {decrypted_sequence}")
if __name__ == "__main__":
    main()
