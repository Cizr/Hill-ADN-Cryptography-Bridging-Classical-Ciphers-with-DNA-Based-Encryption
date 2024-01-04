# **Hill-ADN Cryptography: A Fusion of Classical Ciphers and DNA Encryption**

## **Introduction**

Welcome to the Hill-ADN Cryptography repository! This innovative project amalgamates traditional Hill cipher methodologies with DNA encryption techniques, aiming to revolutionize data security. By harnessing the intricate properties of DNA sequences, this approach offers a robust and advanced encryption mechanism.

![Project Phases Flowchart](https://github.com/Cizr/Hill-ADN-Cryptography-Bridging-Classical-Ciphers-with-DNA-Based-Encryption/assets/100844208/b4d253cb-a53f-4550-8b9a-8c8e51e784e3)

## **Encryption Overview**

### **1. Binary and DNA Transformation**

- **Binary Conversion**: The plaintext is initially converted into an 8-bit binary representation.
- **DNA Encoding**: Subsequently, the binary data is mapped to DNA sequences using the rule: `A(00), C(01), G(10), T(11)`.
  
  **Example**: 'foo' in binary becomes `01100110 01101111 01101111`, and this is further encoded as the DNA sequence `CGCGCGTTCGTT`.

### **2. Amino Acid Translation**

![Amino Acid Mapping](https://github.com/Cizr/Hill-ADN-Cryptography-Bridging-Classical-Ciphers-with-DNA-Based-Encryption/assets/100844208/99812733-cc98-4b4b-86ef-f4f0fe5119e3)

- **Mapping to Amino Acids**: The DNA sequences, composed of triplets (e.g., CGT), are translated into amino acids based on the genetic codon table. For instance, CGT corresponds to the amino acid Arginine (R).

## **Hill Cipher Encryption with Amino Acid Sequence**

- **Input**: The amino acid sequence, derived from the previous step, serves as the input plaintext for the Hill cipher encryption.
- **Key Setup**: A key matrix is established based on the provided encryption key, tailored to the encryption algorithm's specific parameters.
- **Encryption Process**: 
  - Divide the amino acid sequence into blocks, each aligned with the key matrix's dimensions.
  - Convert each block into a numerical matrix using the amino acid-to-numeric mapping.
  - Perform matrix multiplication between the key matrix and the blocks, followed by modulo operations to fit within the desired alphabet size (e.g., 26 for English alphabets).
  - The resulting matrices form the ciphertext, which is then translated back into an amino acid sequence.

## **Decryption Overview**

### **Hill-ADN Decryption**

- **Input**: The encrypted amino acid sequence undergoes decryption using the Hill cipher.
- **Decryption Steps**: 
  - Utilize the Hill cipher decryption technique to decrypt the encrypted amino acid sequence.
  - The decrypted result yields the original amino acid sequence as the output.

## **Conclusion**

The Hill-ADN Cryptography repository presents a novel encryption paradigm by integrating classical Hill cipher techniques with DNA-based encoding. By leveraging the intricacies of biological structures, this innovative approach aspires to advance data security and encryption methodologies.

