def dna_complement(dna_seq):
    """
    Computes the complement of a DNA sequence.

    :param dna_seq: str, the DNA sequence
    :return: str, the complementary DNA sequence
    """
    complement_table = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'
    }
    complement_dna = ''.join([complement_table[nuc] for nuc in dna_seq.upper()])
    return complement_dna


def dna_to_mrna(dna_seq):
    """
    Translates a DNA sequence into an mRNA sequence.
    Replaces 'T' with 'U' in the DNA sequence to form mRNA.

    :param dna_seq: str, the DNA sequence (must be a multiple of 3)
    :return: str, the corresponding mRNA sequence
    """
    if len(dna_seq) % 3 != 0:
        raise ValueError("DNA sequence length must be a multiple of 3")
    return dna_seq.upper().replace('T', 'U')


def translate_mrna_to_protein(mrna_seq):
    """
    Translates an mRNA sequence into an amino acid sequence (protein).

    :param mrna_seq: str, the mRNA sequence (must be a multiple of 3)
    :return: str, the corresponding amino acid sequence
    """
    codon_table = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGU': 'S', 'AGC': 'S',
        'AGA': 'R', 'AGG': 'R', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    if len(mrna_seq) % 3 != 0:
        raise ValueError("mRNA sequence length must be a multiple of 3")

    protein = []
    for i in range(0, len(mrna_seq), 3):
        codon = mrna_seq[i:i + 3]
        amino_acid = codon_table.get(codon, '?')  # '?' for unknown codons
        protein.append(amino_acid)

    return '-'.join(protein)


# Task 2 Start
def codon_frequency_for_aminoacid(amino_acids):
    """
    Calculates all possible mRNA sequences from given amino acids
    and the frequency of each codon encoding those amino acids.

    :param amino_acids: str, a string of amino acids (max length 3)
    :return: list of tuples, each containing an mRNA sequence and its codon frequencies
    """
    codon_table = {
        'F': ['UUU', 'UUC'], 'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
        'I': ['AUU', 'AUC', 'AUA'], 'M': ['AUG'], 'V': ['GUU', 'GUC', 'GUA', 'GUG'],
        'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'], 'P': ['CCU', 'CCC', 'CCA', 'CCG'],
        'T': ['ACU', 'ACC', 'ACA', 'ACG'], 'A': ['GCU', 'GCC', 'GCA', 'GCG'],
        'Y': ['UAU', 'UAC'], '*': ['UAA', 'UAG', 'UGA'], 'H': ['CAU', 'CAC'],
        'Q': ['CAA', 'CAG'], 'N': ['AAU', 'AAC'], 'K': ['AAA', 'AAG'],
        'D': ['GAU', 'GAC'], 'E': ['GAA', 'GAG'], 'C': ['UGU', 'UGC'],
        'W': ['UGG'], 'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'G': ['GGU', 'GGC', 'GGA', 'GGG']
    }

    if len(amino_acids) > 3:
        raise ValueError("Input amino acid sequence should be at most 3 characters long")

    # Generate all combinations of codons for the given amino acids
    from itertools import product
    codon_combinations = [codon_table[aa] for aa in amino_acids]
    all_mrna_sequences = list(product(*codon_combinations))

    results = []
    combined_codon_freq = {}

    for mrna_tuple in all_mrna_sequences:
        # Join codons to form an mRNA sequence
        mrna_seq = ''.join(mrna_tuple)
        results.append(mrna_seq)

        # Count the frequency of each codon
        for codon in mrna_tuple:
            if codon in combined_codon_freq:
                combined_codon_freq[codon] += 1
            else:
                combined_codon_freq[codon] = 1

    return results, combined_codon_freq


# Task 1: DNA -> mRNA -> Amino Acid Sequence
def task_1_process_dna_sequence(dna_sequence):
    print("\nThis is Task 1: Translating DNA to mRNA and Protein Sequence\n")

    # Input DNA sequence
    print(f"Input DNA: {dna_sequence}")

    # Complement DNA sequence
    complement_dna = dna_complement(dna_sequence)
    print(f"Complement DNA: {complement_dna}")

    # Convert DNA to mRNA
    mrna_sequence = dna_to_mrna(dna_sequence)
    print(f"mRNA: {mrna_sequence}")

    # Convert mRNA to Amino Acid Sequence
    protein_sequence = translate_mrna_to_protein(mrna_sequence)
    print(f"Amino Acid Sequence: {protein_sequence}")


# Task 2: Codon Frequency Calculation for Amino Acid
def task_2_codon_frequency(amino_acids):
    print("\nThis is Task 2: Calculating Codon Frequencies for Amino Acids\n")

    # Calculate all possible mRNA sequences and their corresponding codon frequencies
    mrna_sequences, codon_freq = codon_frequency_for_aminoacid(amino_acids)

    print(f"Input Amino Acid: {amino_acids}")
    for mrna in mrna_sequences:
        print(f"mRNA = {mrna}")

    print(f"Codon Frequencies: {codon_freq}")


# Example Execution for Task 1 and Task 2
if __name__ == "__main__":
    # Example for Task 1
    # dna_sequence = "TTACGA"
    dna_sequence = input("Please enter DNA sequence: ")
    task_1_process_dna_sequence(dna_sequence)

    # Example for Task 2
    # amino_acids = "WYW"
    amino_acids = input("Please enter amino acids: ")
    task_2_codon_frequency(amino_acids)
