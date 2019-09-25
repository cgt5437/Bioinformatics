
'''
translate - a subroutine that translates an RNA sequence into its
corresponding amino acid sequence

Parameter: an RNA sequence that includes only one line of nucleotides
formatted through readfasta.py. Must include only the characters: 'A', 'C',
'G', and 'U'.
'''
def translate(sequence):
    # Create an empty amino acid sequence list to copy translations into
    aasequence = ''

    # Establish a code of trios of nucleotides and their corresponding single-letter amino acids
    code = {'UCA' : 'S',    # Serine
            'UCC' : 'S',    # Serine
            'UCG' : 'S',    # Serine
            'UCU' : 'S',    # Serine
            'UUC' : 'F',    # Phenylalanine
            'UUU' : 'F',    # Phenylalanine
            'UUA' : 'L',    # Leucine
            'UUG' : 'L',    # Leucine
            'UAC' : 'Y',    # Tyrosine
            'UAU' : 'Y',    # Tyrosine
            'UAA' : '_',    # Stop
            'UAG' : '_',    # Stop
            'UGC' : 'C',    # Cysteine
            'UGU' : 'C',    # Cysteine
            'UGA' : '_',    # Stop
            'UGG' : 'W',    # Tryptophan
            'CUA' : 'L',    # Leucine
            'CUC' : 'L',    # Leucine
            'CUG' : 'L',    # Leucine
            'CUU' : 'L',    # Leucine
            'CCA' : 'P',    # Proline
            'CCC' : 'P',    # Proline
            'CCG' : 'P',    # Proline
            'CCU' : 'P',    # Proline
            'CAC' : 'H',    # Histidine
            'CAU' : 'H',    # Histidine
            'CAA' : 'Q',    # Glutamine
            'CAG' : 'Q',    # Glutamine
            'CGA' : 'R',    # Arginine
            'CGC' : 'R',    # Arginine
            'CGG' : 'R',    # Arginine
            'CGU' : 'R',    # Arginine
            'AUA' : 'I',    # Isoleucine
            'AUC' : 'I',    # Isoleucine
            'AUU' : 'I',    # Isoleucine
            'AUG' : 'M',    # Methionine
            'ACA' : 'T',    # Threonine
            'ACC' : 'T',    # Threonine
            'ACG' : 'T',    # Threonine
            'ACU' : 'T',    # Threonine
            'AAC' : 'N',    # Asparagine
            'AAU' : 'N',    # Asparagine
            'AAA' : 'K',    # Lysine
            'AAG' : 'K',    # Lysine
            'AGC' : 'S',    # Serine
            'AGU' : 'S',    # Serine
            'AGA' : 'R',    # Arginine
            'AGG' : 'R',    # Arginine
            'GUA' : 'V',    # Valine
            'GUC' : 'V',    # Valine
            'GUG' : 'V',    # Valine
            'GUU' : 'V',    # Valine
            'GCA' : 'A',    # Alanine
            'GCC' : 'A',    # Alanine
            'GCG' : 'A',    # Alanine
            'GCU' : 'A',    # Alanine
            'GAC' : 'D',    # Aspartic Acid
            'GAU' : 'D',    # Aspartic Acid
            'GAA' : 'E',    # Glutamic Acid
            'GAG' : 'E',    # Glutamic Acid
            'GGA' : 'G',    # Glycine
            'GGC' : 'G',    # Glycine
            'GGG' : 'G',    # Glycine
            'GGU' : 'G'}    # Glycine

    # Go through each trio of nucleotides in the RNA sequence
    for i in range(0, len(sequence), 3):
        # Add the corresponding single-letter amino acid to the amino acid sequence
        aasequence += code[sequence[i:i+3]]

        # If the stop codon is reached, the translation is finished
        if code[sequence[i:i+3]] == '_':
            break

    # Return the resulting amino acid sequence
    return aasequence


    