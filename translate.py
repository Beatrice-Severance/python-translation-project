#! /usr/bin/env python3

import sys

def translate_sequence(rna_sequence, genetic_code):
    """Translates a sequence of RNA into a sequence of amino acids.

    Translates `rna_sequence` into string of amino acids, according to the
    `genetic_code` given as a dict. Translation begins at the first position of
    the `rna_sequence` and continues until the first stop codon is encountered
    or the end of `rna_sequence` is reached.

    If `rna_sequence` is less than 3 bases long, or starts with a stop codon,
    an empty string is returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    str
        A string of the translated amino acids.
    """
    rna = rna_sequence.upper()
    translation = ""
    if len(rna) >= 3:
        if len(rna)%3 == 1:
            rna = rna[:-1]
        for input in range (0, len(rna), 3):
            codon = rna[input:input + 3]
            translation+= genetic_code[codon]
        if "*" in str(translation):
            separate = "*"
            translation = translation.split(separate, 1)[0]
            translation = translation.replace("*", "")
    return translation

def get_all_translations(rna_sequence, genetic_code):
    """Get a list of all amino acid sequences encoded by an RNA sequence.

    All three reading frames of `rna_sequence` are scanned from 'left' to
    'right', and the generation of a sequence of amino acids is started
    whenever the start codon 'AUG' is found. The `rna_sequence` is assumed to
    be in the correct orientation (i.e., no reverse and/or complement of the
    sequence is explored).

    The function returns a list of all possible amino acid sequences that
    are encoded by `rna_sequence`.

    If no amino acids can be translated from `rna_sequence`, an empty list is
    returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    list
        A list of strings; each string is an sequence of amino acids encoded by
        `rna_sequence`.
    """
    rna_seqs = rna_sequence.upper()
    aa_list = []
    start_pos = 0
    def translate(start_pos, rna_seqs, genetic_code):
        proteins = ''
        for input in range(start_pos, len(rna_seqs), 3):
            codon = rna_seqs[input:input + 3]
            if codon in ['UAG', 'UAA', 'UGA'] or len(codon) != 3:
                break
            else: proteins += genetic_code[codon]
        return proteins
    while start_pos < len(rna_seqs):
        start_codon = rna_seqs[start_pos:start_pos + 3]
        if start_codon == 'AUG':
            translation = translate(start_pos, rna_seqs, genetic_code)
            aa_list.append(translation)
        start_pos += 1
    return aa_list
    """rna = rna_sequence.upper()
    translation = []
    start = rna.find("AUG")
    starttrans = rna[int(start):]
    for input in range (0, len(starttrans), 3):
        if starttrans[input:input + 3] in genetic_code:
            translation+= genetic_code[starttrans[input:input + 3]]
        if start == 0:
            translation = []
    return translation"""

def get_reverse(sequence):
    """Reverse orientation of `sequence`.

    Returns a string with `sequence` in the reverse order.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> get_reverse('AUGC')
    'CGUA'
    >>> get_reverse('ATGC')
    'CGTA'
    """
    sequence = list(sequence.upper())
    print("The sequence is: ", sequence)
    reverseseq = sequence [::-1]
    reverse = ("".join(reverseseq))
    return reverse

def get_complement(sequence):
    """Get the complement of a `sequence` of nucleotides.

    Returns a string with the complementary sequence of `sequence`.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> get_reverse('AUGC')
    'UACG'
    >>> get_reverse('ATGC')
    'TACG'
    """
    sequence = list(sequence.upper())
    print("The sequence is: ", sequence)
    complement = ""
    for input in sequence:
        if input == "A":
            complement += "U"
        if input == "U":
            complement += "A"
        if input == "C":
            complement += "G"
        if input == "G":
            complement += "C"
    return complement

def reverse_and_complement(sequence):
    """Get the reversed and complemented form of a `sequence` of nucleotides.

    Returns a string that is the reversed and complemented sequence
    of `sequence`.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> reverse_and_complement('AUGC')
    'GCAU'
    >>> reverse_and_complement('ATGC')
    'GCAT'
    """
    sequence = list(sequence.upper())
    print("The sequence is: ", sequence)
    reverseseq = sequence [::-1]
    reverse = ("".join(reverseseq))
    print(reverse)
    complement = ""
    for input in reverse:
        if input == "A":
            complement += "U"
        if input == "U":
            complement += "A"
        if input == "C":
            complement += "G"
        if input == "G":
            complement += "C"
    return complement

def get_longest_peptide(rna_sequence, genetic_code):
    """Get the longest peptide encoded by an RNA sequence.

    Explore six reading frames of `rna_sequence` (the three reading frames of
    `rna_sequence`, and the three reading frames of the reverse and complement
    of `rna_sequence`) and return (as a string) the longest sequence of amino
    acids that it encodes, according to the `genetic_code`.

    If no amino acids can be translated from `rna_sequence` nor its reverse and
    complement, an empty string is returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    str
        A string of the longest sequence of amino acids encoded by
        `rna_sequence`.
    """
    peps = get_all_translations(rna_sequence = rna_sequence, genetic_code = genetic_code)
    rev_and_comp = reverse_and_complement(rna_sequence)
    rev_tran = get_all_translations(rna_sequence = rev_and_comp, genetic_code = genetic_code)
    peps += rev_tran
    if not peps:
        return ""
    if len(peps) < 2:
        return peps[0]
    most_bases = -1
    longest_pep = -1
    for pep, aa_list in enumerate(peps):
        if len(aa_list) > most_bases:
            longest_pep = pep
            most_bases = len(aa_list)
    return peps[longest_pep]
    """rna_seqs = rna_sequence.upper()
    aa_list = []
    start_pos = 0
    def translate(start, rna, genetic_code):
        proteins = ""
        for input in range(start, len(rna), 3):
            codon = rna[input:input + 3]
            if codon in ["UAG", "UAA", "UGA"] or len(codon) != 3:
                break
            else: proteins += genetic_code[codon]
        return proteins
    def reverse_and_complement(rna_seqs):
        rna_seqs = list(rna_seqs.upper())
        rev_seqs = rna_seqs[::-1]
        reversed = ("".join(rev_seqs))
        complement = ""
        for i in reversed:
            if i == "A":
                complement += "U"
            if i == "U":
                complement += "A"
            if i == "C":
                complement += "G"
            if i == "G":
                complement += "C"
        return complement
    while start_pos < len(rna_seqs):
        start_codon = rna_seqs[start_pos:start_pos + 3]
        if start_codon == 'AUG':
            translation = translate(start_pos, rna_seqs, genetic_code)
            aa_list.append(translation)
            rev_and_comp = reverse_and_complement(rna_seqs)
            translation_back = translate(start_pos, rev_and_comp, genetic_code)
            aa_list.append(translation_back)
        start_pos += 1
    check = 'M'
    aa_list_start = [idx for idx in aa_list if idx[0].lower() == check.lower()]
    longest_aa = ''
    if aa_list_start == []:
        longest_aa == ''
    if aa_list_start != []:
        longest_aa = max((aa_list_start), key=len)
    return longest_aa"""

if __name__ == '__main__':
    genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}
    rna_seq = ("AUG"
            "UAC"
            "UGG"
            "CAC"
            "GCU"
            "ACU"
            "GCU"
            "CCA"
            "UAU"
            "ACU"
            "CAC"
            "CAG"
            "AAU"
            "AUC"
            "AGU"
            "ACA"
            "GCG")
    longest_peptide = get_longest_peptide(rna_sequence = rna_seq,
            genetic_code = genetic_code)
    assert isinstance(longest_peptide, str), "Oops: the longest peptide is {0}, not a string".format(longest_peptide)
    message = "The longest peptide encoded by\n\t'{0}'\nis\n\t'{1}'\n".format(
            rna_seq,
            longest_peptide)
    sys.stdout.write(message)
    if longest_peptide == "MYWHATAPYTHQNISTA":
        sys.stdout.write("Indeed.\n")
