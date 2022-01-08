from scripts import *


# sequence a cyclopeptide based on a given experimental mass spectrum
def cyclopeptide_sequencing(spectrum):
    final = []
    f2 = []
    candidates = expand([], spectrum)
    while candidates:
        for peptide in candidates:
            if mass(peptide) == spectrum[-1]:
                if cyclospectrum(peptide, get_aa_mass()) == spectrum and convert_to_masses(peptide) not in final \
                   and peptide not in f2:
                    final.append(convert_to_masses(peptide))
                    f2.append(peptide)
        candidates = expand(candidates, spectrum)
    return final, f2


# output leader peptide based on a constructed leaderboard from mass spectrum
def leaderboard_sequencing(spectrum, n):
    board = ['']
    leader = ''
    while board:
        board = expand(board, spectrum)
        for peptide in board:
            if mass(peptide) == spectrum[-1]:
                if score(peptide, spectrum) > score(leader, spectrum):
                    leader = peptide
            if mass(peptide) > spectrum[-1]:
                board.pop(board.index(peptide))
        board = trim(board, spectrum, n)
    return leader, convert_to_masses(leader)


# sequence a cyclopeptide using convolution
def convolution_sequencing(m, n, spectrum):
    con = convolution(spectrum)
    con_dict = {}
    for c in set(filter(lambda c: 57 <= c <= 200, con)):
        num_c = con.count(c)
        if num_c in con_dict:
            con_dict[num_c].append(c)
        else:
            con_dict[num_c] = [c]
