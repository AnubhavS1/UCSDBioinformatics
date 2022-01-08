from collections import defaultdict, Counter
from itertools import product


# find all kmers in text
def composition(text, k):
    kmers = list()
    for i in range(len(text) - k + 1):
        kmers.append(text[i:i+k])
    #kmers.sort()
    return kmers


# form a DNA string based on pieces of the string given in correct order
def overlap(pieces):
    final = ''
    final += pieces[0]
    for i in range(1, len(pieces)):
        final += pieces[i][-1]
    return final


# form an adjacency graph for a list of kmers
def overlap_graph(kmers):
    adjacency_list = defaultdict(set)
    for kmer in kmers:
        adjacency_list[kmer[:-1]].add(kmer)
    for kmer in kmers:
        suffixes = adjacency_list[kmer[1:]]
        if suffixes:
            print(kmer + " -> " + ",".join(suffixes))


# form deBruijin graph
def deBruijin_graph(sequences):
    d = dict()
    for seq in sorted(sequences):
        if seq[:-1] in d:
            d[seq[:-1]].append(seq[1:]) 
        else:
            d[seq[:-1]] = [seq[1:]]
    graph = [' -> '.join([item[0], ','.join(item[1])]) for item in d.items()]
    return d


# assemble a Eulerian cycle based on a Eulerian graph
def eulerian_cycle(graph):
    current = list(graph)[0]
    cycle = [current]
    while True:
        # add the node current was pointing to into cycle
        cycle.append(graph[current][0])
        # remove the node we just add to cycle from graph[current]
        if len(graph[current]) != 1:
            graph[current] = graph[current][1:]
        else:
            del graph[current]
        # make the node we just added to cycle current
        if cycle[-1] in graph:
            current = cycle[-1]
        else:
            break
    # iterate over cycles from unexplored nodes
    while len(graph) > 0:
        for i in range(len(cycle)):
            if cycle[i] in graph:
                current = cycle[i]
                temp = [current]
                # repeat process from previous loop on the unexplored regions
                while True:
                    temp.append(graph[current][0])
                    if len(graph[current]) != 1:
                        graph[current] = graph[current][1:]
                    else:
                        del graph[current]
                    x = temp[-1]
                    if x in graph.keys():
                        current = x
                    else:
                        break
                cycle = cycle[:i] + temp + cycle[i+1:]
                break
    return cycle


# make an Eulerian path based on a given graph
def eulerian_path(graph):
    # determine for which node do the input_edges != outputs_edges
    global u_to, u_from
    nodes = []   # set for all the nodes, in case one of the dicts doesn't contain it
    in_edges = {}
    for value in graph.values():
        for node in value:
            if node in in_edges:
                in_edges[node].append(1)
            else:
                in_edges[node] = [1]
                nodes.append(node)

    out_edges = {}
    for key in graph.keys():
        out_edges[key] = [len(graph[key])]
        if key not in nodes:
            nodes.append(key)

    for node in nodes:
        if node not in in_edges:
            in_edges[node] = [0]
        if node not in out_edges:
            out_edges[node] = [0]
        if in_edges[node] > out_edges[node]:
            u_from = node
        elif in_edges[node] < out_edges[node]:
            u_to = node

    # add an edge between the unbalanced to and from
    if u_from in graph:
        graph[u_from].append(u_to)
    else:
        graph[u_from] = [u_to]
    # get Eulerian cycle for the graph
    cycle = eulerian_cycle(graph)
    # check if an edge links both the unbalance from/to nodes
    break_point = list(filter(lambda i: cycle[i:i+2] == [u_from, u_to], range(len(cycle) - 1)))[0]

    return cycle[break_point+1:]+cycle[1:break_point+1]


def universal_string(k):
    d = {}
    t = [''.join(elem) for elem in product('01', repeat=k)]
    for kmer in t:
        if kmer[:-1] in d and kmer[:-1] != kmer[1:]:
            d[kmer[:-1]].append(kmer[1:])
        else:
            if kmer[:-1] != kmer[1:]:
                d[kmer[:-1]] = [kmer[1:]]
    return ''.join(eulerian_cycle(d))


# build debruijn graph based on kmer pairs seperated by d bp (XXXXXX|XXXXXX)
def paired_deBruijn(pair_dict):
    graph = {}
    for key, value in pair_dict.items():
        for k, v in pair_dict.items():
            if key == v:
                if v in graph and value not in graph:
                    graph[v].append(value)
                else:
                    graph[k] = [v]
                    graph[v] = [value]
    return graph


# construct a string from a list of gapped patterns
def string_from_gapped_patterns(pairs, k, d):
    first = pairs[0][0] + ''.join([pair[0][-1] for pair in pairs[1:]])
    second = pairs[0][1] + ''.join([pair[1][-1] for pair in pairs[1:]])
    if first[k+d:] == second[:-k-d]:
        return first + second[-(k+d):]
    return 'None'


# put the data from the codon table into a dictionary for use
t = open('RNA_codon_table_1.txt', 'r')
codons = {}
tmp = [line.strip().split() for line in t]
for translation in tmp:
    try:
        codons[translation[0]] = translation[1]
    except IndexError:
        pass
t.close()


# translate RNA into amino acid string
def translate_RNA(rna):
    amino_acids = ''
    for i in range(0, len(rna), 3):
        for c in codons.keys():
            codon = rna[i:i+3]
            if codon == c:
                amino_acids += codons[c]
    return amino_acids


# get reverse complement
def rev_comp(dna):
    d = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G'}
    return ''.join(d[base] for base in dna[::-1])


# find all substrings in a DNA string that encode for the desired peptide
def peptide_encoding(dna, peptide):
    substrings = list()
    for i in range(len(dna) - (3*len(peptide)) + 1):
        sub = dna[i:3*(len(peptide))+i].replace('T', 'U')
        reverse = rev_comp(dna[i:3*(len(peptide))+i]).replace('T', 'U')
        if translate_RNA(sub) == peptide or translate_RNA(reverse) == peptide:
            substrings.append(dna[i:3*len(peptide)+i])
    return substrings


# construct a cyclic spectrum for a peptide chain
def cyclospectrum(peptide, mass_dict):
    prefix_mass = [0 for i in range(len(peptide) + 1)]
    for i in range(1, len(peptide) + 1):
        prefix_mass[i] = prefix_mass[i-1] + int(mass_dict[peptide[i-1]])
    spectrum = [0]
    for i in range(len(peptide)):
        for j in range(i+1, len(peptide)+1):
            spectrum.append(prefix_mass[j] - prefix_mass[i])
            if i > 0 and j < len(peptide):
                spectrum.append(prefix_mass[-1] - (prefix_mass[j] - prefix_mass[i]))
    return sorted(spectrum)


# get the linear spectrum of a peptide
def linear_spectrum(peptide):
    mass_dict = get_aa_mass()
    prefix_mass = [0 for i in range(len(peptide) + 1)]
    for i in range(1, len(peptide) + 1):
        prefix_mass[i] = prefix_mass[i-1] + int(mass_dict[peptide[i-1]])
    spectrum = [0]
    for i in range(len(peptide)):
        for j in range(i+1, len(peptide)+1):
            spectrum.append(prefix_mass[j] - prefix_mass[i])
    return sorted(spectrum)


# count the number of peptides that have a given mass
def count_mass(mass, mass_dict, alpha):
    if mass == 0:
        return 1, mass_dict
    if mass < 57:
        return 0, mass_dict
    if mass in mass_dict:
        return mass_dict[mass], mass_dict
    n = 0
    for a in alpha:
        k, mass_dict = count_mass(mass - int(alpha[a]), mass_dict, alpha)
        n += k
    mass_dict[mass] = n
    return n, mass_dict


# sequence a cyclopeptide from an experimental spectrum ################################################################
def get_aa_mass():
    aa = {}
    with open('integer_mass_table.txt') as f:
        for line in f:
            aa[line.split()[0]] = line.strip().split()[1]
    return aa


# check if amino acids are inconsistent with the spectrum and remove if so
def inconsistent(peptide, spectrum):
    spec = spectrum[:]
    linear = linear_spectrum(peptide)
    for p in linear:
        if p not in spec:
            return True
        spec.remove(p)
    return False


# find the mass of a peptide
def mass(peptide):
    m = 0
    for aa in peptide[:]:
        m += int(get_aa_mass()[aa])
    return m


# return array of masses of amino acids from a peptide
def convert_to_masses(peptide):
    c = []
    for aa in peptide:
        c.append(get_aa_mass()[aa])
    return c


# expand an array of peptides to include all peptides that are consistent with a mass spectrum
def expand(peptides, spectrum):
    expanded_peptides = []
    if not peptides:
        return list(get_aa_mass().keys())
    for peptide in peptides:
        for key in get_aa_mass().keys():
            if not inconsistent(peptide + key, spectrum):
                expanded_peptides.append(peptide + key)
    return expanded_peptides


# Compute the number of masses shared with the experimental and theoretical spectrums of a peptide
def score(peptide, spectrum):
    s = 0
    spec = cyclospectrum(peptide, get_aa_mass())
    tcount = Counter(spec)
    ecount = Counter(spectrum)
    for mass in ecount.keys():
        if mass in tcount.keys():
            s += min(ecount[mass], tcount[mass])
    return s


def linear_score(peptide, spectrum):
    s = 0
    spec = linear_spectrum(peptide)
    tcount = Counter(spec)
    ecount = Counter(spectrum)
    for mass in ecount.keys():
        if mass in tcount.keys():
            s += min(ecount[mass], tcount[mass])
    return s


# trim a leader board to include the peptides with the highest score
def trim(board, spectrum, n):
    linear_scores = []
    for i in range(len(board)):
        linear_scores.append(linear_score(board[i], spectrum))
    board = [x for _,x in sorted(zip(linear_scores, board), reverse=True)]
    linear_scores = sorted(linear_scores, reverse=True)
    for i in range(n+1, len(board)):
        if linear_scores[i] < linear_scores[n]:
            if i % 2 == 1:
                del board[i-1:]
            else:
                del board[i:]
            return board
    return board


# output the convolution list of a spectrum
def convolution(spectrum):
    return [i-j for i in spectrum for j in spectrum if i - j > 0]
