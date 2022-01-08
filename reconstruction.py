from scripts import *


# reconstruct a string based on a list of pieces
def string_reconstruction(patterns, k):
    dB = deBruijin_graph(patterns)
    path = eulerian_path(dB)
    string = ''
    for node in path:
        if node == path[0]:
            string += node
        else:
            string += node[-1]
    return string


def pair_comp_graph(d, pairs):
    pair_dict = {}
    for elem in pairs:
        pre = (elem[0][:-1], elem[1][:-1])
        suff = (elem[0][1:], elem[1][1:])
        pair_dict[pre] = suff
    return pair_dict


# reconstruct a string from pairs of kmers
def pair_reconstruction(k, d, pairs):
    pair_graph = pair_comp_graph(d, pairs)
    dB = paired_deBruijn(pair_graph)
    path = eulerian_path(dB)
    return string_from_gapped_patterns(path, k, d)