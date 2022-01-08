import HammingDistance

store = list()
nucleotides = {'A', 'C', 'G', 'T'}


def neighbors(pattern, d):
    # pattern = pattern[0:len(pattern) - 1]
    if d == 0:
        return pattern
    # base case
    if len(pattern) == 1:
        return {'A', 'C', 'G', 'T'}
    neighborhood = list()
    suffixNeighbors = neighbors(suffix(pattern), d)

    for text in suffixNeighbors:
        suff = suffix(pattern)
        dist = HammingDistance.HammingDistance(suff, text)
        if dist < int(d):
            for nucleotide in {'A', 'C', 'G', 'T'}:
                neighborhood.append(nucleotide + text)
        else:
            neighborhood.append(first(pattern) + text)
    return neighborhood


def suffix(pattern):
    return pattern[1:len(pattern)]


def first(pattern):
    return pattern[0]