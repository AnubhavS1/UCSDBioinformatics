from HammingDistance import HammingDistance
from Neighborhood import neighbors


def d(pattern, dna):
    pattern = pattern[0:len(pattern)]
    k = len(pattern)
    distance = 0
    for text in dna:
        hamming = float("inf")
        for i in range(len(text)-k+1):
            Pattern = text[i:i+k]
            if hamming > HammingDistance(pattern, Pattern):
                hamming = HammingDistance(pattern, Pattern)
        distance += hamming
    return distance


def MedianString(dna, k):
    kmers = list()
    k = int(k)
    # add all kmers to list kmers
    for text in dna:
        line = text.strip("\n")
        for i in range(len(line) - k + 1):
            pattern = line[i:i + k]
            for neighbor in neighbors(pattern, 1):
                if neighbor not in kmers:
                    kmers.append(neighbor)

    distance = float("inf")
    medians = list()
    # find pattern that minimizes distance
    for kmer in kmers:
        test = d(kmer, dna)
        if distance > test:
            distance = d(kmer, dna)
            medians.append(kmer)
    return medians

