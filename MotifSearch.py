from random import *
from profile import *


def greedy_motif_search(dna, k, t):
    BestMotifs = []
    for i in range(t):
        BestMotifs.append(dna[i][0:k])
    for i in range(len(dna[0]) - k + 1):
        motifs = [dna[0][i:i + k]]
        for j in range(1, t):
            p = profile(motifs[0:j])
            motifs.append(ProfileMostProbable(dna[j], k, p))
        if score(motifs) < score(BestMotifs):
            BestMotifs = motifs
    return BestMotifs


# get list of randomized motifs from each string in DNA
def random_motifs(dna, k, t):
    m = list()
    for i in range(t):
        r = randint(0, len(dna[i])-k)
        m.append(dna[i][r:r + k])
    return m


def random_motif_search(dna, k, t):
    m = random_motifs(dna, k, t)
    BestMotifs = m
    while True:
        p = profile(m)
        m = motifs(p, dna)
        if score(m) < score(BestMotifs):
            BestMotifs = m[:]
        else:
            return BestMotifs


# generate list of motifs where each motif corresponds to each line of dna
def motifs(p, dna):
    m = list()
    t = len(dna)
    k = len(p["A"])
    for i in range(t):
        m.append(ProfileMostProbable(dna[i], k, p))
    return m


def repeat(dna, k, t):
    b = float("inf")
    best = list()
    for i in range(1000):
        m = random_motif_search(dna, k, t)
        s = score(m)
        if s < b:
            b = s
            best = m
    return best


def gibbs_sampler(dna, k, t, N):
    m = random_motifs(dna, k, t)
    best = m
    for i in range(N):
        i = randint(0, t-1)
        m.pop(i)
        p = profile(m)
        m.insert(i, ProfileMostProbable(dna[i], k, p))
        if score(m) < score(best):
            best = m[:]
    return best


def repeat2(dna, k, t, N):
    b = float("inf")
    best = list()
    for i in range(100):
        m = gibbs_sampler(dna, k, t, N)
        s = score(m)
        if s < b:
            b = s
            best = m
    return best
