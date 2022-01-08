import skew
from HammingDistance import HammingDistance
import ApproxPattern
import FreqWords
import Neighborhood
import MotifEnumeration
from MedianString import *
from profile import *
from MotifSearch import *


file = open("test.txt", "r")

def testSkew():
    gen = file.readline()
    for elem in skew.skew(gen):
        print(elem, end=" ")


def testHamm():
    p1 = file.readline()
    q1 = file.readline()
    print(HammingDistance(p1, q1))


def testApprox():
    pat = file.readline()
    genome = file.readline()
    D = file.readline()
    for elem in ApproxPattern.approx(pat, genome, D):
        print(elem, end=" ")
    print()
    print(ApproxPattern.approxCount(pat, genome, D))


def testMismatch():
    genome = file.readline()
    temp = file.readline()
    K = temp.split()[0]
    D = temp.split()[1]
    for key in FreqWords.mismatches(genome, K, D):
        print(key, end=" ")


def testNeighbor():
    pattern = file.readline()
    D = file.readline()
    for elem in Neighborhood.neighbors(pattern, D):
        print(elem, end="\n")


def testMotif():
    temp = file.readline()
    k = temp.split()[0]
    d = temp.split()[1]
    DNA = list()
    for line in file:
        if line[len(line)-1] == '\n':
            DNA.append(line[0:len(line) - 1])
        else:
            DNA.append(line[0:len(line)])
    for elem in MotifEnumeration.MotifEnumeration(DNA, k, d):
        print(elem, end=" ")


def testDistance():
    patt = file.readline()
    DNA = file.readline().split()
    print(d(patt, DNA))


def testMedian():
    K = file.readline()
    DNA = list()
    for line in file:
        DNA.append(line.strip("\n"))
    for elem in MedianString(DNA, K):
        print(elem, end=" ")
    #print(MedianString(DNA, K))


def testProfile():
    DNA = file.readline().strip("\n")
    K = file.readline().strip("\n")
    matrix = {'A': None, 'C' : None, 'G' : None, 'T' : None}
    for nuc in matrix:
        elems = file.readline().split()
        for i in range(len(elems)):
            elems[i] = float(elems[i])
        matrix[nuc] = elems
    #print(ProfileMostProbable(DNA, K, matrix))
    print(prob("GAGCTA", matrix))


def testGreedy():
    temp = file.readline()
    k = temp.split()[0]
    t = temp.split()[1]
    DNA = {}
    for i in range(int(t)):
        DNA[i] = file.readline().strip("\n")
    for elem in greedy_motif_search(DNA, int(k), int(t)):
        print(elem, end="\n")


def testRandom():
    temp = file.readline().split()
    k = int(temp[0])
    t = int(temp[1])
    DNA = {}
    for i in range(t):
        DNA[i] = file.readline().strip("\n")
    for elem in repeat(DNA, k, t):
        print(elem, end="\n")


def testGibbs():
    temp = file.readline().split()
    k = int(temp[0])
    t = int(temp[1])
    N = int(temp[2])
    DNA = {}
    for i in range(t):
        DNA[i] = file.readline().strip("\n")
    for elem in repeat2(DNA, k, t, N):
        print(elem, end="\n")


# testSkew()
testHamm()
# testApprox()
# testMismatch()
# testNeighbor()
# testMotif()
# testDistance()
# testMedian()
# testProfile()
# testGreedy()
# testRandom()
# testGibbs()
