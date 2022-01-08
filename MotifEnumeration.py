import Neighborhood
import HammingDistance


def MotifEnumeration(dna, k, d):
    patterns = set()
    test = list()
    k = int(k)
    d = int(d)
    for sequence in dna:
        for i in range(len(sequence) - k + 1):
            pattern = sequence[i:i+k]
            #for j in range(len(sequence) - k + 1):
            for neighbor in Neighborhood.neighbors(pattern, d):
                test = sequence[i:i+k]
                if HammingDistance.HammingDistance(neighbor, test) <= d:
                    if shared(dna, test, d):
                        patterns.add(test)
                    if shared(dna, neighbor, d):
                        patterns.add(neighbor)
    return patterns


def shared(dna, test, d):
    d = int(d)
    count = 0
    for line in dna:
        for i in range(len(line) - len(test) + 1):
            t = line[i:i+len(test)]
            ham = HammingDistance.HammingDistance(t, test)
            if ham <= d:
                count += 1
                break
    if count == len(dna):
        return True
    return False

