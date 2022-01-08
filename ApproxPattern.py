# match pattern with approximate patterns in the genome/text
import HammingDistance


def approx(pattern, text, d):
    indices = list()
    pattern = pattern[0:len(pattern)-1]
    text = text[0:len(text)-1]
    for i in range(len(text) - len(pattern) + 1):
        if HammingDistance.HammingDistance(pattern, text[i:i+len(pattern)]) <= int(d):
            indices.append(i)
    return indices


def approxCount(pattern, text, d):
    count = 0
    for i in range(len(text) - len(pattern) + 1):
        test = HammingDistance.HammingDistance(pattern, text[i:i+len(pattern)-1])
        if test <= int(d):
            count += 1
    return count
    #return len(approx(pattern, text, d))
