import PatternCount
import Neighborhood
import ComputingFrequencies


def freqWords(text, k):
    patterns = list()
    counts = list()
    for i in range(len(text) - k):
        pattern = text[i : i+k]
        counts.append(PatternCount.count(text, pattern))
    maxCount = 0
    for val in counts:
        if val > maxCount:
            maxCount = val
    for i in range(len(text) - k):
        if counts[i] == maxCount:
            patterns.append(text[i : i+k])
    #remove dups
    final = list()
    for elem in patterns:
        if patterns.count(elem) == 1:
            final.append(elem)
        else:
            if final.count(elem) == 0:
                final.append(elem)
    return final


def mismatches(text, k, d):
    patterns = set()
    neighborhoods = [Neighborhood.neighbors(text[i:i+int(k)], d) for i in range(len(text) - int(k) + 1)]
    hoodsList = [elems for list in neighborhoods for elems in list]
    index = [ComputingFrequencies.PatternToNumber(hoodsList[i]) for i in range(len(hoodsList))]
    counts = [1] * len(index)
    sortedIndex = sorted(index)

    for i in range(len(neighborhoods) - 1):
        if sortedIndex[i] == sortedIndex[i+1]:
            counts[i+1] = counts[i] + 1
    maxCount = max(counts)
    for i in range(len(hoodsList)):
        if counts[i] == maxCount:
            patterns.add(ComputingFrequencies.NumberToPattern(sortedIndex[i], k))
    return patterns

