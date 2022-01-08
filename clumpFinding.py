import PatternCount
import FreqWords


def clump(text, k, L, t):
    patterns = list()
    for i in range(len(text) - k):
        p = FreqWords.freqWords(text[i:i+L], k)
        for element in p:
            count = PatternCount.count(text[i:i+L], element)
            if count >= t and patterns.count(element) == 0:
                patterns.append(element)
    return patterns


file = open("EColi.txt", "r")
genome = file.readline()
l2 = file.readline().split(' ')
kp = int(l2[0])
Lp = int(l2[1])
tp = int(l2[2])
for elem in clump(genome, kp, Lp, tp):
    print(elem, end=" ")
