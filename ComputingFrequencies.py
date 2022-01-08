def PatternToNumber(pattern):
    temp = ''
    for char in pattern:
        if char == 'A':
            temp += '0'
        elif char == 'C':
            temp += '1'
        elif char == 'G':
            temp += '2'
        elif char == 'T':
            temp += '3'
    return int(temp, 4)


def NumberToPattern(index, k):
    temp = list()
    for i in range(int(k)):
        temp.append(index % 4)
        index = int(index / 4)
    final = ''
    for i in range(len(temp)):
        if temp[len(temp)-1-i] == 0:
            final += 'A'
        elif temp[len(temp)-1-i] == 1:
            final += 'C'
        elif temp[len(temp)-1-i] == 2:
            final += 'G'
        elif temp[len(temp)-1-i] == 3:
            final += 'T'
    return final


def ComputingFrequencies(genome, k):
    freqs = list()
    x = 4**int(k)
    for i in range(x):
        freqs.append(0)
    for i in range(len(genome)-int(k) + 1):
        pattern = genome[i:i+int(k)]
        j = PatternToNumber(pattern)
        freqs[j] += 1
    return freqs

