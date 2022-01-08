def matching (pattern, genome):
    pattLen = len(pattern)
    genLen = len(genome)
    output = list()

    for i in range(genLen - pattLen):
        test = genome[i:i+pattLen-1]
        if test == pattern[0:pattLen-1]:
            output.append(i)

    return output


file = open("Vibrio_cholerae.txt", "r")
pat = file.readline()
gen = file.readline()
final = matching(pat, gen)
for i in final:
    print(i, end=" ")
