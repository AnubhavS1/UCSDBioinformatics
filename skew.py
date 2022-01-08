# finds skew(Genome): #G - #C


def skew(genome):
    skews = list()
    skews.append(0)
    for i in range(len(genome)):
        if genome[i] == 'G':
            skews.append(skews[i]+1)
        elif genome[i] == 'C':
            skews.append(skews[i]-1)
        else:
            skews.append(skews[i])
    minVal = min(skews)
    minSkews = list()
    for i in range(len(skews)):
        if skews[i] == minVal:
            minSkews.append(i)
    return minSkews