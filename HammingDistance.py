# find number of mismatches between two strings
def HammingDistance(p, q):
    counter = 0
    for i in range(len(q)):
        if p[i] != q[i]:
            counter += 1
    return counter