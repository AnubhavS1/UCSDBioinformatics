import numpy as np
from collections import defaultdict


# get change based on denomination of coins
def dpchange(money, coins):
    min_coins = [0]
    for i in range(1, money + 1):
        min_coins.append(float('Inf'))
        for coin in coins:
            if i >= coin:
                if min_coins[i-coin] + 1 < min_coins[i]:
                    min_coins[i] = min_coins[i-coin] + 1
    return min_coins[money]


# find the length of the longest path in a weighted graph
def manhattan_tourist(n, m, down, right):
    s = np.zeros((n+1, m+1), dtype=int)
    for i in range(1, n+1):
        s[i][0] = s[i-1][0] + down[i-1][0]
    for i in range(1, m+1):
        s[0][i] = s[0][i-1] + right[0][i-1]
    for i in range(1, n+1):
        for j in range(1, m+1):
            s[i][j] = max(s[i-1][j] + down[i-1][j], s[i][j-1] + right[i][j-1])
    return s[n][m]


# construct a backtrack instructions to output LCS
def LCSbacktrack(v, w):
    # record score
    s = np.zeros((len(v)+1, len(w)+1), dtype=int)
    # record backtrack through the graph
    backtrack = np.zeros((len(v)+1, len(w)+1), dtype=int)
    # 1 - up, 2 - left, 3 - diagonal
    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            m = 0
            if v[i-1] == w[i-1]:
                m = 1
            s[i, j] = max(s[i - 1, j], s[i, j - 1], s[i - 1, j - 1] + m)
            # else:
            #     s[i, j] = max(s[i-1, j], s[i, j-1])
            if s[i, j] == s[i-1, j]:
                backtrack[i, j] = 1
            elif s[i, j] == s[i, j-1]:
                backtrack[i, j] = 2
            elif s[i, j] == s[i-1, j-1] + m:
                backtrack[i, j] = 3
    return backtrack


def outputLCS(backtrack, v, i, j):
    final = ''
    while i != 0 and j != 0:
        if backtrack[i, j] == 1:
            i -= 1
        elif backtrack[i, j] == 2:
            j -= 1
        else:
            final += v[i-1]
            i -= 1
            j -= 1
    return final[::-1]


# Calculate the alignment score of two peptides
indel = 5
blosum = {
    'A': {'A': 4, 'C': 0, 'E': -1, 'D': -2, 'G': 0, 'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1,
          'N': -2, 'Q': -1, 'P': -1, 'S': 1, 'R': -1, 'T': 0, 'W': -3, 'V': 0, 'Y': -2},
    'C': {'A': 0, 'C': 9, 'E': -4, 'D': -3, 'G': -3, 'F': -2, 'I': -1, 'H': -3, 'K': -3, 'M': -1, 'L': -1,
          'N': -3, 'Q': -3, 'P': -3, 'S': -1, 'R': -3, 'T': -1, 'W': -2, 'V': -1, 'Y': -2},
    'E': {'A': -1, 'C': -4, 'E': 5, 'D': 2, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 1, 'M': -2, 'L': -3,
          'N': 0, 'Q': 2, 'P': -1, 'S': 0, 'R': 0, 'T': -1, 'W': -3, 'V': -2, 'Y': -2},
    'D': {'A': -2, 'C': -3, 'E': 2, 'D': 6, 'G': -1, 'F': -3, 'I': -3, 'H': -1, 'K': -1, 'M': -3, 'L': -4,
          'N': 1, 'Q': 0, 'P': -1, 'S': 0, 'R': -2, 'T': -1, 'W': -4, 'V': -3, 'Y': -3},
    'G': {'A': 0, 'C': -3, 'E': -2, 'D': -1, 'G': 6, 'F': -3, 'I': -4, 'H': -2, 'K': -2, 'M': -3, 'L': -4,
          'N': 0, 'Q': -2, 'P': -2, 'S': 0, 'R': -2, 'T': -2, 'W': -2, 'V': -3, 'Y': -3},
    'F': {'A': -2, 'C': -2, 'E': -3, 'D': -3, 'G': -3, 'F': 6, 'I': 0, 'H': -1, 'K': -3, 'M': 0, 'L': 0,
          'N': -3, 'Q': -3, 'P': -4, 'S': -2, 'R': -3, 'T': -2, 'W': 1, 'V': -1, 'Y': 3},
    'I': {'A': -1, 'C': -1, 'E': -3, 'D': -3, 'G': -4, 'F': 0, 'I': 4, 'H': -3, 'K': -3, 'M': 1, 'L': 2,
          'N': -3, 'Q': -3, 'P': -3, 'S': -2, 'R': -3, 'T': -1, 'W': -3, 'V': 3, 'Y': -1},
    'H': {'A': -2, 'C': -3, 'E': 0, 'D': -1, 'G': -2, 'F': -1, 'I': -3, 'H': 8, 'K': -1, 'M': -2, 'L': -3,
          'N': 1, 'Q': 0, 'P': -2, 'S': -1, 'R': 0, 'T': -2, 'W': -2, 'V': -3, 'Y': 2},
    'K': {'A': -1, 'C': -3, 'E': 1, 'D': -1, 'G': -2, 'F': -3, 'I': -3, 'H': -1, 'K': 5, 'M': -1, 'L': -2,
          'N': 0, 'Q': 1, 'P': -1, 'S': 0, 'R': 2, 'T': -1, 'W': -3, 'V': -2, 'Y': -2},
    'M': {'A': -1, 'C': -1, 'E': -2, 'D': -3, 'G': -3, 'F': 0, 'I': 1, 'H': -2, 'K': -1, 'M': 5, 'L': 2,
          'N': -2, 'Q': 0, 'P': -2, 'S': -1, 'R': -1, 'T': -1, 'W': -1, 'V': 1, 'Y': -1},
    'L': {'A': -1, 'C': -1, 'E': -3, 'D': -4, 'G': -4, 'F': 0, 'I': 2, 'H': -3, 'K': -2, 'M': 2, 'L': 4,
          'N': -3, 'Q': -2, 'P': -3, 'S': -2, 'R': -2, 'T': -1, 'W': -2, 'V': 1, 'Y': -1},
    'N': {'A': -2, 'C': -3, 'E': 0, 'D': 1, 'G': 0, 'F': -3, 'I': -3, 'H': 1, 'K': 0, 'M': -2, 'L': -3,
          'N': 6, 'Q': 0, 'P': -2, 'S': 1, 'R': 0, 'T': 0, 'W': -4, 'V': -3, 'Y': -2},
    'Q': {'A': -1, 'C': -3, 'E': 2, 'D': 0, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 1, 'M': 0, 'L': -2,
          'N': 0, 'Q': 5, 'P': -1, 'S': 0, 'R': 1, 'T': -1, 'W': -2, 'V': -2, 'Y': -1},
    'P': {'A': -1, 'C': -3, 'E': -1, 'D': -1, 'G': -2, 'F': -4, 'I': -3, 'H': -2, 'K': -1, 'M': -2, 'L': -3,
          'N': -2, 'Q': -1, 'P': 7, 'S': -1, 'R': -2, 'T': -1, 'W': -4, 'V': -2, 'Y': -3},
    'S': {'A': 1, 'C': -1, 'E': 0, 'D': 0, 'G': 0, 'F': -2, 'I': -2, 'H': -1, 'K': 0, 'M': -1, 'L': -2,
          'N': 1, 'Q': 0, 'P': -1, 'S': 4, 'R': -1, 'T': 1, 'W': -3, 'V': -2, 'Y': -2},
    'R': {'A': -1, 'C': -3, 'E': 0, 'D': -2, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 2, 'M': -1, 'L': -2,
          'N': 0, 'Q': 1, 'P': -2, 'S': -1, 'R': 5, 'T': -1, 'W': -3, 'V': -3, 'Y': -2},
    'T': {'A': 0, 'C': -1, 'E': -1, 'D': -1, 'G': -2, 'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1,
          'N': 0, 'Q': -1, 'P': -1, 'S': 1, 'R': -1, 'T': 5, 'W': -2, 'V': 0, 'Y': -2},
    'W': {'A': -3, 'C': -2, 'E': -3, 'D': -4, 'G': -2, 'F': 1, 'I': -3, 'H': -2, 'K': -3, 'M': -1, 'L': -2,
          'N': -4, 'Q': -2, 'P': -4, 'S': -3, 'R': -3, 'T': -2, 'W': 11, 'V': -3, 'Y': 2},
    'V': {'A': 0, 'C': -1, 'E': -2, 'D': -3, 'G': -3, 'F': -1, 'I': 3, 'H': -3, 'K': -2, 'M': 1, 'L': 1,
          'N': -3, 'Q': -2, 'P': -2, 'S': -2, 'R': -3, 'T': 0, 'W': -3, 'V': 4, 'Y': -1},
    'Y': {'A': -2, 'C': -2, 'E': -2, 'D': -3, 'G': -3, 'F': 3, 'I': -1, 'H': 2, 'K': -2, 'M': -1, 'L': -1,
          'N': -2, 'Q': -1, 'P': -3, 'S': -2, 'R': -2, 'T': -2, 'W': 2, 'V': -1, 'Y': 7}}


# return score matrix, and a aligned sequences
def global_alignment(v, w):
    # construct score and backtrack matrix
    s = np.zeros((len(v) + 1, len(w) + 1), dtype=int)
    backtrack = np.zeros((len(v) + 1, len(w) + 1), dtype=int)
    # 1 - up, 2 - left, 3 - diagonal
    for i in range(1, len(v)+1):
        s[i, 0] = s[i-1, 0] - indel
        backtrack[i, 0] = 1
    for i in range(1, len(w)+1):
        s[0, i] = s[0, i-1] - indel
        backtrack[0, i] = 2
    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            up = s[i-1, j] - indel
            left = s[i, j-1] - indel
            diagonal = s[i-1, j-1] + blosum[v[i-1]][w[j-1]]
            s[i, j] = max(up, left, diagonal)
            if diagonal >= up:
                if diagonal >= left:
                    backtrack[i, j] = 3
                else:
                    backtrack[i, j] = 2
            else:
                if up >= left:
                    backtrack[i, j] = 1
                else:
                    backtrack[i, j] = 2
    # construct aligned sequences
    v_aligned, w_aligned = v, w
    insert_indel = lambda word, x: word[:x] + '-' + word[x:]
    n, m = len(v), len(w)
    while True:
        if backtrack[n, m] == 0:
            break
        if backtrack[n, m] == 1:
            n -= 1
            w_aligned = insert_indel(w_aligned, m)
        elif backtrack[n, m] == 2:
            m -= 1
            v_aligned = insert_indel(v_aligned, n)
        else:
            n -= 1
            m -= 1
    for r in range(n):
        w_aligned = insert_indel(w_aligned, 0)
    for r in range(m):
        v_aligned = insert_indel(v_aligned, 0)
    return str(s[len(v), len(w)]), v_aligned, w_aligned


# return max score and aligned sequences, taking into account gap penalties
def affineGap_alignment(v, w):
    # gap opening and extension penalties
    op = -11
    ep = -1

    n = len(v)
    m = len(w)

    s = np.zeros((3, n + 1, m + 1), dtype=int)
    backtrack = np.zeros((3, n + 1, m + 1), dtype=int)
    # fill in the sides of the score graph
    for i in range(1, n+1):
        s[0, i, 0] = op + ep * (i - 1)
        s[1, i, 0] = op + ep * (i - 1)
        s[2, i, 0] = 10 * op
    for i in range(1, m+1):
        s[0, 0, i] = 10 * op
        s[1, 0, i] = op + ep * (i - 1)
        s[2, 0, i] = op + ep * (i - 1)
    # fill in the rest of score and backtrack graphs
    for i in range(1, n+1):
        for j in range(1, m+1):
            bottom = [s[0, i-1, j] + ep, s[1, i-1, j] + op]
            s[0, i, j] = max(bottom)
            backtrack[0, i, j] = bottom.index(s[0, i, j])

            top = [s[2, i, j-1] + ep, s[1, i, j-1] + op]
            s[2, i, j] = max(top)
            backtrack[2, i, j] = top.index(s[2, i, j])

            mid = [s[0, i, j], s[1, i-1, j-1] + blosum[v[i-1]][w[j-1]], s[2, i, j]]
            s[1, i, j] = max(mid)
            backtrack[1, i, j] = mid.index(s[1, i, j])

    scores = [s[0, n, m], s[1, n, m], s[2, n, m]]
    matrix = scores.index(max(scores))
    # construct aligned sequences
    insert_indel = lambda string, i: string[:i] + '-' + string[i:]
    v_aligned, w_aligned = v, w
    while n * m != 0:
        if matrix == 0:
            if backtrack[0, n, m] == 1:
                matrix = 1
            n -= 1
            w_aligned = insert_indel(w_aligned, m)
        elif matrix == 1:
            if backtrack[1, n, m] == 0:
                matrix = 0
            elif backtrack[1, n, m] == 2:
                matrix = 2
            else:
                n -= 1
                m -= 1
        else:
            if backtrack[2, n, m] == 1:
                matrix = 1
            m -= 1
            v_aligned = insert_indel(v_aligned, n)
    for i in range(n):
        w_aligned = insert_indel(w_aligned, 0)
    for j in range(m):
        v_aligned = insert_indel(v_aligned, 0)

    return str(max(scores)), v_aligned, w_aligned


# return the middle edge of the matrix formed between two strings
def middle_edge(v, w):
    n = len(v)
    m = len(w)

    s = np.zeros((n + 1, m + 1), dtype=int)

    for i in range(1, n+1):
        s[i, 0] = -indel*i
    for i in range(1, m+1):
        s[0, i] = -indel*i
    for i in range(1, n+1):
        for j in range(1, int(m/2) + 1):
            tmp = s[i-1, j-1] + blosum[v[i-1]][w[j-1]]
            s[i, j] = max(s[i-1, j]-indel, s[i, j-1]-indel, tmp)

    mid_col1, mid_col2 = [], []
    for i in range(n+1):
        mid_col1.append(s[i, int(m/2)])
        mid_col2.append(s[i, int(m/2-1)])

    mid_node1 = max(mid_col1)
    mid_node2 = max(mid_col2)
    x1, x2 = 0, 0
    for i in range(n+1):
        if s[i, int(m/2)] == mid_node1:
            x1 = i
        if s[i, int(m/2-1)] == mid_node2:
            x2 = i
    return (x2+1, int(m/2)), (x1+1, int(m/2+1))


# return the length of the LCS of three strings and their alignments
def multipleLCS(u, v, w):
    l, n, m = len(u), len(v), len(w)

    s = np.zeros((l+1, n+1, m+1), dtype=int)
    backtrack = np.zeros((l+1, n+1, m+1), dtype=int)

    for i in range(1, l+1):
        for j in range(1, n+1):
            for k in range(1, m+1):
                if u[i-1] == v[j-1] == w[k-1]:
                    scores = [s[i-1, j, k], s[i, j-1, k], s[i, j, k-1], s[i-1, j-1, k], s[i-1, j, k-1], s[i, j-1, k-1], s[i-1, j-1, k-1]+1]
                    backtrack[i, j, k], s[i, j, k] = max(enumerate(scores, 1), key=lambda x: x[1])
                else:
                    scores = [s[i-1, j, k], s[i, j-1, k], s[i, j, k-1], s[i-1, j-1, k], s[i-1, j, k-1], s[i, j-1, k-1], s[i-1, j-1, k-1]]
                    backtrack[i, j, k], s[i, j, k] = max(enumerate(scores, 1), key=lambda x: x[1])

    insert_indel = lambda word, x: word[:x] + '-' + word[x:]
    aligned = [u, v, w]

    while l * n * m != 0:
        if backtrack[l, n, m] == 1:
            l -= 1
            aligned[1] = insert_indel(aligned[1], n)
            aligned[2] = insert_indel(aligned[2], m)
        elif backtrack[l, n, m] == 2:
            n -= 1
            aligned[0] = insert_indel(aligned[0], l)
            aligned[2] = insert_indel(aligned[2], m)
        elif backtrack[l, n, m] == 3:
            m -= 1
            aligned[0] = insert_indel(aligned[0], l)
            aligned[1] = insert_indel(aligned[1], n)
        elif backtrack[l,n, m] == 4:
            l -= 1
            n -= 1
            aligned[2] = insert_indel(aligned[2], m)
        elif backtrack[l, n, m] == 5:
            l -= 1
            m -= 1
            aligned[1] = insert_indel(aligned[1], n)
        elif backtrack[l, n, m] == 6:
            n -= 1
            m -= 1
            aligned[0] = insert_indel(aligned[0], l)
        else:
            l -= 1
            n -= 1
            m -= 1

    while len(aligned[0]) != max(len(aligned[0]), len(aligned[1]), len(aligned[2])):
        aligned[0] = insert_indel(aligned[0], 0)
    while len(aligned[1]) != max(len(aligned[0]), len(aligned[1]), len(aligned[2])):
        aligned[1] = insert_indel(aligned[1], 0)
    while len(aligned[2]) != max(len(aligned[0]), len(aligned[1]), len(aligned[2])):
        aligned[2] = insert_indel(aligned[2], 0)

    return str(s[len(u), len(v), len(w)]), aligned


# greedy permutation algorithm
def reverse(rev):
    tmp = rev[::-1]
    tmp = [-int(i) for i in tmp]
    for i in range(len(tmp)):
        t = tmp[i]
        if tmp[i] < 0:
            tmp[i] = str(t)
        else:
            tmp[i] = '+' + str(t)
    return tmp


def greedy_sort(p):
    permutations = []

    for pos in range(len(p)):
        if int(p[pos]) != pos:
            nums = [abs(int(i)) for i in p[pos:]]
            ind = nums.index(pos + 1) + pos
            test = p[pos:ind+1]
            if pos+1 != abs(int(p[-1])):
                p[pos:ind+1] = reverse(p[pos:ind+1])
            else:
                p[pos:] = reverse(p[pos:])
            permutations.append(p.copy())
            if int(p[pos]) < 0:
                tmp = abs(int(p[pos]))
                p[pos] = '+' + str(tmp)
                permutations.append(p.copy())
    return permutations


# count number of breakpoints in a permutation
def num_breakpoints(p):
    counter = 0
    for i in range(len(p)-1):
        if p[i] + 1 != p[i+1]:
            counter += 1
    return counter


# get reverse complement
def rev_comp(dna):
    d = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G'}
    return ''.join(d[base] for base in dna[::-1])


# find shared kmers between two dna strands
def shared_kmers(k, a, b):
    a_dict = defaultdict(list)
    for i in range(len(a) - k + 1):
        a_dict[a[i:i+k]].append(i)
    b_dict = defaultdict(list)
    for j in range(len(b) - k + 1):
        b_dict[b[j:j+k]].append(j)
    b_revcomp = defaultdict(list)
    revcomp = rev_comp(b)
    for j in range(len(b) - k + 1):
        b_revcomp[revcomp[j:j+k]].append(j)
    shared = []
    for key in a_dict.keys():
        if key in b_dict and (a_dict[key], b_dict[key]) not in shared:
            for num1 in a_dict[key]:
                for num2 in b_dict[key]:
                    shared.append((num1, num2))
        if key in b_revcomp and (a_dict[key], b_revcomp[key]) not in shared:
            for num1 in a_dict[key]:
                for num2 in b_revcomp[key]:
                    shared.append((num1, num2))
    return shared

