from scripts import *
from LongestPath import *
from TwoBreak import *
import numpy as np

file = open("test.txt", 'r')
output = open("out.txt", 'w')


def test_dpchange():
    money = int(file.readline().strip())
    coins = [int(x) for x in file.readline().split(",")]
    output.write(str(dpchange(money, coins)))


def test_manhattan():
    n, m = file.readline().strip().split()
    matrices = file.read().split('-')

    r = matrices[1].split('\n')[1:]
    right = []
    for i in r:
        right.append(i.split())
    right = np.array(right, dtype=int)

    d = matrices[0].split('\n')[:-1]
    down = []
    for i in d:
        down.append(i.split())
    down = np.array(down, dtype=int)

    output.write(str(manhattan_tourist(int(n), int(m), down, right)))


def test_outputLCS():
    v = file.readline().strip()
    w = file.readline()
    backtrack = LCSbacktrack(v, w)
    output.write(outputLCS(backtrack, v, len(v), len(w)))


def test_longestPath():
    content = [line.strip() for line in file.readlines()]
    start = int(content[0])
    end = int(content[1])
    graph = {}
    for elem in content[2:]:
        tmp = elem.split('->')
        key = int(tmp[0])
        val = (int(tmp[1].split(':')[0]), int(tmp[1].split(':')[1]))
        if key not in graph.keys():
            graph[key] = {val[0]: val[1]}
        else:
            graph[key][val[0]] = val[1]
    output.write(str(longest_path(graph, start, end)[0]) + '\n' + longest_path(graph, start, end)[1])


def test_galignment():
    v = file.readline().strip()
    w = file.readline().strip()
    output.write(global_alignment(v, w)[0] + '\n')
    output.write(global_alignment(v, w)[1] + '\n')
    output.write(global_alignment(v, w)[2])


def test_affineGap_alignment():
    v = file.readline().strip()
    w = file.readline().strip()
    a = affineGap_alignment(v, w)
    output.write(a[0] + '\n')
    output.write(a[1] + '\n')
    output.write(a[2])


def test_ME():
    v = file.readline().strip()
    w = file.readline().strip()
    coords = middle_edge(v, w)
    output.write(str(coords[0]) + ' ' + str(coords[1]))


def test_multipleLCS():
    u = file.readline().strip()
    v = file.readline().strip()
    w = file.readline().strip()
    lcs = multipleLCS(u, v, w)
    output.write(lcs[0] + '\n' + lcs[1][0] + '\n' + lcs[1][1] + '\n' + lcs[1][2])


def test_greedy():
    p = file.readline().strip().split()
    g = greedy_sort(p)
    # for elem in g:
    #     if g.count(elem) > 0:
    #         g.remove(g[g.index(elem)])
    for elem in g:
        output.write('(')
        for num in elem:
            if num != elem[len(elem)-1]:
                output.write(num + ' ')
            else:
                output.write(num)
        output.write(')')
        output.write('\n')

    seen = set()
    with open('out.txt') as f:
        for line in f:
            if line in seen:
                print('(' + line + ')')
            else:
                seen.add(line)


def test_breakpoints():
    p = [0]
    for num in file.readline().split():
        p.append(int(num))
    p.append(len(p))
    print(num_breakpoints(p))


def test_twoBreakDist():
    p, q = [line.strip().lstrip('(').rstrip(')').split(')(') for line in file.readlines()]
    P, Q = [], []
    for i in range(len(p)):
        P.append([])
        for num in p[i].split():
            P[i].append(int(num))
    for i in range(len(q)):
        Q.append([])
        for num in q[i].split():
            Q[i].append(int(num))
    print(two_break_dist(P, Q))


def test_graphtogenome():
    tmp = file.readline().split('), (')
    tmp[0] = tmp[0].lstrip('(')
    tmp[len(tmp)-1] = tmp[len(tmp)-1].rstrip(')')
    str = []
    for tup in tmp:
        for num in tup.split(', '):
            str.append(int(num))
    print(graphtogenome(str))
    # for block in graphtogenome(str):
    #     output.write('(')
    #     for num in block:
    #         if num > 0:
    #             output.write()


def test_sharedkmers():
    k = int(file.readline().strip())
    a = file.readline().strip()
    b = file.readline().strip()
    tmp = shared_kmers(k, a, b)
    shared = sorted(tmp, key=lambda x: x[1])
    for c in shared:
        output.write(str(c) + '\n')


# test_dpchange()
# test_manhattan()
# test_outputLCS()
# test_longestPath()
# test_galignment()
# test_affineGap_alignment()
# test_ME()
# test_multipleLCS()
# test_greedy()
# test_breakpoints()
# test_twoBreakDist()
# test_graphtogenome()
test_sharedkmers()
